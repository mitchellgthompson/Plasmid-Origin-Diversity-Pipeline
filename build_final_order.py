"""Build the final 250 kb synthesis/cloning order with RB-TnSeq barcode cassettes.

Starting from the combined reference CSV, applies:
  1. Ian Blaby constraint-fail list (Twist review, 2026-04-20)
  2. Exact-sequence deduplication
  3. Rep-family exclusions (Gemini, MobT, Phg_2220_C)
  4. Tier-priority greedy fit (Tier 3 > 2 > 4 > 5) to a 250,000 bp TOTAL order
     budget that accounts for the 56 bp barcode cassette added to every origin
  5. A deterministic 20 bp barcode per origin, flanked by the canonical
     Deutschbauer-lab RB-TnSeq priming sites (Wetmore et al. 2015, mBio)

Cassette structure, top strand 5'->3', appended to the 3' end of each origin:
    A (1 bp spacer) - U1 (18 bp) - N20 (20 bp) - U2 (18 bp) = 57 bp total
    U1 = GATGTCCACGAGGTCTCT   (BarSeq_P2 primer-binding site)
    U2 = CGTACGCTGCAGGTCGAC   (BarSeq_P1 primer-binding site, top strand)

The 1 bp 'A' spacer sits between the origin 3' end and U1 to prevent a G
from U1 extending any G-homopolymer at the origin tail into a >=6 bp run
(which would fail Twist's G/C homopolymer constraint).  BarSeq amplification
is unaffected because the BarSeq_P2 primer's 3' gene-specific tail lands
entirely on U1 regardless of what sits 5' of U1.

Barcode design:
  - 20 bp, GC 40-60%
  - No homopolymer run of 4 or more
  - Pairwise Hamming distance >= 3
  - No U1/U2/revcomp substring of length >= 8
  - Not found (forward or revcomp) in the target origin sequence
  - Deterministic RNG seed for reproducibility

Produces `final_order.csv` with columns added:
  - production_method  (twist_synthesis | genomic_pcr_amplification)
  - barcode_N20        (20 bp barcode)
  - barcoded_sequence  (origin + U1 + N20 + U2 — what Twist synthesizes)

Pools:
  - repABC_primary       -> twist_synthesis      (26 origins)
  - repABC_alt_cloning   -> genomic_pcr_amplification (7 pTi-family origins
                           that failed Twist constraints; amplified from
                           genomic/plasmid template)
  - diversity_PLSDB      -> twist_synthesis      (41 origins)
"""
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
IN_CSV = ROOT / "combined_synthesis_origins.csv"
OUT_DIR = ROOT / "final_order"
OUT_DIR.mkdir(exist_ok=True)
OUT_CSV = OUT_DIR / "final_order.csv"
OUT_FASTA = OUT_DIR / "final_order_barcoded.fasta"

BUDGET_BP = 250_000
RNG_SEED = 20260420  # date of Ian's review; frozen for reproducibility

U1 = "GATGTCCACGAGGTCTCT"  # BarSeq_P2 priming site (upstream of barcode)
U2 = "CGTACGCTGCAGGTCGAC"  # BarSeq_P1 priming site (downstream of barcode)
SPACER = "A"                # 1 bp spacer between origin and U1
BARCODE_BP = 20
CASSETTE_TOTAL_BP = len(SPACER) + len(U1) + BARCODE_BP + len(U2)  # 57

TWIST_CLONAL_GENE_MAX_BP = 5_000  # origins exceeding this after cassette are flagged

REPABC_ALT_CLONING = {
    "C58_AE007871.2_pTi_Type_I.a_repABC",
    "1D1524_plasmid2_pTi_Type_IV.a_repABC",
    "AMS008_plasmid3_pTi_Type_IV.b_repABC",
    "CG474_CG474_3_pTi_Type_I.c_repABC",
    "Q15_94_contig_5_pTi_Type_VI_repABC",
    "CA75_95_contig_3_pTi_Type_IV.c_repABC",
    "Atu_K84_pAtK84b_CP000630.1_repABC",
}

# Origins dropped because the NCBI assembly contains ambiguous bases (N) that
# cannot be resolved by homology — Twist rejects redundant bases. The greedy
# fit will backfill with the next-best candidate from the same tier system.
AMBIGUOUS_BASE_DROPS = {
    "1908110475",  # Defluviicoccus vanus (NZ_CP053924.1) — 6 N's in RepC CDS
                   # at codons 182, 291, 300; no clear blastp consensus.
}

# Origins dropped because their Rep family lacks established bacterial-host
# precedence — keeps the library to Rep families known to function in
# bacterial hosts. Greedy will backfill the freed budget with the next-best
# bacterial candidate.
NO_BACTERIAL_PRECEDENCE_DROPS = {
    "2270922916",  # Haloarcula marina (RNA_helicase__RCR, archaeal halophile)
}

DIVERSITY_FAILS = {
    "1884895362","1928226770","1008893761","2708435457","1608334067","187736698",
    "1054801980","2428992113","636800953","2218211345","1848374391","870972790",
    "2285137328","1908983349","1635401334","169546397","190410555","2486783363",
    "2270919113","1148871448","1888256587","2664483544","2674559851","2674557724",
    "1724206600","2711208204","1950532025","1869459527","2623834384","209809836",
    "1883370430","2452734037","2272232194",
}

EXCLUDED_FAMILIES = {"MobT", "Phg_2220_C", "Gemini_AL1__RCR"}

TIER_RANK = {
    "Tier 3: Cross-phylum BHR": 0,
    "Tier 2: Alphaproteobacteria": 1,
    "Tier 4: Phylum-exclusive": 2,
    "Tier 5: Diversity fill": 3,
}


def revcomp(s):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(table)[::-1]


def gc_frac(s):
    return (s.count("G") + s.count("C")) / len(s)


def has_homopolymer(s, n=4):
    for b in "ACGT":
        if b * n in s:
            return True
    return False


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


U_KMERS = set()
for flank in (U1, U2, revcomp(U1), revcomp(U2)):
    for i in range(len(flank) - 7):
        U_KMERS.add(flank[i:i + 8])


def contains_flank_kmer(s, k=8):
    for i in range(len(s) - k + 1):
        if s[i:i + k] in U_KMERS:
            return True
    return False


def generate_barcodes(n, origin_seqs, min_hd=3, seed=RNG_SEED):
    """Generate `n` unique 20 bp barcodes satisfying design constraints.

    Constraints:
      - GC fraction in [0.40, 0.60]
      - No homopolymer of length >=4
      - No 8-mer overlap with U1/U2 (forward or revcomp)
      - Does not appear (forward or revcomp) in any origin sequence
      - Pairwise Hamming distance >= min_hd
    """
    rng = random.Random(seed)
    pool = []
    attempts = 0
    # Uppercase origin sequences once for searching
    origin_seqs_upper = [s.upper() for s in origin_seqs]
    while len(pool) < n:
        attempts += 1
        if attempts > 2_000_000:
            raise RuntimeError(f"Barcode generation stalled after {attempts} attempts ({len(pool)}/{n})")
        bc = "".join(rng.choices("ACGT", k=BARCODE_BP))
        if not 0.40 <= gc_frac(bc) <= 0.60:
            continue
        if has_homopolymer(bc, 4):
            continue
        if contains_flank_kmer(bc, 8):
            continue
        bc_rc = revcomp(bc)
        if any(bc in s or bc_rc in s for s in origin_seqs_upper):
            continue
        if any(hamming(bc, other) < min_hd for other in pool):
            continue
        # avoid near-palindrome barcodes that could fold on themselves
        if hamming(bc, bc_rc) < min_hd:
            continue
        pool.append(bc)
    return pool


def score_key(r):
    try:
        return -float(r["functional_or_synthesis_score"] or 0)
    except ValueError:
        return 0.0


def main():
    rows = list(csv.DictReader(IN_CSV.open()))

    # Split repABC primary into Twist vs alt-cloning pools
    repabc_twist, repabc_alt = [], []
    for r in rows:
        if r["source_tier"] != "repABC_primary":
            continue
        if r["id"] in REPABC_ALT_CLONING:
            repabc_alt.append(dict(r, source_tier="repABC_alt_cloning",
                                   production_method="genomic_pcr_amplification"))
        else:
            repabc_twist.append(dict(r, production_method="twist_synthesis"))

    # Diversity pool: drop fails, drop excluded families, dedup exact sequence
    div_pool = [r for r in rows if r["source_tier"] == "diversity_PLSDB"]
    div_pool = [r for r in div_pool if r["id"] not in DIVERSITY_FAILS]
    div_pool = [r for r in div_pool if r["id"] not in AMBIGUOUS_BASE_DROPS]
    div_pool = [r for r in div_pool if r["id"] not in NO_BACTERIAL_PRECEDENCE_DROPS]
    div_pool = [r for r in div_pool if r["rep_class_or_family"] not in EXCLUDED_FAMILIES]
    by_seq = defaultdict(list)
    for r in div_pool:
        by_seq[r["sequence"]].append(r)
    div_pool = [bucket[0] for bucket in by_seq.values()]

    # Greedy priority fit using EFFECTIVE size = origin_len + cassette_total
    # so that total order bp (origins + all cassettes) <= BUDGET_BP.
    fixed_order_bp = sum(len(r["sequence"]) + CASSETTE_TOTAL_BP for r in repabc_twist) \
                   + sum(len(r["sequence"]) + CASSETTE_TOTAL_BP for r in repabc_alt)
    div_budget = BUDGET_BP - fixed_order_bp

    ordered = sorted(div_pool, key=lambda r: (TIER_RANK.get(r["priority_or_category"], 9),
                                              score_key(r), len(r["sequence"])))
    picked, order_bp = [], 0
    for r in ordered:
        eff = len(r["sequence"]) + CASSETTE_TOTAL_BP
        if order_bp + eff <= div_budget:
            picked.append(dict(r, production_method="twist_synthesis"))
            order_bp += eff

    # Single smallest-remaining top-up: allow a modest overshoot (up to 1 kb over
    # the 250 kb total budget) so the library lands at the intended ~73-origin
    # size even when constraint-fail or ambiguous-base drops free a too-narrow gap.
    OVERSHOOT_TOLERANCE_BP = 1_500
    picked_ids = {p["id"] for p in picked}
    remaining = sorted([r for r in ordered if r["id"] not in picked_ids],
                       key=lambda r: len(r["sequence"]))
    for cand in remaining:
        eff = len(cand["sequence"]) + CASSETTE_TOTAL_BP
        if order_bp + eff <= div_budget + OVERSHOOT_TOLERANCE_BP:
            picked.append(dict(cand, production_method="twist_synthesis"))
            order_bp += eff
            break

    # Pre-barcode collision check: make sure no origin already contains U1/U2/revcomps.
    # If so, flag the origin — BarSeq would double-prime.
    flank_variants = [U1, U2, revcomp(U1), revcomp(U2)]
    flagged = []
    final_rows = repabc_twist + repabc_alt + picked
    for r in final_rows:
        s = r["sequence"].upper()
        hits = [v for v in flank_variants if v in s]
        if hits:
            flagged.append((r["id"], hits))

    # Generate barcodes
    origin_seqs = [r["sequence"] for r in final_rows]
    barcodes = generate_barcodes(len(final_rows), origin_seqs)

    # Append cassette to each origin; flag oversize + junction issues
    for r, bc in zip(final_rows, barcodes):
        cassette = SPACER + U1 + bc + U2
        origin = r["sequence"].upper()
        r["barcode_N20"] = bc
        r["barcoded_sequence"] = origin + cassette

        r["oversize_flag"] = "YES" if len(r["barcoded_sequence"]) > TWIST_CLONAL_GENE_MAX_BP else ""

        # Only flag if the cassette extends an origin-terminal homopolymer
        # across the junction. Pre-existing internal homopolymers in the origin
        # were already reviewed by Ian — don't re-flag them.
        tail_base = origin[-1]
        origin_tail_run = 1
        for i in range(len(origin) - 2, -1, -1):
            if origin[i] == tail_base:
                origin_tail_run += 1
            else:
                break
        cassette_head_run = 1 if cassette[0] == tail_base else 0
        for i in range(1, len(cassette)):
            if cassette[i] == tail_base and cassette_head_run == i:
                cassette_head_run += 1
            else:
                break
        crossing_run = origin_tail_run + cassette_head_run if cassette_head_run > 0 else 0
        if crossing_run == 0:
            r["junction_homopolymer_flag"] = ""
        elif (tail_base in "GC" and crossing_run >= 6) or (tail_base in "AT" and crossing_run >= 10):
            r["junction_homopolymer_flag"] = f"{crossing_run}x{tail_base} (extended by cassette)"
        else:
            r["junction_homopolymer_flag"] = ""

    # Write CSV
    out_cols = ["source_tier", "production_method", "id", "species_or_strain",
                "length_bp", "GC_content", "rep_class_or_family",
                "functional_or_synthesis_score", "priority_or_category",
                "est_cost_usd", "ease_score", "taxonomy_class", "taxonomy_family",
                "notes", "sequence", "barcode_N20", "barcoded_sequence",
                "oversize_flag", "junction_homopolymer_flag"]
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols)
        w.writeheader()
        w.writerows(final_rows)

    # Emit a FASTA of the barcoded sequences (what Twist/alt-clone will produce)
    with OUT_FASTA.open("w") as fh:
        for r in final_rows:
            flags = []
            if r["oversize_flag"]:
                flags.append("oversize")
            if r["junction_homopolymer_flag"]:
                flags.append(f"junction_{r['junction_homopolymer_flag'].split()[0]}")
            flag_str = f" flags={','.join(flags)}" if flags else ""
            header = (f">{r['id']} method={r['production_method']} "
                      f"barcode={r['barcode_N20']} origin_len={len(r['sequence'])}bp "
                      f"total_len={len(r['barcoded_sequence'])}bp{flag_str}")
            fh.write(header + "\n")
            seq = r["barcoded_sequence"]
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    # --- Report ---
    by_method = Counter(r["production_method"] for r in final_rows)
    origin_bp = sum(len(r["sequence"]) for r in final_rows)
    cassette_bp = CASSETTE_TOTAL_BP * len(final_rows)
    total_order_bp = sum(len(r["barcoded_sequence"]) for r in final_rows)

    print(f"Wrote {OUT_CSV.name}")
    print(f"  {len(final_rows)} origins: {by_method.get('twist_synthesis',0)} twist_synthesis + "
          f"{by_method.get('genomic_pcr_amplification',0)} genomic_pcr_amplification")
    print(f"  origin bp:    {origin_bp:>7,}")
    print(f"  cassette bp:  {cassette_bp:>7,}  ({len(final_rows)} x {CASSETTE_TOTAL_BP} bp)")
    print(f"  total order:  {total_order_bp:>7,}  (budget {BUDGET_BP:,}, Δ{total_order_bp-BUDGET_BP:+d})")
    print()
    print("  Pool breakdown:")
    for label, pool in [("repABC primary (Twist)", repabc_twist),
                         ("repABC alt_cloning (genomic PCR)", repabc_alt),
                         ("diversity (Twist)", picked)]:
        obp = sum(len(r["sequence"]) for r in pool)
        ord_bp = obp + CASSETTE_TOTAL_BP * len(pool)
        print(f"    {label:<35} n={len(pool):>3}  origin={obp:>7,} bp  order={ord_bp:>7,} bp")
    print()
    print("  Diversity tier mix:")
    for t, n in Counter(p["priority_or_category"] for p in picked).most_common():
        print(f"    {t}: {n}")
    print()
    print("  Barcode constraints satisfied:")
    print(f"    n={len(barcodes)}  len={BARCODE_BP}  min pairwise Hamming={min(hamming(a,b) for i,a in enumerate(barcodes) for b in barcodes[i+1:]):.0f}")
    print(f"    GC range: {min(gc_frac(b) for b in barcodes):.2f} - {max(gc_frac(b) for b in barcodes):.2f}")
    print()
    if flagged:
        print(f"  WARNING: {len(flagged)} origin(s) already contain a U1/U2 priming site — BarSeq would double-prime:")
        for oid, hits in flagged:
            print(f"    {oid}: {hits}")
    else:
        print("  No origin contains U1/U2 or their revcomps — no double-priming risk.")

    oversize = [r["id"] for r in final_rows if r["oversize_flag"]]
    homopoly = [(r["id"], r["junction_homopolymer_flag"]) for r in final_rows if r["junction_homopolymer_flag"]]
    print()
    print(f"  oversize_flag (>{TWIST_CLONAL_GENE_MAX_BP:,} bp after cassette): {len(oversize)}")
    for oid in oversize:
        bar_len = next(len(r["barcoded_sequence"]) for r in final_rows if r["id"] == oid)
        print(f"    {oid}: {bar_len:,} bp")
    print(f"  junction_homopolymer_flag (>=6 G/C or >=10 A/T at junction): {len(homopoly)}")
    for oid, h in homopoly:
        print(f"    {oid}: {h}")


if __name__ == "__main__":
    main()
