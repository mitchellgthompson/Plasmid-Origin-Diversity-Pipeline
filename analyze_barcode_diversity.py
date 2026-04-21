"""Audit the diversity of the 20 bp RB-TnSeq barcodes in `final_order.csv`.

Computes metrics that together establish whether the barcodes are robust for
BarSeq demultiplexing:

  1. Pairwise Hamming distance distribution
     - Hamming distance d = number of mismatched positions in aligned 20-mers.
     - A barcode set with min(d) >= 3 can detect any 2-error read and correct
       any single-error read. min(d) >= 4 allows correction of 1-2 errors.

  2. Pairwise Levenshtein (edit) distance distribution
     - Accounts for indels as well as substitutions.
     - Important because RB-TnSeq can introduce indels during sequencing.

  3. Positional base composition
     - Each of the 20 positions should have ~25 % A/C/G/T so that Illumina
       base-calling is unbiased across the cycle (low-diversity positions
       cause phasing errors).

  4. GC-content distribution.

  5. Self / cross reverse-complement check (no barcode should equal the
     revcomp of any other barcode or of itself — would collapse under BarSeq).

  6. Shared k-mer overlap — fraction of 8-mers shared between any two
     barcodes. High overlap reduces demux robustness.

  7. Comparison against random 20-mer baseline at the same set size.

Writes `barcode_diversity_report.txt` and `barcode_diversity_report.pdf`
(pairwise-distance histograms + positional composition heatmap).
"""
import csv
import math
import random
import statistics
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "final_order"
OUT_DIR.mkdir(exist_ok=True)
IN_CSV = OUT_DIR / "final_order.csv"
TXT_OUT = OUT_DIR / "barcode_diversity_report.txt"
PDF_OUT = OUT_DIR / "barcode_diversity_report.pdf"

BARCODE_BP = 20


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def levenshtein(a, b):
    """Iterative DP; a, b are same-length short strings here (20 bp), so O(n^2) is fine."""
    n, m = len(a), len(b)
    prev = list(range(m + 1))
    curr = [0] * (m + 1)
    for i in range(1, n + 1):
        curr[0] = i
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev, curr = curr, prev
    return prev[m]


def random_baseline(n, bp=20, seed=1):
    rng = random.Random(seed)
    return ["".join(rng.choices("ACGT", k=bp)) for _ in range(n)]


def shannon_entropy(counter, total):
    h = 0.0
    for v in counter.values():
        p = v / total
        if p > 0:
            h -= p * math.log2(p)
    return h


def summarise_pairwise(values, label):
    mn = min(values); mx = max(values); md = statistics.median(values)
    mean = statistics.mean(values); sd = statistics.pstdev(values)
    return {
        "label": label, "min": mn, "max": mx, "median": md, "mean": mean, "sd": sd,
        "count": len(values),
    }


def main():
    rows = list(csv.DictReader(IN_CSV.open()))
    barcodes = [r["barcode_N20"] for r in rows]
    n = len(barcodes)
    assert all(len(b) == BARCODE_BP for b in barcodes), "Barcode length mismatch"
    assert len(set(barcodes)) == n, "Non-unique barcodes detected"

    # 1. Pairwise Hamming + Levenshtein
    pairs = list(combinations(range(n), 2))
    ham = [hamming(barcodes[i], barcodes[j]) for i, j in pairs]
    lev = [levenshtein(barcodes[i], barcodes[j]) for i, j in pairs]

    ham_stats = summarise_pairwise(ham, "pairwise Hamming distance (bp)")
    lev_stats = summarise_pairwise(lev, "pairwise Levenshtein distance (bp)")

    # 1b. Identify any "near-neighbour" pairs (lowest-distance pairs)
    ham_min = min(ham)
    low_pairs = [(i, j, hamming(barcodes[i], barcodes[j]))
                 for i, j in pairs if hamming(barcodes[i], barcodes[j]) == ham_min]

    # 2. Baseline: random 20-mers, same set size
    base = random_baseline(n, bp=BARCODE_BP, seed=1)
    base_ham = [hamming(base[i], base[j]) for i, j in pairs]

    # 3. Positional base composition
    pos_counts = [Counter() for _ in range(BARCODE_BP)]
    for b in barcodes:
        for i, ch in enumerate(b):
            pos_counts[i][ch] += 1
    pos_entropy = [shannon_entropy(c, n) for c in pos_counts]  # max = 2 bits for 4 bases
    pos_max_freq = [max(c.values()) / n for c in pos_counts]

    # 4. GC distribution
    gcs = [(b.count("G") + b.count("C")) / BARCODE_BP for b in barcodes]

    # 5. RevComp collisions
    rc_set = {revcomp(b) for b in barcodes}
    self_palindrome = [b for b in barcodes if revcomp(b) == b]
    rc_cross = [b for b in barcodes if revcomp(b) in set(barcodes) and revcomp(b) != b]

    # 6. Shared 8-mer overlap between barcode pairs
    K = 8
    kmers_per_bc = [set(b[i:i + K] for i in range(BARCODE_BP - K + 1)) for b in barcodes]
    shared_k = []
    for i, j in pairs:
        shared_k.append(len(kmers_per_bc[i] & kmers_per_bc[j]))
    shared_k_stats = summarise_pairwise(shared_k, f"shared {K}-mers between barcode pairs")

    # 7. Write text report
    lines = []
    push = lines.append
    push("=" * 72)
    push("RB-TnSeq BARCODE DIVERSITY AUDIT")
    push("=" * 72)
    push(f"Source:       {IN_CSV.name}")
    push(f"Barcodes:     n = {n}, length = {BARCODE_BP} bp")
    push(f"Uniqueness:   all {n} barcodes distinct")
    push("")
    push("--- Pairwise Hamming distance ---")
    push(f"  min:    {ham_stats['min']} bp   (target >= 3)")
    push(f"  median: {ham_stats['median']}")
    push(f"  mean:   {ham_stats['mean']:.2f}")
    push(f"  max:    {ham_stats['max']}")
    push(f"  sd:     {ham_stats['sd']:.2f}")
    push(f"  Hamming distance histogram (counts):")
    ham_hist = Counter(ham)
    for d in sorted(ham_hist):
        bar = "#" * min(60, ham_hist[d] // max(1, len(ham) // 300))
        push(f"    d={d:>2}:  {ham_hist[d]:>5}  {bar}")
    push("")
    push(f"  Random 20-mer baseline (n={n}): mean={statistics.mean(base_ham):.2f}, min={min(base_ham)}")
    push(f"  Empirical vs random: {'OK' if ham_stats['min'] >= 3 and ham_stats['mean'] >= statistics.mean(base_ham) - 1 else 'WARNING'}")
    push("")
    push(f"  Closest pairs at Hamming = {ham_min} bp ({len(low_pairs)} pair(s)):")
    for i, j, d in low_pairs[:10]:
        push(f"    {barcodes[i]}  <-> {barcodes[j]}   (origins: {rows[i]['id']}  <-> {rows[j]['id']})")
    push("")
    push("--- Pairwise Levenshtein (edit) distance ---")
    push(f"  min:    {lev_stats['min']}   median: {lev_stats['median']}   mean: {lev_stats['mean']:.2f}   max: {lev_stats['max']}")
    push("  Levenshtein <= Hamming; a Levenshtein-min of 3 means no pair can be reached by fewer than 3 edits.")
    push("")
    push("--- Positional base composition (% A / C / G / T per column) ---")
    push("  pos  A%    C%    G%    T%     max    entropy (bits, max 2.00)")
    for i, c in enumerate(pos_counts):
        a = c.get('A',0)/n*100; cc = c.get('C',0)/n*100
        g = c.get('G',0)/n*100; t = c.get('T',0)/n*100
        push(f"  {i+1:>3}  {a:>4.1f}  {cc:>4.1f}  {g:>4.1f}  {t:>4.1f}   {pos_max_freq[i]*100:>4.1f}   {pos_entropy[i]:.3f}")
    push(f"  Median positional entropy: {statistics.median(pos_entropy):.3f}")
    push(f"  Median positional max-base frequency: {statistics.median(pos_max_freq)*100:.1f} %")
    push(f"  Guidance: position entropy should be close to 2.00 bits; max-base freq << 40 % avoids Illumina phasing bias.")
    push("")
    push("--- GC content ---")
    push(f"  range: {min(gcs)*100:.0f} % - {max(gcs)*100:.0f} %    mean: {statistics.mean(gcs)*100:.1f} %")
    push(f"  All within design constraint [40 %, 60 %]: {all(0.40 <= g <= 0.60 for g in gcs)}")
    push("")
    push("--- Reverse-complement collisions ---")
    push(f"  Self-palindrome barcodes: {len(self_palindrome)}")
    push(f"  Cross-revcomp pairs:      {len(rc_cross)//2}  (would collapse under BarSeq if > 0)")
    push("")
    push(f"--- Shared {K}-mer overlap between barcode pairs ---")
    push(f"  min: {shared_k_stats['min']}   median: {shared_k_stats['median']}   mean: {shared_k_stats['mean']:.2f}   max: {shared_k_stats['max']}")
    push(f"  Pairs sharing >= 5 consecutive 8-mers would indicate low diversity; observed max = {shared_k_stats['max']}.")
    push("")
    push("--- Verdict ---")
    pass_checks = []
    pass_checks.append(("min Hamming >= 3", ham_stats['min'] >= 3))
    pass_checks.append(("min Levenshtein >= 3", lev_stats['min'] >= 3))
    pass_checks.append(("all barcodes unique", len(set(barcodes)) == n))
    pass_checks.append(("no revcomp collisions", len(self_palindrome) == 0 and len(rc_cross) == 0))
    pass_checks.append(("positional entropy >= 1.7 bits (4-base balance)", all(e >= 1.7 for e in pos_entropy)))
    pass_checks.append(("positional max-base freq < 40 %", all(f < 0.40 for f in pos_max_freq)))
    pass_checks.append(("GC in [40, 60] %", all(0.40 <= g <= 0.60 for g in gcs)))
    for label, ok in pass_checks:
        push(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    all_pass = all(ok for _, ok in pass_checks)
    push("")
    push(f"OVERALL: {'PASS - barcodes are significantly dissimilar' if all_pass else 'FAIL - see checks above'}")
    TXT_OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {TXT_OUT.name}")

    # 8. PDF: histograms + positional composition heatmap
    with PdfPages(str(PDF_OUT)) as pdf:
        fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))

        # Hamming histogram
        ax = axs[0, 0]
        bins = range(0, BARCODE_BP + 2)
        ax.hist(ham, bins=bins, color="#1565C0", alpha=0.85, label="Empirical")
        ax.hist(base_ham, bins=bins, color="#B0B0B0", alpha=0.55, label="Random baseline (n=73)")
        ax.axvline(3, color="red", linestyle="--", lw=1, label="min target (>=3)")
        ax.set_xlabel("Hamming distance (bp)")
        ax.set_ylabel("# barcode pairs")
        ax.set_title("Pairwise Hamming distance")
        ax.legend(fontsize=8)

        # Levenshtein histogram
        ax = axs[0, 1]
        ax.hist(lev, bins=bins, color="#6A1B9A", alpha=0.85)
        ax.axvline(3, color="red", linestyle="--", lw=1)
        ax.set_xlabel("Levenshtein distance")
        ax.set_ylabel("# barcode pairs")
        ax.set_title("Pairwise Levenshtein (edit) distance")

        # Positional composition heatmap
        ax = axs[1, 0]
        bases = ["A", "C", "G", "T"]
        grid = np.array([[c.get(b, 0) / n for b in bases] for c in pos_counts])
        im = ax.imshow(grid.T, aspect="auto", cmap="Blues", vmin=0.0, vmax=0.5)
        ax.set_xticks(range(BARCODE_BP)); ax.set_xticklabels(range(1, BARCODE_BP + 1), fontsize=7)
        ax.set_yticks(range(4)); ax.set_yticklabels(bases)
        ax.set_xlabel("barcode position"); ax.set_title("Positional base frequency")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="freq")

        # GC distribution
        ax = axs[1, 1]
        ax.hist(gcs, bins=np.linspace(0.3, 0.7, 17), color="#2E7D32", alpha=0.85)
        ax.axvline(0.40, color="red", linestyle="--", lw=1)
        ax.axvline(0.60, color="red", linestyle="--", lw=1)
        ax.set_xlabel("GC fraction"); ax.set_ylabel("# barcodes")
        ax.set_title("GC-content distribution")

        fig.suptitle(f"RB-TnSeq barcode diversity audit — n={n}", fontsize=13, weight="bold")
        plt.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote {PDF_OUT.name}")


if __name__ == "__main__":
    main()
