#!/usr/bin/env python3
"""
Convert pipeline origin FASTA sequences to individual annotated GenBank files.

After ORF prediction and functional classification, any origin where a
functionally important ORF (rep, partition, mobilization, recombinase)
sits within 200 bp of the window edge is EXTENDED using the full plasmid
sequence from PLSDB, then re-annotated.

Reads:
  - origin_sequences.fasta          (from pipeline output)
  - diversity_library_metadata.csv  (from pipeline output)
  - RIP HMM profile                 (for ORF classification)
  - PLSDB sequences.fasta           (full plasmid seqs for window extension)

Writes:
  - genbanks/<Accession>.gb   (one per origin, with annotated ORFs)

Usage:
  python origins_to_genbank.py \\
      --fasta      results_plsdb/origin_sequences.fasta \\
      --csv        results_plsdb/diversity_library_metadata.csv \\
      --rip-hmm    RIPs/RIP.hmm \\
      --plsdb-fasta plsdb2025/sequences.fasta \\
      --outdir     results_plsdb/genbanks
"""

import argparse
import os
import shutil
import subprocess
import tempfile
from datetime import date
from pathlib import Path

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# ── Constants ────────────────────────────────────────────────────────────────
BUFFER_BP = 200

# ── CDS classification keywords (mirrored from pipeline) ─────────────────────
REP_KW = {"replic", "initiator", "trfa", "repa", "repb", "repc",
           "helix-turn-helix", "dna-binding", "rep protein", "plasmid replication"}
PAR_KW = {"partiti", "parb", "para", "spo0j", "sopa", "sopb",
           "segregation", "plasmid stabiliz", "noc", "parg", "ribbon-helix"}
MOB_KW = {"mobili", "relax", "transfer", "conjugal", "trai", "traa",
           "trbb", "trbc", "mob protein", "type iv"}
REC_KW = {"transposase", "integrase", "recombinas", "resolvase",
           "is element", "insertion sequence"}
AMR_KW = {"resistance", "beta-lactam", "aminoglycoside", "tetracycline",
           "sulfonamide", "chloramphenicol", "fluoroquinolone",
           "macrolide", "rifampin", "colistin", "vancomycin"}

# Categories considered "functionally important" for buffer enforcement
IMPORTANT_CATS = {"rep", "partition", "mobilization", "recombinase"}

CATEGORY_COLORS = {
    "rep":          "#E53935",
    "partition":    "#1565C0",
    "mobilization": "#2E7D32",
    "recombinase":  "#7B1FA2",
    "amr":          "#FF6F00",
    "structural":   "#FF8F00",
    "hypothetical": "#9E9E9E",
    "other":        "#78909C",
}


def cds_category(product: str) -> str:
    pl = product.lower()
    if any(k in pl for k in REP_KW):    return "rep"
    if any(k in pl for k in PAR_KW):    return "partition"
    if any(k in pl for k in MOB_KW):    return "mobilization"
    if any(k in pl for k in REC_KW):    return "recombinase"
    if any(k in pl for k in AMR_KW):    return "amr"
    if "hypothetical" in pl or "unknown" in pl: return "hypothetical"
    return "other"


def run_hmmsearch(proteins_fasta: str, hmm_path: str,
                  evalue: float = 1e-5) -> dict[str, list[dict]]:
    if not shutil.which("hmmsearch"):
        print("WARNING: hmmsearch not on PATH -- ORFs will not be classified by HMM")
        return {}
    tbl = proteins_fasta + ".tbl"
    cmd = ["hmmsearch", "--tblout", tbl, "--noali",
           "-E", str(evalue), hmm_path, proteins_fasta]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    hits: dict[str, list[dict]] = {}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            orf_id   = parts[0]
            hmm_name = parts[2]
            evalue_v = float(parts[4])
            score_v  = float(parts[5])
            acc = "_".join(orf_id.split("_")[:-1])
            hits.setdefault(acc, []).append({
                "orf_id": orf_id, "name": hmm_name,
                "evalue": evalue_v, "score": score_v,
            })
    os.unlink(tbl)
    return hits


def predict_orfs(seq: str, acc: str) -> list[dict]:
    """Predict ORFs with pyrodigal for a single sequence."""
    finder = pyrodigal.GeneFinder(meta=True)
    genes = finder.find_genes(seq.encode())
    orfs = []
    for k, g in enumerate(genes):
        aa = g.translate().rstrip("*")
        if len(aa) >= 30:
            orfs.append({
                "id": f"{acc}_orf{k+1}",
                "start": g.begin, "end": g.end,
                "strand": 1 if g.strand == 1 else -1,
                "aa": aa,
            })
    return orfs


def classify_orfs_rip(all_orfs: dict[str, list[dict]], hmm_path: str) -> dict[str, str]:
    """Run hmmsearch against RIP.hmm, return {orf_id: hmm_name}."""
    orf_hmm: dict[str, str] = {}
    if not hmm_path or not os.path.isfile(hmm_path):
        return orf_hmm
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                      delete=False, prefix="gb_orfs_") as tmp:
        prot_path = tmp.name
        for acc, orfs in all_orfs.items():
            for orf in orfs:
                tmp.write(f">{orf['id']}\n{orf['aa']}\n")
    raw = run_hmmsearch(prot_path, hmm_path)
    os.unlink(prot_path)
    for acc_hits in raw.values():
        for hit in acc_hits:
            orf_hmm[hit["orf_id"]] = hit["name"]
    return orf_hmm


def classify_orfs_pfam(all_orfs: dict[str, list[dict]],
                        pfam_hmm: str) -> dict[str, tuple[str, str, str]]:
    """
    Run hmmscan against Pfam-A.hmm.
    Returns {orf_id: (pfam_name, pfam_acc, description)}.
    """
    if not pfam_hmm or not os.path.isfile(pfam_hmm):
        return {}
    if not shutil.which("hmmscan"):
        print("WARNING: hmmscan not on PATH")
        return {}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                      delete=False, prefix="pfam_orfs_") as tmp:
        prot_path = tmp.name
        for acc, orfs in all_orfs.items():
            for orf in orfs:
                tmp.write(f">{orf['id']}\n{orf['aa']}\n")

    domtbl = prot_path + ".domtbl"
    cmd = ["hmmscan", "--domtblout", domtbl, "--noali",
           "-E", "1e-3", "--cpu", "4", pfam_hmm, prot_path]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    os.unlink(prot_path)

    hits: dict[str, tuple[str, str, str]] = {}
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 23:
                continue
            pfam_name = p[0]
            pfam_acc  = p[1]
            orf_id    = p[3]
            evalue    = float(p[6])
            desc      = " ".join(p[22:])
            if orf_id not in hits or evalue < hits[orf_id][3]:
                hits[orf_id] = (pfam_name, pfam_acc, desc, evalue)
    os.unlink(domtbl)

    return {k: (v[0], v[1], v[2]) for k, v in hits.items()}


# Pfam families that are origin-relevant (replication, partitioning, regulation)
ORIGIN_PFAM = {
    # Replication initiators
    "Rep_trans", "Rep3_N", "RPA", "RepL", "RepSA", "Rep_1", "Rep_2", "Rep_3",
    "IncFII_repA", "Gemini_AL1", "DnaA_N", "AAA_31", "AAA_25", "RP-C",
    "IncW_repA", "Replicase", "DnaB_2", "Viral_rep", "Phage_rep_O",
    "DNA_primase_S", "Prim-Pol", "RCR", "Viral_Rep", "ColE2_RepA",
    # Partitioning
    "ParBc", "SopB_HTH", "ParB", "CbiA", "SpoVG",
    "Plasmid_parti", "HTH_Tnp_1",
    # Origin-binding / regulation
    "HTH_3", "HTH_13", "HTH_31", "HTH_XRE", "KilA-N", "IclR_C",
    "Peptidase_S24",  # LexA-like repressors often in origin regions
    # Stability / toxin-antitoxin (often in origin regions)
    "RelE", "VapC", "MazF", "CcdB", "Toxin_higB", "ParE_toxin",
    "AbrB", "MazE", "CcdA", "Antitoxin",
}


def pfam_is_origin_relevant(pfam_name: str) -> bool:
    """Check if a Pfam family is plausibly origin-relevant."""
    if pfam_name in ORIGIN_PFAM:
        return True
    pn = pfam_name.lower()
    return any(k in pn for k in {"rep", "part", "segr", "initi", "prim",
                                  "hth", "helix_turn", "plasmid", "kila",
                                  "antitoxin", "toxin", "stab"})


def label_orfs(all_orfs: dict[str, list[dict]],
               rip_hmm: dict[str, str],
               pfam_hits: dict[str, tuple[str, str, str]]):
    """Assign product name and category to each ORF using RIP HMM + Pfam-A."""
    for acc, orfs in all_orfs.items():
        for orf in orfs:
            oid = orf["id"]

            # Priority 1: RIP HMM hit → Rep protein
            if oid in rip_hmm:
                hmm_name = rip_hmm[oid]
                cat = cds_category(hmm_name)
                if cat == "other":
                    cat = "rep"
                orf["product"]  = hmm_name
                orf["category"] = cat
                orf["origin_relevant"] = True
                continue

            # Priority 2: Pfam-A hit → real functional annotation
            if oid in pfam_hits:
                pfam_name, pfam_acc, desc = pfam_hits[oid]
                # Use description as product name, pfam_name for category
                product = f"{pfam_name}: {desc}" if desc else pfam_name
                cat = cds_category(product)
                if cat == "other":
                    # Check specific Pfam families
                    if pfam_is_origin_relevant(pfam_name):
                        cat = "structural"
                orf["product"]  = product
                orf["category"] = cat
                orf["pfam"]     = pfam_name
                orf["pfam_acc"] = pfam_acc
                orf["origin_relevant"] = pfam_is_origin_relevant(pfam_name)
                continue

            # Priority 3: No hit → hypothetical
            aa_len = len(orf["aa"])
            orf["product"]  = f"hypothetical protein ({aa_len} aa)"
            orf["category"] = "hypothetical"
            orf["origin_relevant"] = False


def trim_to_origin_core(seqs: dict[str, str],
                         all_orfs: dict[str, list[dict]],
                         buffer: int = BUFFER_BP) -> dict[str, str]:
    """
    Trim origins to only include the core origin region:
    from the outermost origin-relevant ORF + buffer on each side.

    Non-origin ORFs (phage, metabolism, AMR, etc.) at the edges are removed.
    Hypothetical proteins adjacent to origin-relevant ORFs are kept.
    """
    trimmed = {}
    for acc, seq in seqs.items():
        orfs = all_orfs.get(acc, [])
        if not orfs:
            trimmed[acc] = seq
            continue

        # Find origin-relevant ORFs
        relevant = [o for o in orfs if o.get("origin_relevant")]
        if not relevant:
            # Keep Rep at minimum (should always exist)
            relevant = [o for o in orfs if o.get("category") == "rep"]
        if not relevant:
            trimmed[acc] = seq
            continue

        # Core region: outermost relevant ORFs + buffer
        core_start = min(o["start"] for o in relevant)
        core_end   = max(o["end"]   for o in relevant)

        # Also include hypothetical ORFs that are BETWEEN relevant ORFs
        # (they may be unannotated origin components)
        for o in orfs:
            if o["start"] >= core_start and o["end"] <= core_end:
                o["origin_relevant"] = True  # inside core → keep

        # Add buffer
        trim_start = max(0, core_start - buffer)
        trim_end   = min(len(seq), core_end + buffer)

        # Don't trim if it would remove less than 500bp from each side
        # (not worth it for small origins)
        if (trim_start > 500) or (len(seq) - trim_end > 500):
            new_seq = seq[trim_start:trim_end]
            # Adjust ORF coordinates
            new_orfs = []
            for o in orfs:
                if o["end"] <= trim_start or o["start"] >= trim_end:
                    continue  # ORF entirely outside trimmed region
                new_o = dict(o)
                new_o["start"] = max(0, o["start"] - trim_start)
                new_o["end"]   = min(len(new_seq), o["end"] - trim_start)
                new_orfs.append(new_o)
            all_orfs[acc] = new_orfs
            trimmed[acc] = new_seq
            if len(new_seq) < len(seq):
                print(f"  {acc}: trimmed {len(seq):,} -> {len(new_seq):,} bp "
                      f"(removed {len(seq)-len(new_seq):,} bp of non-origin flanking)")
        else:
            trimmed[acc] = seq

    return trimmed


def check_edge_buffer(orfs: list[dict], seq_len: int,
                      buffer: int = BUFFER_BP) -> tuple[int, int]:
    """
    Check if any functionally important ORF is within `buffer` bp of the
    sequence edges.  Returns (need_left, need_right) — how many extra bp
    are needed on each side (0 if no extension required).
    """
    need_left = 0
    need_right = 0
    for orf in orfs:
        if orf.get("category") not in IMPORTANT_CATS:
            continue
        # Left edge: ORF starts within buffer of position 0
        if orf["start"] < buffer:
            shortfall = buffer - orf["start"]
            need_left = max(need_left, shortfall)
        # Right edge: ORF ends within buffer of the sequence end
        if orf["end"] > seq_len - buffer:
            shortfall = buffer - (seq_len - orf["end"])
            need_right = max(need_right, shortfall)
    return need_left, need_right


def extend_origin(acc: str, origin_seq: str, meta: dict,
                  plasmid_seqs: dict[str, str],
                  need_left: int, need_right: int) -> str:
    """
    Extend the origin window using the full plasmid sequence.
    Handles circular wrapping: if the origin is at position 0 and needs
    left extension, wraps to the end of the plasmid.
    Returns the extended sequence.
    """
    full_seq = plasmid_seqs.get(acc, "")
    if not full_seq:
        print(f"    WARNING: no full plasmid for {acc}, cannot extend")
        return origin_seq

    ws = int(meta.get("Origin_start_bp", 0))
    we = int(meta.get("Origin_end_bp", len(origin_seq)))
    plen = len(full_seq)

    # Circular extension: wrap around if at boundaries
    left_ext = ""
    right_ext = ""

    if need_left > 0:
        if ws >= need_left:
            left_ext = full_seq[ws - need_left: ws]
        else:
            # Wrap: take from end of plasmid + whatever is before ws
            shortfall = need_left - ws
            left_ext = full_seq[plen - shortfall:] + full_seq[:ws]

    if need_right > 0:
        if we + need_right <= plen:
            right_ext = full_seq[we: we + need_right]
        else:
            # Wrap: take remainder + from start of plasmid
            shortfall = need_right - (plen - we)
            right_ext = full_seq[we:] + full_seq[:shortfall]

    return left_ext + origin_seq + right_ext


def build_genbank_record(acc: str, seq: str, orfs: list[dict],
                          meta: dict) -> SeqRecord:
    """Build an annotated SeqRecord for one origin."""
    record = SeqRecord(
        Seq(seq),
        id=acc,
        name=acc.replace(".", "_")[:16],
        description=(f"Plasmid origin region | {meta.get('Species', 'unknown')} | "
                     f"RIP: {meta.get('Origin_RIP_type', '?')} | "
                     f"{len(seq):,} bp"),
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["data_file_division"] = "SYN"
    record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()
    record.annotations["source"] = str(meta.get("Species", "unknown organism"))
    record.annotations["organism"] = str(meta.get("Species", "unknown organism"))
    record.annotations["taxonomy"] = [
        str(v) if pd.notna(v) else ""
        for v in [meta.get("Phylum", ""), meta.get("Class", ""),
                  meta.get("Order", ""), meta.get("Family", ""),
                  meta.get("Genus", "")]
    ]

    # Source feature
    source = SeqFeature(
        location=FeatureLocation(0, len(seq)),
        type="source",
        qualifiers={
            "organism":    [str(meta.get("Species", "unknown organism"))],
            "mol_type":    ["genomic DNA"],
            "plasmid":     [acc],
            "note":        [f"Origin region from PLSDB 2025 | "
                            f"RIP type: {meta.get('Origin_RIP_type', '?')} | "
                            f"Tier: {meta.get('Selection_tier', '?')}"],
        },
    )
    record.features.append(source)

    # CDS features for every ORF
    for orf in orfs:
        strand = orf["strand"]
        loc = FeatureLocation(orf["start"], orf["end"], strand=strand)
        cat = orf.get("category", "other")
        color = CATEGORY_COLORS.get(cat, "#78909C")
        notes = [f"category={cat}; color={color}"]
        if orf.get("pfam"):
            notes.append(f"Pfam={orf['pfam']} ({orf.get('pfam_acc', '')})")
        if orf.get("origin_relevant"):
            notes.append("origin_relevant=yes")
        qualifiers = {
            "product":      [orf["product"]],
            "locus_tag":    [orf["id"]],
            "translation":  [orf["aa"]],
            "codon_start":  [1],
            "function":     [cat],
            "note":         notes,
        }
        if orf.get("pfam_acc"):
            qualifiers["db_xref"] = [f"Pfam:{orf['pfam_acc']}"]
        # Add ApE/SnapGene-compatible color
        qualifiers["ApEinfo_fwdcolor"] = [color]
        qualifiers["ApEinfo_revcolor"] = [color]
        feature = SeqFeature(location=loc, type="CDS", qualifiers=qualifiers)
        record.features.append(feature)

    return record


def load_plsdb_plasmids(fasta_path: str, needed_accs: set[str]) -> dict[str, str]:
    """Load only the needed plasmid sequences from the PLSDB FASTA."""
    print(f"  Loading full plasmid sequences from {fasta_path} ...")
    plasmids = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        acc = rec.id.split()[0]
        if acc in needed_accs:
            plasmids[acc] = str(rec.seq).upper()
        else:
            acc_nv = acc.rsplit(".", 1)[0]
            if acc_nv in needed_accs:
                plasmids[acc_nv] = str(rec.seq).upper()
        if len(plasmids) >= len(needed_accs):
            break
    print(f"  Loaded {len(plasmids)} plasmid sequences")
    return plasmids


def main():
    p = argparse.ArgumentParser(
        description="Convert origin FASTA to annotated GenBank files with Pfam annotation and trimming")
    p.add_argument("--fasta",       required=True, help="Path to origin_sequences.fasta")
    p.add_argument("--csv",         required=True, help="Path to diversity_library_metadata.csv")
    p.add_argument("--rip-hmm",     default="RIPs/RIP.hmm", help="Path to RIP HMM profile")
    p.add_argument("--pfam-hmm",    default="Pfam-A.hmm", help="Path to Pfam-A.hmm (for ORF annotation)")
    p.add_argument("--plsdb-fasta", default=None,
                   help="Path to PLSDB sequences.fasta (needed for window extension)")
    p.add_argument("--outdir",      default="results_plsdb/genbanks", help="Output directory")
    p.add_argument("--buffer",      type=int, default=200,
                   help="Minimum bp between window edge and important ORF (default: 200)")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    buf = args.buffer

    # Load origin sequences
    print(f"Loading sequences from {args.fasta} ...")
    seqs = {}
    for rec in SeqIO.parse(args.fasta, "fasta"):
        seqs[rec.id] = str(rec.seq).upper()
    print(f"  {len(seqs)} origin sequences")

    # Load metadata
    print(f"Loading metadata from {args.csv} ...")
    df = pd.read_csv(args.csv)
    meta_lookup = df.set_index("Accession").to_dict("index")
    print(f"  {len(df)} rows")

    # Pass 1: predict ORFs
    print("Pass 1: Predicting ORFs ...")
    all_orfs = {acc: predict_orfs(seq, acc) for acc, seq in seqs.items()}

    # Classify with RIP HMM (Rep proteins)
    print("  Running RIP HMM classification ...")
    rip_hits = classify_orfs_rip(all_orfs, args.rip_hmm)

    # Classify with Pfam-A (all other functions)
    print("  Running Pfam-A classification ...")
    pfam_hits = classify_orfs_pfam(all_orfs, args.pfam_hmm)

    # Label all ORFs
    label_orfs(all_orfs, rip_hits, pfam_hits)

    total_orfs = sum(len(v) for v in all_orfs.values())
    n_rep      = sum(1 for orfs in all_orfs.values() for o in orfs if o.get("category") == "rep")
    n_pfam     = sum(1 for orfs in all_orfs.values() for o in orfs if o.get("pfam"))
    n_hyp      = sum(1 for orfs in all_orfs.values() for o in orfs if o.get("category") == "hypothetical")
    n_relevant = sum(1 for orfs in all_orfs.values() for o in orfs if o.get("origin_relevant"))
    print(f"  {total_orfs} ORFs total: {n_rep} Rep, {n_pfam} Pfam-annotated, "
          f"{n_hyp} hypothetical, {n_relevant} origin-relevant")

    # Trim origins: remove non-origin flanking ORFs from 6kb windows
    print("\nTrimming origins to origin-relevant core ...")
    seqs = trim_to_origin_core(seqs, all_orfs, buf)

    important  = sum(1 for orfs in all_orfs.values()
                     for o in orfs if o.get("category") in IMPORTANT_CATS)

    # Check which origins need extension
    need_extension: dict[str, tuple[int, int]] = {}
    for acc, orfs in all_orfs.items():
        nl, nr = check_edge_buffer(orfs, len(seqs[acc]), buf)
        if nl > 0 or nr > 0:
            need_extension[acc] = (nl, nr)

    if need_extension:
        print(f"\n{len(need_extension)} origins need window extension for {buf}bp ORF buffer:")
        for acc, (nl, nr) in sorted(need_extension.items()):
            print(f"  {acc}: +{nl}bp left, +{nr}bp right")

        if args.plsdb_fasta and os.path.isfile(args.plsdb_fasta):
            # Load all needed plasmid sequences (may need more in later rounds)
            all_needed = set(seqs.keys())
            plasmid_seqs = load_plsdb_plasmids(args.plsdb_fasta, all_needed)

            # Iteratively extend until no important ORF is within buffer of edge
            for round_num in range(1, 6):  # max 5 rounds
                extended_any = False
                for acc, (nl, nr) in need_extension.items():
                    meta = meta_lookup.get(acc, {})
                    new_seq = extend_origin(acc, seqs[acc], meta, plasmid_seqs, nl, nr)
                    if len(new_seq) != len(seqs[acc]):
                        if round_num == 1:
                            print(f"  {acc}: {len(seqs[acc]):,} -> {len(new_seq):,} bp")
                        seqs[acc] = new_seq
                        # Update meta for next round's extend_origin
                        # After circular extension, track new absolute positions
                        ws = int(meta.get("Origin_start_bp", 0))
                        we = int(meta.get("Origin_end_bp", 0))
                        plen = len(plasmid_seqs.get(acc, ""))
                        meta["Origin_start_bp"] = (ws - nl) % plen if plen else max(0, ws - nl)
                        meta["Origin_end_bp"]   = (we + nr) % plen if plen else we + nr
                        meta_lookup[acc] = meta
                        extended_any = True

                # Re-predict and re-classify on all origins (extended ones changed)
                print(f"\nPass {round_num + 1}: Re-predicting ORFs on extended origins ...")
                for acc in need_extension:
                    all_orfs[acc] = predict_orfs(seqs[acc], acc)
                rip_hits = classify_orfs_rip(all_orfs, args.rip_hmm)
                pfam_hits = classify_orfs_pfam(all_orfs, args.pfam_hmm)
                label_orfs(all_orfs, rip_hits, pfam_hits)

                # Check again
                need_extension = {}
                for acc, orfs in all_orfs.items():
                    nl, nr = check_edge_buffer(orfs, len(seqs[acc]), buf)
                    if nl > 0 or nr > 0:
                        need_extension[acc] = (nl, nr)

                if not need_extension:
                    print(f"  All origins now have >= {buf}bp buffer. Done after {round_num} round(s).")
                    break
                print(f"  {len(need_extension)} origins still need extension ...")
                if not extended_any:
                    print(f"  Cannot extend further (at plasmid boundaries). Accepting as-is.")
                    break
            else:
                if need_extension:
                    print(f"  WARNING: {len(need_extension)} origins still within {buf}bp "
                          f"after 5 rounds (at plasmid boundaries)")
        else:
            print(f"\n  WARNING: --plsdb-fasta not provided or not found.")
            print(f"  Cannot extend windows. Provide --plsdb-fasta to enable extension.")
    else:
        print(f"\nAll origins have >= {buf}bp buffer around important ORFs. No extension needed.")

    # Write GenBank files
    print(f"\nWriting GenBank files to {outdir}/ ...")
    written = 0
    for acc, seq in seqs.items():
        orfs = all_orfs.get(acc, [])
        meta = meta_lookup.get(acc, {})
        record = build_genbank_record(acc, seq, orfs, meta)
        gb_path = outdir / f"{acc}.gb"
        with open(gb_path, "w") as fh:
            SeqIO.write(record, fh, "genbank")
        written += 1

    print(f"\nDone: {written} GenBank files written to {outdir}/")
    n_cats = {}
    for orfs in all_orfs.values():
        for o in orfs:
            cat = o.get("category", "other")
            n_cats[cat] = n_cats.get(cat, 0) + 1
    print("ORF categories across all origins:")
    for cat in sorted(n_cats, key=lambda c: -n_cats[c]):
        print(f"  {cat:15s}: {n_cats[cat]:4d}  {CATEGORY_COLORS.get(cat, '')}")

    # ── Regenerate FASTA and PDF with extended sequences + classified ORFs ──
    results_dir = Path(args.fasta).parent

    # Update FASTA
    fasta_out = results_dir / "origin_sequences.fasta"
    print(f"\nUpdating {fasta_out} with extended sequences ...")
    with open(fasta_out, "w") as fh:
        for acc, seq in seqs.items():
            fh.write(f">{acc}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    print(f"  {len(seqs)} sequences written")

    # Regenerate PDF with classified ORFs
    try:
        import sys, importlib
        # Import pipeline functions
        sys.path.insert(0, str(Path(__file__).parent))
        import plasmid_origin_pipeline as pp
        importlib.reload(pp)

        # Build gb_map from our classified ORFs (same format pipeline expects)
        gb_map = {}
        for acc, orfs in all_orfs.items():
            gb_map[acc] = [{
                "product":  o["product"],
                "start":    o["start"],
                "end":      o["end"],
                "strand":   o["strand"],
                "category": o.get("category", "other"),
            } for o in orfs]

        # Load validation data
        val_csv = results_dir / "diversity_library_metadata.csv"
        lib_df  = pd.read_csv(val_csv)

        # Run validation to get val_df
        lib_seqs_val = {acc: seqs[acc] for acc in lib_df["Accession"] if acc in seqs}
        all_orfs_pipe = pp.predict_orfs_for_library(lib_seqs_val)
        val_df = pp.validate_library(lib_df, lib_seqs_val, {}, args.rip_hmm)

        # Update Origin_span_bp in lib_df to match extended sequences
        lib_df["Origin_span_bp"] = lib_df["Accession"].map(
            lambda a: len(seqs.get(a, "")))

        pdf_path = results_dir / "origin_plasmid_maps.pdf"
        print(f"\nRegenerating {pdf_path} with classified ORFs ...")
        pp.generate_pdf(lib_df, seqs, {}, all_orfs_pipe, gb_map, val_df, str(pdf_path))

        # Regenerate overview
        overview_path = results_dir / "diversity_library_overview.png"
        pp.generate_overview(lib_df, val_df, str(overview_path))

    except Exception as e:
        import traceback
        print(f"\nWARNING: Could not regenerate PDF: {e}")
        traceback.print_exc()
        print("GenBank files are still valid. PDF can be regenerated by re-running the pipeline.")


if __name__ == "__main__":
    main()
