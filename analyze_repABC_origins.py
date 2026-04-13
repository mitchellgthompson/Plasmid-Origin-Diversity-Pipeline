#!/usr/bin/env python3
"""
Comprehensive repABC origin analysis pipeline.

Parses all .fna/.gbk files from ARC_repABC_loci_fna_gbk/, integrates
parS data, computes functional/taxonomic metadata, runs Twist synthesis
assessment, generates diversity statistics, and prioritises origins into
primary (~125 kb) and secondary synthesis tiers.

Usage:
    python3 analyze_repABC_origins.py
"""

import os, re, glob, math, warnings
from pathlib import Path
from collections import Counter, defaultdict
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

warnings.filterwarnings("ignore")

# ── Feature colour palette ───────────────────────────────────────────────────
FEATURE_COLORS = {
    "repA":         "#E53935",   # Red
    "repB":         "#1565C0",   # Blue
    "repC":         "#2E7D32",   # Green
    "rep_origin":   "#FF8F00",   # Amber
    "ctRNA":        "#7B1FA2",   # Purple
    "S-element":    "#00897B",   # Teal
    "parS":         "#F06292",   # Pink
    "other_CDS":    "#78909C",   # Steel gray
}

# ── paths ────────────────────────────────────────────────────────────────────
BASE     = Path(__file__).resolve().parent
ARC_DIR  = BASE / "ARC_repABC_loci_fna_gbk"
PARS_XLS = ARC_DIR / "parS.table.xlsx"
OUT_DIR  = BASE / "results_repABC"
OUT_DIR.mkdir(exist_ok=True)


# ═════════════════════════════════════════════════════════════════════════════
# 1.  PARSE GBK + FNA FILES
# ═════════════════════════════════════════════════════════════════════════════
print("─── Step 1: Parsing GBK and FNA files ───")

fna_files = sorted(glob.glob(str(ARC_DIR / "*.fna")))
gbk_files = sorted(glob.glob(str(ARC_DIR / "*.gbk")))
print(f"  Found {len(fna_files)} .fna and {len(gbk_files)} .gbk files")

# Build lookup: stem -> paths
stems = set()
for f in fna_files:
    stems.add(Path(f).stem)

rows = []
sequences = {}    # locus_id -> DNA string
gbk_features = {} # locus_id -> list of feature dicts for drawing

for fna_path in fna_files:
    stem = Path(fna_path).stem  # e.g. C58_AE007871.2_pTi_Type_I.a_repABC
    gbk_path = str(Path(fna_path).with_suffix(".gbk"))

    # ── read sequence ────────────────────────────────────────────────────
    rec_fna = next(SeqIO.parse(fna_path, "fasta"))
    seq = str(rec_fna.seq).upper()
    seq_len = len(seq)
    gc = gc_fraction(rec_fna.seq)

    # ── parse filename for strain / replicon / pTi type ──────────────────
    # Strip trailing _repABC
    base = stem.replace("_repABC", "")

    # Detect replicon type from filename
    replicon_type = "unknown"
    pti_type = None
    pri_type = None

    if "_chromid" in base:
        replicon_type = "chromid"
    elif "_pAt" in base:
        replicon_type = "pAt"
    elif "_pAg" in base:
        replicon_type = "pAg"
    elif "_pRi_" in base:
        replicon_type = "pRi"
        m = re.search(r"_pRi_(Type_\w+)", base)
        if m:
            pri_type = m.group(1).replace("_", " ")
    elif "_pTi_" in base:
        replicon_type = "pTi"
        m = re.search(r"_pTi_(Type_[\w.]+)", base)
        if m:
            pti_type = m.group(1).replace("_", " ")

    # Detect rep copy number (rep_1, rep_2)
    rep_copy = None
    m_rep = re.search(r"_rep_(\d+)_repABC", stem)
    if m_rep:
        rep_copy = int(m_rep.group(1))

    # Extract strain name (everything before the first accession or replicon tag)
    # Heuristic: take prefix before first _CP, _NZ_, _AE, _CA, _contig_, _chromid, _pAt, _pAg, _plasmid
    strain = re.split(
        r"_(?:CP|NZ_|AE|CA|CAICSX|JAAM|U\d{5}|contig_|chromid|pAt|pAg|plasmid|CG\d{3}_\d|CG\d{3}_CG)",
        base
    )[0]

    # ── parse GBK for functional annotations ─────────────────────────────
    has_repA = has_repB = has_repC = False
    has_rep_origin = False
    has_ctRNA = has_S_element = False
    repA_pseudo = False
    repA_len = repB_len = repC_len = 0
    repA_product = repB_product = repC_product = ""
    repA_uniref = repB_uniref = repC_uniref = ""
    n_cds = 0
    rep_origin_len = 0
    organism = ""
    definition = ""

    if os.path.exists(gbk_path):
        rec_gbk = next(SeqIO.parse(gbk_path, "genbank"))
        organism = rec_gbk.annotations.get("organism", "")
        definition = rec_gbk.description or ""

        for feat in rec_gbk.features:
            if feat.type == "CDS":
                n_cds += 1
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                pseudo = "pseudogene" in feat.qualifiers
                cds_len = len(feat.location)

                # Extract UniRef50 from db_xref
                uniref = ""
                for xr in feat.qualifiers.get("db_xref", []):
                    if "UniRef50" in xr:
                        uniref = xr.split(":")[-1] if ":" in xr else xr

                if gene == "repA":
                    has_repA = True
                    repA_len = cds_len
                    repA_product = product
                    repA_pseudo = pseudo
                    repA_uniref = uniref
                elif gene == "repB":
                    has_repB = True
                    repB_len = cds_len
                    repB_product = product
                    repB_uniref = uniref
                elif gene == "repC":
                    has_repC = True
                    repC_len = cds_len
                    repC_product = product
                    repC_uniref = uniref

            elif feat.type == "rep_origin":
                has_rep_origin = True
                rep_origin_len = len(feat.location)

            elif feat.type == "ncRNA":
                ncclass = feat.qualifiers.get("ncRNA_class", [""])[0]
                note = " ".join(feat.qualifiers.get("note", []))
                if "ctRNA" in note or "ctRNA" in str(feat.qualifiers):
                    has_ctRNA = True

            elif feat.type == "regulatory":
                note = " ".join(feat.qualifiers.get("note", []))
                reg_class = feat.qualifiers.get("regulatory_class", [""])[0]
                if "S-element" in note or "S-element" in reg_class:
                    has_S_element = True

    # ── collect drawable features ────────────────────────────────────────
    feats_draw = []
    if os.path.exists(gbk_path):
        rec_gbk2 = next(SeqIO.parse(gbk_path, "genbank"))
        for feat in rec_gbk2.features:
            start = int(feat.location.start)
            end   = int(feat.location.end)
            strand = feat.location.strand or 1

            if feat.type == "CDS":
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                pseudo = "pseudogene" in feat.qualifiers
                if gene == "repA":
                    cat = "repA"
                    label = "repA (pseudo)" if pseudo else "repA"
                elif gene == "repB":
                    cat = "repB"
                    label = "repB"
                elif gene == "repC":
                    cat = "repC"
                    label = "repC"
                else:
                    cat = "other_CDS"
                    label = product[:30] if product else gene[:30] if gene else "CDS"
                feats_draw.append({"start": start, "end": end, "strand": strand,
                                   "category": cat, "label": label})

            elif feat.type == "rep_origin":
                feats_draw.append({"start": start, "end": end, "strand": 0,
                                   "category": "rep_origin", "label": "oriV"})

            elif feat.type == "ncRNA":
                note = " ".join(feat.qualifiers.get("note", []))
                if "ctRNA" in note or "ctRNA" in str(feat.qualifiers):
                    feats_draw.append({"start": start, "end": end, "strand": strand,
                                       "category": "ctRNA", "label": "ctRNA"})

            elif feat.type == "regulatory":
                note = " ".join(feat.qualifiers.get("note", []))
                if "S-element" in note:
                    feats_draw.append({"start": start, "end": end, "strand": strand,
                                       "category": "S-element", "label": "S-element"})

    locus_id = stem
    gbk_features[locus_id] = feats_draw

    # ── completeness flag ────────────────────────────────────────────────
    repABC_complete = has_repA and has_repB and has_repC
    functional_elements = sum([has_rep_origin, has_ctRNA, has_S_element])

    sequences[locus_id] = seq

    rows.append({
        "locus_id": locus_id,
        "strain": strain,
        "replicon_type": replicon_type,
        "pTi_type": pti_type,
        "pRi_type": pri_type,
        "rep_copy": rep_copy,
        "definition": definition,
        "length_bp": seq_len,
        "GC_content": round(gc, 4),
        "n_CDS": n_cds,
        "has_repA": has_repA,
        "has_repB": has_repB,
        "has_repC": has_repC,
        "repABC_complete": repABC_complete,
        "repA_pseudogene": repA_pseudo,
        "repA_len_bp": repA_len,
        "repB_len_bp": repB_len,
        "repC_len_bp": repC_len,
        "repA_product": repA_product,
        "repB_product": repB_product,
        "repC_product": repC_product,
        "repA_UniRef50": repA_uniref,
        "repB_UniRef50": repB_uniref,
        "repC_UniRef50": repC_uniref,
        "has_rep_origin": has_rep_origin,
        "rep_origin_len_bp": rep_origin_len,
        "has_ctRNA": has_ctRNA,
        "has_S_element": has_S_element,
        "n_regulatory_elements": functional_elements,
    })

df = pd.DataFrame(rows)
print(f"  Parsed {len(df)} repABC loci")
print(f"  Complete repABC: {df['repABC_complete'].sum()}/{len(df)}")
print(f"  Replicon types: {dict(df['replicon_type'].value_counts())}")


# ═════════════════════════════════════════════════════════════════════════════
# 2.  INTEGRATE parS TABLE
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 2: Integrating parS data ───")

if PARS_XLS.exists():
    pars_df = pd.read_excel(PARS_XLS)
    print(f"  parS table: {len(pars_df)} rows, columns: {list(pars_df.columns)}")

    # The parS table likely has a column matching locus names
    # Try to match on the first column or a column containing locus identifiers
    pars_cols = list(pars_df.columns)
    print(f"  First 3 rows:\n{pars_df.head(3).to_string()}")

    # Determine matching column - try first column
    pars_id_col = pars_cols[0]

    # Filter to actual parS hits (rows where parS filter is not NaN)
    pars_hits = pars_df.dropna(subset=["parS filter"]).copy()
    print(f"  parS hits (non-NaN filter): {len(pars_hits)}")

    # Count parS sites per locus
    pars_per_locus = pars_hits.groupby(pars_id_col).agg(
        n_parS=pd.NamedAgg(column=pars_id_col, aggfunc="count"),
    ).reset_index()

    # Also collect parS motifs if available
    motif_col = None
    for c in pars_cols:
        if any(k in c.lower() for k in ["motif", "sequence", "pattern", "pars"]):
            if c != pars_id_col:
                motif_col = c
                break

    if motif_col is None and len(pars_cols) > 1:
        # Try the last column as motif/sequence
        motif_col = pars_cols[-1]

    if motif_col:
        motif_agg = pars_hits.groupby(pars_id_col)[motif_col].apply(
            lambda x: "; ".join(str(v) for v in x.unique())
        ).reset_index()
        motif_agg.columns = [pars_id_col, "parS_motifs"]
        pars_per_locus = pars_per_locus.merge(motif_agg, on=pars_id_col, how="left")

    # Also track parS filter types per locus
    filter_agg = pars_hits.groupby(pars_id_col)["parS filter"].apply(
        lambda x: "; ".join(str(v) for v in x.unique())
    ).reset_index()
    filter_agg.columns = [pars_id_col, "parS_filter_types"]
    pars_per_locus = pars_per_locus.merge(filter_agg, on=pars_id_col, how="left")

    # Collect start/end positions if available
    start_col = end_col = strand_col = None
    for c in pars_cols:
        cl = c.lower()
        if "start" in cl: start_col = c
        elif "end" in cl: end_col = c
        elif "strand" in cl: strand_col = c

    if start_col and end_col:
        pos_agg = pars_hits.dropna(subset=[start_col, end_col]).groupby(pars_id_col).apply(
            lambda g: "; ".join(
                f"{int(s)}-{int(e)}" for s, e in zip(g[start_col], g[end_col])
            )
        ).reset_index()
        pos_agg.columns = [pars_id_col, "parS_positions"]
        pars_per_locus = pars_per_locus.merge(pos_agg, on=pars_id_col, how="left")

    pars_per_locus = pars_per_locus.rename(columns={pars_id_col: "locus_id"})

    # Merge - try exact match first, then fuzzy
    pre_len = len(df)
    df = df.merge(pars_per_locus, on="locus_id", how="left")

    # If no matches, try matching on stem without _repABC
    matched = df["n_parS"].notna().sum()
    if matched == 0:
        print("  Direct match failed, trying fuzzy locus matching...")
        df = df.drop(columns=[c for c in pars_per_locus.columns if c != "locus_id"], errors="ignore")

        # Try matching parS IDs to our locus_ids
        pars_ids = set(pars_per_locus["locus_id"].astype(str))
        locus_ids = set(df["locus_id"])

        # Build mapping: try substring matching
        mapping = {}
        for lid in locus_ids:
            lid_base = lid.replace("_repABC", "")
            for pid in pars_ids:
                pid_clean = str(pid).strip()
                if pid_clean == lid or pid_clean == lid_base or lid_base in pid_clean or pid_clean in lid_base:
                    mapping[lid] = pid_clean
                    break

        if mapping:
            pars_per_locus_mapped = pars_per_locus.copy()
            reverse_map = {v: k for k, v in mapping.items()}
            pars_per_locus_mapped["locus_id"] = pars_per_locus_mapped["locus_id"].astype(str).map(
                lambda x: reverse_map.get(x.strip(), x)
            )
            df = df.merge(pars_per_locus_mapped, on="locus_id", how="left")
            matched = df["n_parS"].notna().sum()

    df["n_parS"] = df["n_parS"].fillna(0).astype(int)
    df["has_parS"] = df["n_parS"] > 0
    print(f"  Matched parS data for {matched}/{len(df)} loci")
    print(f"  Loci with parS sites: {df['has_parS'].sum()}")
else:
    print("  WARNING: parS.table.xlsx not found")
    df["n_parS"] = 0
    df["has_parS"] = False


# ═════════════════════════════════════════════════════════════════════════════
# 3.  TAXONOMY (from DEFINITION lines and strain names)
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 3: Assigning taxonomy ───")

# Manual taxonomy from known strains in the Agrobacterium/Rhizobium complex
TAXONOMY = {
    "C58":       ("Agrobacterium fabrum",       "Rhizobiaceae", "Alphaproteobacteria"),
    "1D1514":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "1D1524":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "15-1187-1-2a": ("Agrobacterium tumefaciens", "Rhizobiaceae", "Alphaproteobacteria"),
    "Atu_K84":   ("Agrobacterium radiobacter",  "Rhizobiaceae", "Alphaproteobacteria"),
    "Atu_B6":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "AMS008":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "K599":      ("Agrobacterium rhizogenes",   "Rhizobiaceae", "Alphaproteobacteria"),
    "NCIB_8196": ("Agrobacterium rhizogenes",   "Rhizobiaceae", "Alphaproteobacteria"),
    "T155_95":   ("Agrobacterium rhizogenes",   "Rhizobiaceae", "Alphaproteobacteria"),
    "NCPPB_1641":("Agrobacterium vitis",        "Rhizobiaceae", "Alphaproteobacteria"),
    "CG474":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "CG507":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "CG958":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "CG1055":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "AB2_73":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "AF3.44":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "B21_90":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "B230_85":   ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "CA75_95":   ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "L51_94":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "Q15_94":    ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "S7_73":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "WS16":      ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "WS163":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "WS192":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "WS198":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "WS287":     ("Agrobacterium tumefaciens",  "Rhizobiaceae", "Alphaproteobacteria"),
    "Rhizobium_etli_CFN_42": ("Rhizobium etli CFN 42", "Rhizobiaceae", "Alphaproteobacteria"),
    "Sinorhizobium_fredii_NGR234": ("Sinorhizobium fredii NGR234", "Rhizobiaceae", "Alphaproteobacteria"),
}

def assign_taxonomy(strain, definition):
    """Look up taxonomy from strain name, falling back to definition."""
    # Try exact match
    if strain in TAXONOMY:
        return TAXONOMY[strain]
    # Try prefix match
    for key in TAXONOMY:
        if strain.startswith(key):
            return TAXONOMY[key]
    # Parse definition for genus
    for genus in ["Agrobacterium", "Rhizobium", "Sinorhizobium", "Mesorhizobium"]:
        if genus.lower() in definition.lower():
            return (definition.split(",")[0].strip(), "Rhizobiaceae", "Alphaproteobacteria")
    return ("Agrobacterium sp.", "Rhizobiaceae", "Alphaproteobacteria")

tax_data = df.apply(lambda r: assign_taxonomy(r["strain"], r["definition"]), axis=1)
df["species"] = [t[0] for t in tax_data]
df["family"]  = [t[1] for t in tax_data]
df["class"]   = [t[2] for t in tax_data]

print(f"  Species: {dict(df['species'].value_counts())}")


# ═════════════════════════════════════════════════════════════════════════════
# 4.  SEQUENCE COMPLEXITY METRICS
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 4: Computing sequence complexity metrics ───")

def compute_complexity(seq):
    """Compute repeat/complexity metrics relevant to synthesis."""
    seq = seq.upper()
    n = len(seq)

    # Homopolymer runs
    max_homopoly = 0
    n_homopoly_8plus = 0
    for m in re.finditer(r"([ACGT])\1{7,}", seq):
        run_len = len(m.group())
        max_homopoly = max(max_homopoly, run_len)
        if run_len >= 8:
            n_homopoly_8plus += 1

    # Direct repeats (≥20 bp)
    n_direct_repeats = 0
    for size in [20, 30, 50]:
        kmers = defaultdict(int)
        for i in range(n - size + 1):
            kmers[seq[i:i+size]] += 1
        n_direct_repeats += sum(1 for v in kmers.values() if v > 1)

    # GC content in sliding windows (detect extreme local GC)
    win = 100
    gc_vals = []
    for i in range(0, n - win + 1, 20):
        w = seq[i:i+win]
        gc_vals.append((w.count("G") + w.count("C")) / len(w))
    gc_min_local = min(gc_vals) if gc_vals else 0
    gc_max_local = max(gc_vals) if gc_vals else 0
    gc_range_local = gc_max_local - gc_min_local

    # Dinucleotide entropy (measure of sequence randomness)
    dinucs = [seq[i:i+2] for i in range(n-1)]
    dc = Counter(dinucs)
    total = sum(dc.values())
    entropy = -sum((c/total) * math.log2(c/total) for c in dc.values() if c > 0)
    max_entropy = math.log2(16)  # 4 bp alphabet, 16 dinucleotides
    complexity_score = entropy / max_entropy  # 0-1, higher = more complex = easier

    return {
        "max_homopolymer_bp": max_homopoly,
        "n_homopolymer_8plus": n_homopoly_8plus,
        "n_direct_repeats_20bp": n_direct_repeats,
        "local_GC_min": round(gc_min_local, 4),
        "local_GC_max": round(gc_max_local, 4),
        "local_GC_range": round(gc_range_local, 4),
        "dinuc_complexity": round(complexity_score, 4),
    }

complexity_rows = []
for _, row in df.iterrows():
    seq = sequences[row["locus_id"]]
    metrics = compute_complexity(seq)
    metrics["locus_id"] = row["locus_id"]
    complexity_rows.append(metrics)

complexity_df = pd.DataFrame(complexity_rows)
df = df.merge(complexity_df, on="locus_id")

print(f"  Max homopolymer range: {df['max_homopolymer_bp'].min()}-{df['max_homopolymer_bp'].max()} bp")
print(f"  Dinucleotide complexity range: {df['dinuc_complexity'].min():.3f}-{df['dinuc_complexity'].max():.3f}")


# ═════════════════════════════════════════════════════════════════════════════
# 5.  TWIST SYNTHESIS ASSESSMENT
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 5: Running Twist synthesis assessment ───")

def _gc_cat(gc):
    if 0.25 <= gc <= 0.65: return "Optimal"
    if 0.20 <= gc <= 0.70: return "Challenging"
    if 0.15 <= gc <= 0.75: return "Difficult"
    return "Very Difficult"

def _strat(bp):
    if bp <= 1800: return "Single fragment"
    if bp <= 5000: return "Clonal gene"
    return "Multi-fragment assembly"

def _frags(bp, f=1800):
    if bp <= 5000: return 1
    return int(np.ceil(bp / (f - 40)))

def _cost(bp, gc):
    c = max(bp*0.07, 99) if bp<=1800 else bp*0.10 if bp<=5000 else _frags(bp)*max(1800*0.07,99)*1.1
    if gc < 0.25 or gc > 0.65: c *= 1.3
    if gc < 0.20 or gc > 0.70: c *= 1.5
    return round(c, 2)

def _ease(gc, bp):
    gp = gc * 100
    gs = 50 if 40<=gp<=60 else 40 if 30<=gp<=65 else 25 if 25<=gp<=70 else 10 if 20<=gp<=75 else 0
    ls = 50 if bp<=1800 else 40 if bp<=3000 else 30 if bp<=5000 else 20 if bp<=6000 else 10
    return gs + ls

def _ease_cat(s):
    return "Easy" if s>=80 else "Moderate" if s>=60 else "Challenging" if s>=40 else "Difficult"

df["Twist_GC_cat"]      = df["GC_content"].apply(_gc_cat)
df["Twist_strategy"]    = df["length_bp"].apply(_strat)
df["Twist_n_fragments"] = df["length_bp"].apply(_frags)
df["Twist_est_cost"]    = df.apply(lambda r: _cost(r["length_bp"], r["GC_content"]), axis=1)
df["Twist_ease_score"]  = df.apply(lambda r: _ease(r["GC_content"], r["length_bp"]), axis=1)
df["Twist_ease_cat"]    = df["Twist_ease_score"].apply(_ease_cat)
df["Twist_has_homopolymer"] = df["locus_id"].map(
    lambda a: bool(re.search(r"([ACGT])\1{10,}", sequences.get(a, "").upper())))

# Enhanced synthesis flag: synthesizable = easy/moderate AND no long homopolymers
df["synthesizable"] = (
    df["Twist_ease_cat"].isin(["Easy", "Moderate"]) &
    ~df["Twist_has_homopolymer"] &
    (df["max_homopolymer_bp"] < 11) &
    (df["n_direct_repeats_20bp"] < 10)
)

print(f"  Synthesis ease: {dict(df['Twist_ease_cat'].value_counts())}")
print(f"  Synthesizable: {df['synthesizable'].sum()}/{len(df)}")
print(f"  Total est. cost: ${df['Twist_est_cost'].sum():,.2f}")
print(f"  Strategy: {dict(df['Twist_strategy'].value_counts())}")


# ═════════════════════════════════════════════════════════════════════════════
# 6.  DIVERSITY METRICS & GROUPING
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 6: Computing diversity metrics ───")

# Assign a "repABC class" for diversity grouping
# Priority: pTi_type > pRi_type > replicon_type + UniRef50 cluster
def assign_rep_class(row):
    """Assign a diversity class to each origin."""
    if row["pTi_type"]:
        return f"pTi {row['pTi_type']}"
    if row["pRi_type"]:
        return f"pRi {row['pRi_type']}"
    if row["replicon_type"] == "chromid":
        return f"chromid_{row['repC_UniRef50'][:20]}" if row["repC_UniRef50"] else "chromid"
    if row["replicon_type"] == "pAt":
        return f"pAt_{row['repC_UniRef50'][:20]}" if row["repC_UniRef50"] else "pAt"
    if row["replicon_type"] == "pAg":
        return "pAg"
    # Unknown: group by RepC UniRef50 cluster
    uniref = row["repC_UniRef50"]
    if uniref:
        return f"other_{uniref[:25]}"
    return f"other_{row['strain']}"

df["rep_class"] = df.apply(assign_rep_class, axis=1)

# Count representatives per class
class_counts = df["rep_class"].value_counts()
print(f"  Unique rep classes: {len(class_counts)}")
print(f"  Top 10 classes:")
for cls, cnt in class_counts.head(10).items():
    print(f"    {cls}: {cnt}")

# Functional completeness score (0-100)
def functional_score(row):
    """Score origin by functional completeness."""
    score = 0
    if row["has_repA"]: score += 20
    if row["has_repB"]: score += 20
    if row["has_repC"]: score += 25
    if row["has_rep_origin"]: score += 15
    if row["has_parS"]: score += 10
    if row["has_ctRNA"]: score += 5
    if row["has_S_element"]: score += 5
    if row["repA_pseudogene"]: score -= 15
    return score

df["functional_score"] = df.apply(functional_score, axis=1)

print(f"  Functional score range: {df['functional_score'].min()}-{df['functional_score'].max()}")
print(f"  Mean functional score: {df['functional_score'].mean():.1f}")


# ═════════════════════════════════════════════════════════════════════════════
# 7.  COMPOSITE RANKING & TIER ASSIGNMENT
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 7: Ranking and tier assignment ───")

# Composite priority score (higher = better candidate)
# Weights: functional completeness (40%), synthesis ease (30%), diversity bonus (20%), compactness (10%)
def composite_score(row):
    func = row["functional_score"] / 100.0  # 0-1
    synth = row["Twist_ease_score"] / 100.0  # 0-1
    compact = max(0, 1.0 - (row["length_bp"] - 3000) / 3000)  # bonus for smaller origins
    compact = min(1.0, max(0, compact))
    complexity = row["dinuc_complexity"]  # 0-1, higher = easier to synthesize
    return round(0.35 * func + 0.30 * synth + 0.15 * compact + 0.20 * complexity, 4)

df["priority_score"] = df.apply(composite_score, axis=1)

# Sort by priority score descending
df = df.sort_values("priority_score", ascending=False).reset_index(drop=True)

# ── Assign tiers ─────────────────────────────────────────────────────────
# PRIMARY: synthesizable, diverse, ~125 kb total
# SECONDARY: everything else (backup)

TARGET_BP = 125_000

# Strategy: prioritise rep classes that are most biologically important,
# then fill to budget. Classes are ranked:
#   Tier A: pTi types (each unique Ti plasmid type)
#   Tier B: pRi types (Ri plasmid types)
#   Tier C: named replicon types (chromid, pAt, pAg)
#   Tier D: other rep classes (by priority_score of best member)

synth_df = df[df["synthesizable"]].copy()

# Build best candidate per class
best_per_class = []
for cls in synth_df["rep_class"].unique():
    candidates = synth_df[synth_df["rep_class"] == cls].sort_values("priority_score", ascending=False)
    best = candidates.iloc[0]
    # Assign class priority tier
    if cls.startswith("pTi"):
        class_tier = 0  # highest
    elif cls.startswith("pRi"):
        class_tier = 1
    elif cls.startswith(("chromid", "pAt", "pAg")):
        class_tier = 2
    else:
        class_tier = 3
    best_per_class.append({
        "locus_id": best["locus_id"],
        "rep_class": cls,
        "class_tier": class_tier,
        "priority_score": best["priority_score"],
        "length_bp": best["length_bp"],
    })

best_per_class = pd.DataFrame(best_per_class).sort_values(
    ["class_tier", "priority_score"], ascending=[True, False]
)

# Greedy fill to ~125 kb
primary_ids = []
primary_bp = 0
classes_in_primary = set()

for _, row in best_per_class.iterrows():
    if primary_bp + row["length_bp"] > TARGET_BP * 1.05:  # allow 5% overshoot
        continue
    primary_ids.append(row["locus_id"])
    primary_bp += row["length_bp"]
    classes_in_primary.add(row["rep_class"])

print(f"  Phase 1 (prioritised classes): {len(primary_ids)} origins, {primary_bp:,} bp "
      f"({len(classes_in_primary)} classes)")

# Phase 2: if under budget, add second-best from underrepresented high-priority classes
if primary_bp < TARGET_BP * 0.95:
    remaining = synth_df[~synth_df["locus_id"].isin(primary_ids)].copy()
    remaining = remaining.sort_values("priority_score", ascending=False)
    for _, row in remaining.iterrows():
        if primary_bp >= TARGET_BP:
            break
        primary_ids.append(row["locus_id"])
        primary_bp += row["length_bp"]
    print(f"  Phase 2 (fill to target): {len(primary_ids)} origins, {primary_bp:,} bp")

# Assign tier
df["synthesis_tier"] = df["locus_id"].apply(
    lambda x: "PRIMARY" if x in primary_ids else "SECONDARY"
)
df["tier_rank"] = 0
for i, lid in enumerate(primary_ids):
    df.loc[df["locus_id"] == lid, "tier_rank"] = i + 1

primary = df[df["synthesis_tier"] == "PRIMARY"]
secondary = df[df["synthesis_tier"] == "SECONDARY"]

print(f"\n  ═══ TIER SUMMARY ═══")
print(f"  PRIMARY:   {len(primary)} origins, {primary['length_bp'].sum():,} bp ({primary['length_bp'].sum()/1000:.1f} kb)")
print(f"             Est. cost: ${primary['Twist_est_cost'].sum():,.2f}")
print(f"             Rep classes: {primary['rep_class'].nunique()}")
print(f"             Strains: {primary['strain'].nunique()}")
print(f"             Mean functional score: {primary['functional_score'].mean():.1f}")
print(f"             Synthesis ease: {dict(primary['Twist_ease_cat'].value_counts())}")
print(f"  SECONDARY: {len(secondary)} origins, {secondary['length_bp'].sum():,} bp ({secondary['length_bp'].sum()/1000:.1f} kb)")
print(f"             Est. cost: ${secondary['Twist_est_cost'].sum():,.2f}")


# ═════════════════════════════════════════════════════════════════════════════
# 8.  DIVERSITY STATISTICS SUMMARY
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 8: Diversity statistics ───")

def diversity_stats(subset, label):
    """Print diversity statistics for a subset."""
    print(f"\n  {label} ({len(subset)} origins, {subset['length_bp'].sum():,} bp)")
    print(f"  {'─'*50}")

    # Replicon types
    print(f"  Replicon types: {dict(subset['replicon_type'].value_counts())}")

    # pTi types covered
    pti = subset[subset["pTi_type"].notna()]["pTi_type"].unique()
    print(f"  pTi types ({len(pti)}): {sorted(pti)}")

    # pRi types covered
    pri = subset[subset["pRi_type"].notna()]["pRi_type"].unique()
    if len(pri) > 0:
        print(f"  pRi types ({len(pri)}): {sorted(pri)}")

    # Rep classes
    print(f"  Unique rep classes: {subset['rep_class'].nunique()}")

    # Strains
    print(f"  Strains represented: {subset['strain'].nunique()}")

    # Species
    print(f"  Species: {dict(subset['species'].value_counts())}")

    # Functional completeness
    print(f"  Complete repABC: {subset['repABC_complete'].sum()}/{len(subset)}")
    print(f"  With rep_origin: {subset['has_rep_origin'].sum()}/{len(subset)}")
    print(f"  With parS: {subset['has_parS'].sum()}/{len(subset)}")
    print(f"  With ctRNA: {subset['has_ctRNA'].sum()}/{len(subset)}")
    print(f"  With S-element: {subset['has_S_element'].sum()}/{len(subset)}")

    # Size distribution
    print(f"  Length: {subset['length_bp'].min():,}-{subset['length_bp'].max():,} bp "
          f"(mean {subset['length_bp'].mean():,.0f}, median {subset['length_bp'].median():,.0f})")

    # GC
    print(f"  GC: {subset['GC_content'].min():.3f}-{subset['GC_content'].max():.3f} "
          f"(mean {subset['GC_content'].mean():.3f})")

    # Synthesis
    print(f"  Synthesis ease: {dict(subset['Twist_ease_cat'].value_counts())}")
    print(f"  Homopolymer flags: {subset['Twist_has_homopolymer'].sum()}")
    print(f"  Est. total cost: ${subset['Twist_est_cost'].sum():,.2f}")

diversity_stats(primary, "PRIMARY SYNTHESIS TIER")
diversity_stats(secondary, "SECONDARY SYNTHESIS TIER")
diversity_stats(df, "ALL repABC ORIGINS")


# ═════════════════════════════════════════════════════════════════════════════
# 9.  WRITE OUTPUTS
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 9: Writing outputs ───")

# Sort: primary first (by tier_rank), then secondary by priority_score
df_out = df.sort_values(
    ["synthesis_tier", "tier_rank", "priority_score"],
    ascending=[True, True, False]
).reset_index(drop=True)

# Full metadata CSV
csv_path = OUT_DIR / "repABC_origins_metadata.csv"
df_out.to_csv(csv_path, index=False)
print(f"  Wrote {csv_path}")

# Primary-only CSV
primary_csv = OUT_DIR / "primary_synthesis_origins.csv"
df_out[df_out["synthesis_tier"] == "PRIMARY"].to_csv(primary_csv, index=False)
print(f"  Wrote {primary_csv}")

# Secondary-only CSV
secondary_csv = OUT_DIR / "secondary_synthesis_origins.csv"
df_out[df_out["synthesis_tier"] == "SECONDARY"].to_csv(secondary_csv, index=False)
print(f"  Wrote {secondary_csv}")

# Primary FASTA
primary_fasta = OUT_DIR / "primary_synthesis_origins.fasta"
with open(primary_fasta, "w") as fh:
    for _, row in df_out[df_out["synthesis_tier"] == "PRIMARY"].iterrows():
        lid = row["locus_id"]
        hdr = (f">{lid} | {row['strain']} | {row['replicon_type']}"
               f" | {row.get('pTi_type','') or row.get('pRi_type','') or ''}"
               f" | {row['length_bp']}bp | GC={row['GC_content']:.3f}"
               f" | ease={row['Twist_ease_score']} | tier_rank={int(row['tier_rank'])}")
        fh.write(hdr + "\n")
        seq = sequences[lid]
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")
print(f"  Wrote {primary_fasta}")

# Secondary FASTA
secondary_fasta = OUT_DIR / "secondary_synthesis_origins.fasta"
with open(secondary_fasta, "w") as fh:
    for _, row in df_out[df_out["synthesis_tier"] == "SECONDARY"].iterrows():
        lid = row["locus_id"]
        hdr = (f">{lid} | {row['strain']} | {row['replicon_type']}"
               f" | {row.get('pTi_type','') or row.get('pRi_type','') or ''}"
               f" | {row['length_bp']}bp | GC={row['GC_content']:.3f}"
               f" | ease={row['Twist_ease_score']}")
        fh.write(hdr + "\n")
        seq = sequences[lid]
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")
print(f"  Wrote {secondary_fasta}")

# Summary statistics text
summary_path = OUT_DIR / "diversity_summary.txt"
with open(summary_path, "w") as fh:
    fh.write("repABC Origin Diversity & Synthesis Prioritization Summary\n")
    fh.write("=" * 60 + "\n\n")

    fh.write(f"Total origins analysed: {len(df)}\n")
    fh.write(f"Primary tier:   {len(primary)} origins, {primary['length_bp'].sum():,} bp\n")
    fh.write(f"Secondary tier: {len(secondary)} origins, {secondary['length_bp'].sum():,} bp\n\n")

    fh.write("PRIMARY TIER COMPOSITION\n")
    fh.write("-" * 40 + "\n")
    for _, row in df_out[df_out["synthesis_tier"] == "PRIMARY"].iterrows():
        fh.write(f"  {int(row['tier_rank']):3d}. {row['locus_id']:<55s} "
                 f"{row['length_bp']:5d} bp  GC={row['GC_content']:.3f}  "
                 f"ease={row['Twist_ease_score']}  ${row['Twist_est_cost']:.0f}  "
                 f"func={row['functional_score']}  class={row['rep_class']}\n")

    fh.write(f"\n  Total: {primary['length_bp'].sum():,} bp  "
             f"Est. cost: ${primary['Twist_est_cost'].sum():,.2f}\n")

    fh.write(f"\nSECONDARY TIER ({len(secondary)} origins)\n")
    fh.write("-" * 40 + "\n")
    for _, row in df_out[df_out["synthesis_tier"] == "SECONDARY"].iterrows():
        fh.write(f"  {row['locus_id']:<55s} "
                 f"{row['length_bp']:5d} bp  GC={row['GC_content']:.3f}  "
                 f"ease={row['Twist_ease_score']}  ${row['Twist_est_cost']:.0f}  "
                 f"func={row['functional_score']}  class={row['rep_class']}\n")

print(f"  Wrote {summary_path}")


# ═════════════════════════════════════════════════════════════════════════════
# 10.  GENERATE PLASMID MAP PDFs
# ═════════════════════════════════════════════════════════════════════════════
print("\n─── Step 10: Generating plasmid map PDFs ───")


def detect_at_rich(seq, window=100, threshold=0.65, step=10):
    """Find AT-rich regions."""
    seq = seq.upper()
    hits = []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i+window]
        at = (w.count("A") + w.count("T")) / len(w)
        if at >= threshold:
            hits.append((i, i + window, at))
    # Merge overlapping
    if not hits:
        return []
    merged = [list(hits[0])]
    for s, e, at in hits[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
            merged[-1][2] = max(merged[-1][2], at)
        else:
            merged.append([s, e, at])
    return [(s, e, at) for s, e, at in merged]


def draw_repABC_map(ax, locus_id, row, seq, feats, pars_positions):
    """Draw a linear map of a repABC origin locus."""
    span = len(seq)
    ax.set_xlim(-span * 0.05, span * 1.18)
    ax.set_ylim(-4.0, 6.0)
    ax.axis("off")

    # ── backbone ─────────────────────────────────────────────────────────
    ax.plot([0, span], [0, 0], "k-", lw=3, solid_capstyle="butt")

    # ── scale ticks ──────────────────────────────────────────────────────
    tick = max(500, (span // 10 // 500) * 500)
    if tick == 0:
        tick = 500
    for p in range(0, span + 1, tick):
        ax.plot([p, p], [-0.2, 0.2], "k-", lw=0.8)
        ax.text(p, -0.5, f"{p/1000:.1f}kb", ha="center", fontsize=6, color="gray")

    # ── AT-rich regions ──────────────────────────────────────────────────
    at_rich = detect_at_rich(seq)
    for s, e, at in at_rich:
        ax.axvspan(s, e, ymin=0.42, ymax=0.58, color="#FFF176", alpha=0.8, zorder=1)
        ax.text((s + e) / 2, 0.28, f"AT\n{at*100:.0f}%",
                ha="center", fontsize=5, color="#F9A825")

    # ── parS sites ───────────────────────────────────────────────────────
    if pars_positions:
        for ps, pe in pars_positions:
            mid = (ps + pe) / 2
            ax.plot([mid, mid], [-0.35, 0.35], color="#F06292", lw=1.5, zorder=5)
            ax.plot(mid, 0.4, marker="v", color="#F06292", markersize=4, zorder=5)

    # ── draw feature arrows ──────────────────────────────────────────────
    y_fwd, y_rev, ah = 0.9, -1.2, 0.55

    def _arrow(st, en, strand, color, label):
        y = y_fwd if strand >= 0 else y_rev
        hd = min(span * 0.025, abs(en - st) * 0.25, 300)
        if strand >= 0:
            body = max(st, en - hd)
            ax.fill_between([st, body], y - ah/2, y + ah/2,
                            color=color, alpha=0.75, lw=0.5, edgecolor="k")
            ax.fill([body, en, body], [y + ah/2, y, y - ah/2],
                    color=color, alpha=0.75, lw=0.5, edgecolor="k")
            lx, ly, va = (st + en) / 2, y + ah/2 + 0.18, "bottom"
        else:
            body = min(en, st + hd)
            ax.fill_between([body, en], y - ah/2, y + ah/2,
                            color=color, alpha=0.75, lw=0.5, edgecolor="k")
            ax.fill([body, st, body], [y + ah/2, y, y - ah/2],
                    color=color, alpha=0.75, lw=0.5, edgecolor="k")
            lx, ly, va = (st + en) / 2, y - ah/2 - 0.15, "top"
        short = label[:28] + "..." if len(label) > 28 else label
        ax.text(lx, ly, short, ha="center", va=va, fontsize=6, style="italic",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.7, ec="none"),
                clip_on=True)

    # Draw rep_origin as a region bar (not arrow)
    for f in feats:
        if f["category"] == "rep_origin":
            ax.fill_between([f["start"], f["end"]], -0.15, 0.15,
                            color=FEATURE_COLORS["rep_origin"], alpha=0.6, zorder=2)
            ax.text((f["start"] + f["end"]) / 2, -0.35, "oriV",
                    ha="center", fontsize=5.5, color="#E65100", fontweight="bold")

    # Draw CDS arrows
    for f in feats:
        if f["category"] in ("repA", "repB", "repC", "other_CDS"):
            _arrow(f["start"], f["end"], f["strand"],
                   FEATURE_COLORS[f["category"]], f["label"])

    # Draw regulatory elements as small markers
    for f in feats:
        if f["category"] in ("ctRNA", "S-element"):
            mid = (f["start"] + f["end"]) / 2
            y_pos = 2.0
            color = FEATURE_COLORS[f["category"]]
            ax.annotate(f["label"], xy=(mid, 0.2), xytext=(mid, y_pos),
                        ha="center", fontsize=5, color=color, fontweight="bold",
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))

    # ── GC content track ─────────────────────────────────────────────────
    win = max(50, span // 100)
    gc_p, gc_v = [], []
    for i in range(0, span - win + 1, win // 2):
        sub = seq[i:i+win].upper()
        gc_v.append((sub.count("G") + sub.count("C")) / len(sub))
        gc_p.append(i + win / 2)
    if gc_v:
        ga = np.array(gc_v)
        rng = ga.max() - ga.min() + 1e-6
        yg = -2.5 + (ga - ga.min()) / rng
        ax.fill_between(gc_p, -2.5, yg, color="#90CAF9", alpha=0.45)
        ax.plot(gc_p, yg, color="#1565C0", lw=0.7)
        ax.text(-span * 0.04, -2.2, "GC%", ha="right", fontsize=6, color="#1565C0")
        ax.text(-span * 0.04, -2.6, f"{ga.mean()*100:.1f}%",
                ha="right", fontsize=5, color="#1565C0")

    # ── title + metadata ─────────────────────────────────────────────────
    rep_type = row.get("pTi_type") or row.get("pRi_type") or row["replicon_type"]
    ax.set_title(
        f"{locus_id.replace('_repABC','')}  |  {rep_type}  |  {row['species']}",
        fontsize=9, fontweight="bold", pad=4
    )

    tier = row.get("synthesis_tier", "?")
    tier_color = "#2E7D32" if tier == "PRIMARY" else "#1565C0"
    rank_str = f"  rank #{int(row['tier_rank'])}" if row.get("tier_rank", 0) > 0 else ""

    completeness = []
    if row["has_repA"]:
        s = "repA"
        if row.get("repA_pseudogene"):
            s += "(pseudo)"
        completeness.append(s)
    if row["has_repB"]: completeness.append("repB")
    if row["has_repC"]: completeness.append("repC")
    comp_str = "+".join(completeness) if completeness else "none"

    elements = []
    if row["has_rep_origin"]: elements.append("oriV")
    if row.get("has_parS"):   elements.append(f"parS({int(row['n_parS'])})")
    if row["has_ctRNA"]:      elements.append("ctRNA")
    if row["has_S_element"]:  elements.append("S-elem")
    elem_str = ", ".join(elements) if elements else "none"

    meta = (
        f"Span: {span:,} bp  |  GC: {row['GC_content']*100:.1f}%  |  "
        f"Tier: {tier}{rank_str}\n"
        f"Genes: {comp_str}  |  Elements: {elem_str}  |  "
        f"Ease: {row['Twist_ease_score']}/100 ({row['Twist_ease_cat']})  |  "
        f"Cost: ${row['Twist_est_cost']:.0f}"
    )
    ax.text(span * 0.5, 5.2, meta, ha="center", va="top", fontsize=6, color=tier_color,
            bbox=dict(boxstyle="round", fc="#F1F8E9", alpha=0.85, ec=tier_color, lw=1.2))

    # ── legend ───────────────────────────────────────────────────────────
    hdls = [mpatches.Patch(color=c, alpha=0.75, label=k.replace("_", " "))
            for k, c in FEATURE_COLORS.items()]
    hdls.append(mpatches.Patch(color="#FFF176", alpha=0.8, label="AT-rich"))
    ax.legend(handles=hdls, loc="upper right", fontsize=5, ncol=2,
              framealpha=0.8, edgecolor="gray")


def generate_repABC_pdf(df_subset, label, path):
    """Generate a multi-page PDF with one map per origin."""
    n = len(df_subset)
    print(f"  Generating {label} PDF ({n} maps)...")

    # Parse parS positions per locus from the metadata
    pars_pos_map = {}
    for _, row in df_subset.iterrows():
        lid = row["locus_id"]
        positions = []
        pos_str = row.get("parS_positions", "")
        if isinstance(pos_str, str) and pos_str.strip():
            for chunk in pos_str.split(";"):
                chunk = chunk.strip()
                if "-" in chunk:
                    parts = chunk.split("-")
                    try:
                        positions.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
        pars_pos_map[lid] = positions

    with PdfPages(str(path)) as pdf:
        # Title page
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.text(0.5, 0.65, f"repABC Origin Maps — {label}",
                ha="center", fontsize=24, fontweight="bold", transform=ax.transAxes)
        ax.text(0.5, 0.50,
                f"{n} origins  |  {df_subset['length_bp'].sum():,} bp total",
                ha="center", fontsize=16, transform=ax.transAxes)
        tier_stats = (
            f"Synthesis ease: {dict(df_subset['Twist_ease_cat'].value_counts())}  |  "
            f"Est. cost: ${df_subset['Twist_est_cost'].sum():,.0f}"
        )
        ax.text(0.5, 0.38, tier_stats,
                ha="center", fontsize=11, transform=ax.transAxes, color="#1B5E20")
        rep_classes = f"Rep classes: {df_subset['rep_class'].nunique()}  |  Strains: {df_subset['strain'].nunique()}"
        ax.text(0.5, 0.28, rep_classes,
                ha="center", fontsize=11, transform=ax.transAxes, color="#1565C0")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for idx, (_, row) in enumerate(df_subset.iterrows()):
            lid = row["locus_id"]
            seq = sequences.get(lid, "")
            if not seq:
                continue
            feats = gbk_features.get(lid, [])
            pars_pos = pars_pos_map.get(lid, [])

            try:
                fig, ax = plt.subplots(figsize=(14, 5))
                draw_repABC_map(ax, lid, row, seq, feats, pars_pos)
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
            except Exception as e:
                fig, ax = plt.subplots(figsize=(14, 5))
                ax.text(0.5, 0.5, f"{lid}: {e}",
                        ha="center", va="center", transform=ax.transAxes, color="red")
                pdf.savefig(fig, bbox_inches="tight")
            finally:
                plt.close(fig)

            if (idx + 1) % 20 == 0:
                print(f"    Maps: {idx+1}/{n}")

    print(f"  Wrote {path}")


# Generate PDFs
primary_pdf_path = OUT_DIR / "primary_repABC_maps.pdf"
generate_repABC_pdf(
    df_out[df_out["synthesis_tier"] == "PRIMARY"],
    "PRIMARY", primary_pdf_path
)

secondary_pdf_path = OUT_DIR / "secondary_repABC_maps.pdf"
generate_repABC_pdf(
    df_out[df_out["synthesis_tier"] == "SECONDARY"],
    "SECONDARY", secondary_pdf_path
)

# Also generate a combined all-origins PDF
all_pdf_path = OUT_DIR / "all_repABC_maps.pdf"
generate_repABC_pdf(df_out, "ALL ORIGINS", all_pdf_path)


print("\n═══ DONE ═══")
print(f"All outputs in: {OUT_DIR}/")
