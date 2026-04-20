"""Combine all repABC (primary + secondary) and PLSDB diversity origins into a single CSV with metadata and sequence."""
import csv
from pathlib import Path

ROOT = Path(__file__).parent
REPABC_PRIMARY_CSV = ROOT / "results_repABC" / "primary_synthesis_origins.csv"
REPABC_PRIMARY_FASTA = ROOT / "results_repABC" / "primary_synthesis_origins.fasta"
REPABC_SECONDARY_CSV = ROOT / "results_repABC" / "secondary_synthesis_origins.csv"
REPABC_SECONDARY_FASTA = ROOT / "results_repABC" / "secondary_synthesis_origins.fasta"
DIV_RANKED_CSV = ROOT / "results_plsdb" / "synthesis_ranked_origins.csv"
DIV_META_CSV = ROOT / "results_plsdb" / "diversity_library_metadata.csv"
OUT_CSV = ROOT / "combined_synthesis_origins.csv"


def parse_fasta(path):
    seqs, cur_id, cur_buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur_id is not None:
                seqs[cur_id] = "".join(cur_buf)
            cur_id = line[1:].split(" ")[0].split("|")[0].strip()
            cur_buf = []
        else:
            cur_buf.append(line.strip())
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_buf)
    return seqs


# 1. repABC primary + secondary
repabc_primary_seqs = parse_fasta(REPABC_PRIMARY_FASTA)
repabc_primary_rows = list(csv.DictReader(REPABC_PRIMARY_CSV.open()))
repabc_secondary_seqs = parse_fasta(REPABC_SECONDARY_FASTA)
repabc_secondary_rows = list(csv.DictReader(REPABC_SECONDARY_CSV.open()))

# 2. PLSDB diversity: full library metadata (sequence lives in Origin_sequence)
div_meta_rows = list(csv.DictReader(DIV_META_CSV.open()))
# Also load the ranked-for-synthesis subset so we can flag/score those 10 as PRIMARY
div_ranked_rows = list(csv.DictReader(DIV_RANKED_CSV.open()))
div_ranked_info = {}
for r in div_ranked_rows:
    parts = r["name"].split(" ")
    acc = parts[0]
    kv = {p.split("=", 1)[0]: p.split("=", 1)[1] for p in parts[1:] if "=" in p}
    div_ranked_info[acc] = kv


# Unified columns
COLS = [
    "source_tier",
    "id",
    "species_or_strain",
    "length_bp",
    "GC_content",
    "rep_class_or_family",
    "functional_or_synthesis_score",
    "priority_or_category",
    "est_cost_usd",
    "ease_score",
    "taxonomy_class",
    "taxonomy_family",
    "notes",
    "sequence",
]

out_rows = []

def add_repabc_rows(rows, seqs, tier_label):
    for r in rows:
        locus = r["locus_id"]
        seq = seqs.get(locus, "")
        out_rows.append({
            "source_tier": tier_label,
            "id": locus,
            "species_or_strain": r.get("species") or r.get("strain", ""),
            "length_bp": r.get("length_bp", ""),
            "GC_content": r.get("GC_content", ""),
            "rep_class_or_family": r.get("rep_class", ""),
            "functional_or_synthesis_score": r.get("functional_score", ""),
            "priority_or_category": f"{r.get('synthesis_tier','')} (rank {r.get('tier_rank','')})",
            "est_cost_usd": r.get("Twist_est_cost", ""),
            "ease_score": r.get("Twist_ease_score", ""),
            "taxonomy_class": r.get("class", ""),
            "taxonomy_family": r.get("family", ""),
            "notes": f"repABC_complete={r.get('repABC_complete','')}; has_parS={r.get('has_parS','')}; pTi_type={r.get('pTi_type','')}; pRi_type={r.get('pRi_type','')}",
            "sequence": seq,
        })


add_repabc_rows(repabc_primary_rows, repabc_primary_seqs, "repABC_primary")

# --- PLSDB diversity (full library, 112 origins across 4 selection tiers) ---
for meta in div_meta_rows:
    acc = meta.get("NUCCORE_UID", "")
    seq = meta.get("Origin_sequence", "")
    ranked = div_ranked_info.get(acc, {})
    species = meta.get("Species", "") or meta.get("TAXONOMY_taxon_name", "")
    out_rows.append({
        "source_tier": "diversity_PLSDB",
        "id": acc,
        "species_or_strain": species.replace("_", " "),
        "length_bp": len(seq) or meta.get("Length_bp", ""),
        "GC_content": meta.get("GC_content", ""),
        "rep_class_or_family": meta.get("Origin_RIP_type", "") or meta.get("Rep_hmm_name", ""),
        "functional_or_synthesis_score": ranked.get("synthesis_score", meta.get("oriv_score", "")),
        "priority_or_category": ranked.get("category", "") or meta.get("Selection_tier", ""),
        "est_cost_usd": "",
        "ease_score": meta.get("Twist_ease_score", ""),
        "taxonomy_class": meta.get("Class", ""),
        "taxonomy_family": meta.get("Family", ""),
        "notes": f"OriVFinder_total_score={meta.get('OriVFinder_total_score','')}; Selection_tier={meta.get('Selection_tier','')}; Replicon_types={meta.get('Replicon_types','')}; synthesis_ranked={'yes' if ranked else 'no'}",
        "sequence": seq,
    })

with OUT_CSV.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    w.writerows(out_rows)

from collections import Counter
tier_counts = Counter(r["source_tier"] for r in out_rows)
total_bp = sum(len(r["sequence"]) for r in out_rows)
missing = [r["id"] for r in out_rows if not r["sequence"]]
print(f"Wrote {OUT_CSV.name}: {len(out_rows)} rows, {total_bp:,} bp total")
for t, n in tier_counts.items():
    print(f"  {t}: {n}")
if missing:
    print(f"WARNING: missing sequence for: {missing}")
