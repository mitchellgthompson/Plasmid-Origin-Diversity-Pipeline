"""Render `final_order.csv` as an annotated linear-map PDF, one page per origin.

Features (CDS arrows, ncRNA boxes, parS / oriV / iteron / AT-rich markers) are
pulled from the per-origin GenBank files produced upstream:

  - repABC origins  -> ARC_repABC_loci_fna_gbk/<id>.gbk  (exact-length match)
  - diversity (PLSDB) -> results_plsdb/genbanks/<Accession>.gb  (substring-
    aligned to the CSV origin sequence; if alignment fails the script falls
    back to pyrodigal ORF prediction on the CSV sequence itself)

The 3' end of every fragment is the 57 bp barcode cassette:
    A  |  U1 (GATGTCCACGAGGTCTCT)  |  N20  |  U2 (CGTACGCTGCAGGTCGAC)
"""
import csv
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "final_order"
OUT_DIR.mkdir(exist_ok=True)
IN_CSV = OUT_DIR / "final_order.csv"
OUT_PDF = OUT_DIR / "final_order_maps.pdf"

REPABC_GBK_DIR = ROOT / "ARC_repABC_loci_fna_gbk"
PLSDB_GB_DIR = ROOT / "results_plsdb" / "genbanks"
DIVERSITY_META = ROOT / "results_plsdb" / "diversity_library_metadata.csv"

SPACER = "A"
U1 = "GATGTCCACGAGGTCTCT"
U2 = "CGTACGCTGCAGGTCGAC"

TIER_COLOR = {
    "repABC_primary": "#1565C0",
    "repABC_alt_cloning": "#C62828",
    "diversity_PLSDB": "#6A1B9A",
}
U1_COLOR = "#FFB74D"
BC_COLOR = "#FFF176"
U2_COLOR = "#4FC3F7"
SPACER_COLOR = "#E0E0E0"

# CDS coloring by product / function keyword
FEATURE_COLORS = {
    "rep":         "#E53935",  # RepA/B/C and Rep_init families
    "partition":   "#43A047",  # parA/parB, partition
    "mobilization":"#FB8C00",  # mob
    "replication": "#1E88E5",  # helicase, topoisomerase
    "recombinase": "#8E24AA",
    "resistance":  "#6D4C41",  # AMR
    "structural":  "#546E7A",
    "hypothetical":"#BDBDBD",
    "other":       "#9E9E9E",
    "ncRNA":       "#00897B",
    "oriv":        "#FFA000",
    "pars":        "#D81B60",
    "s_element":   "#5E35B1",
    "iteron":      "#F06292",
    "at_rich":     "#FFF59D",
}

REP_KEYWORDS = re.compile(r"\b(rep[abc]?|rep_init|replicat|replication protein|repA|repB|repC)\b", re.I)
PAR_KEYWORDS = re.compile(r"\b(par[abc]?|partition|chromosome-partitioning|parS)\b", re.I)
MOB_KEYWORDS = re.compile(r"\b(mob|relaxase|mobilization|conjugat|traI|traG)\b", re.I)
RES_KEYWORDS = re.compile(r"\b(antibiotic|resistance|beta-lactam|aminoglycoside|tetracycline)\b", re.I)
REPL_KEYWORDS = re.compile(r"\b(helicase|topoisomerase|primase|polymerase)\b", re.I)
RECO_KEYWORDS = re.compile(r"\b(recombinase|integrase|resolvase|serine|tyrosine)\b", re.I)
HYP_KEYWORDS = re.compile(r"\b(hypothetical|unknown|putative)\b", re.I)


def classify_cds(product):
    p = (product or "").lower()
    if REP_KEYWORDS.search(p): return "rep"
    if PAR_KEYWORDS.search(p): return "partition"
    if MOB_KEYWORDS.search(p): return "mobilization"
    if REPL_KEYWORDS.search(p): return "replication"
    if RECO_KEYWORDS.search(p): return "recombinase"
    if RES_KEYWORDS.search(p): return "resistance"
    if HYP_KEYWORDS.search(p): return "hypothetical"
    return "other"


def classify_misc(label):
    l = (label or "").lower()
    if "ori" in l and "v" in l: return "oriv"
    if "pars" in l:             return "pars"
    if "s_element" in l or "s-element" in l: return "s_element"
    if "iteron" in l:           return "iteron"
    if "at" in l and "rich" in l: return "at_rich"
    if "ncrna" in l or "ctrna" in l: return "ncRNA"
    return None


def load_features_from_gb(path, target_seq):
    """Parse GenBank features and map them onto target_seq.

    Returns a list of (start, end, strand, category, label) tuples using 0-based
    coordinates on target_seq, or None if no alignment is possible.
    """
    try:
        rec = SeqIO.read(path, "genbank")
    except Exception:
        return None

    gb_seq = str(rec.seq).upper()
    t = target_seq.upper()

    # Determine coordinate offset between GenBank seq and target_seq.
    offset = 0; truncate_start = 0; truncate_end = len(t)
    if gb_seq == t:
        offset = 0
    elif gb_seq in t:
        offset = t.index(gb_seq)
    elif t in gb_seq:
        offset = -gb_seq.index(t)
        truncate_end = len(t)
    else:
        return None

    feats = []
    for f in rec.features:
        try:
            s = int(f.location.start); e = int(f.location.end)
        except (TypeError, ValueError):
            continue
        s += offset; e += offset
        if e < 0 or s >= len(t):
            continue
        s = max(0, s); e = min(len(t), e)
        strand = 1 if (getattr(f.location, "strand", 1) or 1) >= 0 else -1
        if f.type == "CDS":
            product = f.qualifiers.get("product", [""])[0]
            cat = classify_cds(product)
            feats.append((s, e, strand, cat, product[:40]))
        elif f.type == "ncRNA":
            product = f.qualifiers.get("product", [""])[0]
            feats.append((s, e, strand, "ncRNA", product[:40] or "ncRNA"))
        elif f.type in ("misc_feature", "regulatory", "rep_origin"):
            label = (f.qualifiers.get("label", [""])[0]
                     or f.qualifiers.get("note", [""])[0]
                     or f.type)
            cat = classify_misc(label) or "other"
            feats.append((s, e, strand, cat, label[:40]))
    return feats


def predict_features_pyrodigal(seq):
    """Fallback ORF prediction when no GenBank alignment is available."""
    try:
        import pyrodigal
    except ImportError:
        return []
    gf = pyrodigal.GeneFinder(meta=True)
    genes = gf.find_genes(seq.encode())
    feats = []
    for g in genes:
        strand = 1 if g.strand > 0 else -1
        feats.append((g.begin - 1, g.end, strand, "other", "predicted ORF"))
    return feats


def features_for_row(row, acc_lookup):
    target = row["sequence"].upper()
    if row["source_tier"] == "diversity_PLSDB":
        acc = acc_lookup.get(row["id"])
        if acc:
            gb = PLSDB_GB_DIR / f"{acc}.gb"
            if gb.exists():
                feats = load_features_from_gb(gb, target)
                if feats is not None:
                    return feats, f".gb ({acc})"
        return predict_features_pyrodigal(target), "pyrodigal fallback"
    else:
        gbk = REPABC_GBK_DIR / f"{row['id']}.gbk"
        if gbk.exists():
            feats = load_features_from_gb(gbk, target)
            if feats is not None:
                return feats, ".gbk"
        return [], "no_annotation"


def draw_arrow(ax, s, e, strand, y0, h, color, ec="black"):
    w = e - s
    if w <= 0:
        return
    head = min(w * 0.25, max(w * 0.12, 40))
    if strand >= 0:
        pts = [(s, y0), (e - head, y0), (e - head, y0 - h * 0.15),
               (e, y0 + h / 2), (e - head, y0 + h + h * 0.15),
               (e - head, y0 + h), (s, y0 + h)]
    else:
        pts = [(e, y0), (s + head, y0), (s + head, y0 - h * 0.15),
               (s, y0 + h / 2), (s + head, y0 + h + h * 0.15),
               (s + head, y0 + h), (e, y0 + h)]
    ax.add_patch(mpatches.Polygon(pts, closed=True, facecolor=color,
                                   edgecolor=ec, lw=0.6, alpha=0.9))


def draw_page(ax_map, ax_meta, row, feats, source):
    origin_len = len(row["sequence"])
    cassette_len = 1 + len(U1) + len(row["barcode_N20"]) + len(U2)
    total = origin_len + cassette_len

    # Backbone line
    ax_map.axhline(0.5, xmin=0.01, xmax=0.99, color="#424242", lw=1.0, zorder=0)
    ax_map.add_patch(mpatches.Rectangle((0, 0.35), origin_len, 0.3,
                                         facecolor="#ECEFF1", edgecolor="#90A4AE", lw=0.5, zorder=1))

    # Features (CDS arrows + ncRNA/misc boxes)
    for s, e, strand, cat, label in feats:
        color = FEATURE_COLORS.get(cat, FEATURE_COLORS["other"])
        if cat in ("ncRNA", "pars", "at_rich", "iteron", "s_element", "oriv"):
            y0, h = (0.72, 0.15) if strand >= 0 else (0.13, 0.15)
            ax_map.add_patch(mpatches.Rectangle((s, y0), e - s, h,
                                                 facecolor=color, edgecolor="black",
                                                 lw=0.4, alpha=0.85, zorder=2))
            if e - s > total * 0.03:
                ax_map.text((s + e) / 2, y0 + h / 2, label, ha="center",
                            va="center", fontsize=5, color="white", weight="bold", zorder=3)
        else:
            y0 = 0.55 if strand >= 0 else 0.30
            draw_arrow(ax_map, s, e, strand, y0, 0.18, color)
            if cat == "rep" and e - s > total * 0.02:
                ax_map.text((s + e) / 2, y0 + 0.09, label, ha="center",
                            va="center", fontsize=5, color="white", weight="bold", zorder=3)

    # Cassette at 3' end
    x = origin_len
    ax_map.add_patch(mpatches.Rectangle((x, 0.38), 1, 0.24, facecolor=SPACER_COLOR,
                                         edgecolor="black", lw=0.4, zorder=2))
    x += 1
    ax_map.add_patch(mpatches.Rectangle((x, 0.38), len(U1), 0.24, facecolor=U1_COLOR,
                                         edgecolor="black", lw=0.4, zorder=2))
    ax_map.text(x + len(U1) / 2, 0.68, "U1", ha="center", fontsize=6, color="#6D4C41")
    x += len(U1)
    ax_map.add_patch(mpatches.Rectangle((x, 0.38), 20, 0.24, facecolor=BC_COLOR,
                                         edgecolor="black", lw=0.4, zorder=2))
    ax_map.text(x + 10, 0.68, "N20", ha="center", fontsize=6, color="#5D4037")
    x += 20
    ax_map.add_patch(mpatches.Rectangle((x, 0.38), len(U2), 0.24, facecolor=U2_COLOR,
                                         edgecolor="black", lw=0.4, zorder=2))
    ax_map.text(x + len(U2) / 2, 0.68, "U2", ha="center", fontsize=6, color="#01579B")

    # Axis
    ax_map.set_xlim(-total * 0.005, total * 1.005)
    ax_map.set_ylim(0, 1)
    ax_map.set_yticks([])
    ax_map.set_xticks([0, origin_len, total])
    ax_map.set_xticklabels(["0", f"{origin_len:,}", f"{total:,}"], fontsize=7)
    ax_map.spines[:].set_visible(False)
    ax_map.set_xlabel("bp (5'→3', top strand)", fontsize=7)

    title = (f"{row['id']}  |  {row['species_or_strain']}  |  "
             f"{row['source_tier']}  |  {row['production_method']}")
    ax_map.set_title(title, fontsize=9, loc="left", weight="bold", pad=4)

    # CDS count + source note
    cds_count = sum(1 for _, _, _, cat, _ in feats
                    if cat not in ("ncRNA", "pars", "at_rich", "iteron", "s_element", "oriv"))
    other_count = len(feats) - cds_count
    ax_map.text(total, 0.98, f"annotations: {cds_count} CDS + {other_count} misc  ({source})",
                ha="right", va="top", fontsize=6, color="#455A64")

    # Metadata block on the second axis
    ax_meta.axis("off")
    flags = []
    if row.get("oversize_flag"):
        flags.append("OVERSIZE >5 kb (Ian to re-review)")
    if row.get("junction_homopolymer_flag"):
        flags.append(f"JUNCTION: {row['junction_homopolymer_flag']}")
    flag_str = "  ".join(flags) if flags else "none"

    meta_text = (
        f"Tier / category:     {row['priority_or_category']}\n"
        f"Rep family:          {row['rep_class_or_family']}\n"
        f"Taxonomy:            {row['taxonomy_class']} / {row['taxonomy_family']}\n"
        f"GC content:          {row['GC_content']}\n"
        f"Barcode (N20):       {row['barcode_N20']}\n"
        f"Cassette (3' tail):  A | {U1} | N20 | {U2}\n"
        f"Ordered fragment:    {len(row['barcoded_sequence']):,} bp (origin {origin_len:,} + cassette 57)\n"
        f"Flags:               {flag_str}"
    )
    ax_meta.text(0.02, 0.95, meta_text, ha="left", va="top", transform=ax_meta.transAxes,
                  fontsize=7, family="monospace",
                  bbox=dict(boxstyle="round,pad=0.4", facecolor="#F5F5F5", edgecolor="#BDBDBD"))

    # Legend
    legend_handles = [mpatches.Patch(facecolor=FEATURE_COLORS[k], edgecolor="black", label=k)
                      for k in ("rep", "partition", "mobilization", "replication",
                                "recombinase", "resistance", "hypothetical", "other",
                                "ncRNA", "oriv", "pars", "s_element", "iteron", "at_rich")]
    ax_meta.legend(handles=legend_handles, loc="center right", ncol=2, fontsize=6,
                    framealpha=0.9, title="feature categories", title_fontsize=7)


def main():
    rows = list(csv.DictReader(IN_CSV.open()))
    acc_lookup = {m["NUCCORE_UID"]: m["Accession"]
                   for m in csv.DictReader(DIVERSITY_META.open())}
    n = len(rows)

    origin_bp = sum(len(r["sequence"]) for r in rows)
    cassette_bp = sum(len(r["barcoded_sequence"]) - len(r["sequence"]) for r in rows)
    order_bp = origin_bp + cassette_bp
    by_method = Counter(r["production_method"] for r in rows)
    by_pool = Counter(r["source_tier"] for r in rows)
    by_tier = Counter(r["priority_or_category"] for r in rows if r["source_tier"] == "diversity_PLSDB")

    source_counts = Counter()

    with PdfPages(str(OUT_PDF)) as pdf:
        # --- Title page ---
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.text(0.5, 0.92, "Final Synthesis / Cloning Order",
                ha="center", fontsize=22, weight="bold", transform=ax.transAxes)
        ax.text(0.5, 0.86, f"{n} origins  |  {order_bp:,} bp total  |  250 kb budget",
                ha="center", fontsize=13, transform=ax.transAxes)

        summary = [
            f"Origin bp:                  {origin_bp:>7,}",
            f"Cassette bp (57 bp each):   {cassette_bp:>7,}",
            f"Total order bp:             {order_bp:>7,}",
            "",
            "Pool breakdown:",
            *[f"  {k:<28}{v:>3} origins" for k, v in by_pool.items()],
            "",
            "Production method:",
            *[f"  {k:<28}{v:>3} origins" for k, v in by_method.items()],
            "",
            "Diversity-tier mix:",
            *[f"  {k:<34}{v:>3}" for k, v in by_tier.items()],
            "",
            "Barcode cassette (top strand 5'->3'):",
            f"  {SPACER} | {U1} | N20 | {U2}    (spacer | U1 | barcode | U2)",
            "",
            "Priming sites verified against pKMW7 (Wetmore et al. 2015, mBio).",
        ]
        ax.text(0.08, 0.76, "\n".join(summary), ha="left", va="top",
                transform=ax.transAxes, fontsize=9, family="monospace")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- One page per origin ---
        order_sort = {"repABC_primary": 0, "repABC_alt_cloning": 1, "diversity_PLSDB": 2}
        sorted_rows = sorted(rows, key=lambda r: (order_sort.get(r["source_tier"], 9), r["id"]))
        for i, row in enumerate(sorted_rows):
            feats, source = features_for_row(row, acc_lookup)
            source_counts[source] += 1
            fig, (ax_map, ax_meta) = plt.subplots(2, 1, figsize=(11, 6.5),
                                                   gridspec_kw={"height_ratios": [1.6, 1]})
            draw_page(ax_map, ax_meta, row, feats, source)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            if (i + 1) % 20 == 0:
                print(f"  rendered {i + 1}/{n}")

    print(f"Wrote {OUT_PDF}")
    print(f"Annotation sources used: {dict(source_counts)}")


if __name__ == "__main__":
    main()
