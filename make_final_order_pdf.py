"""Render `final_order.csv` as a PDF, one linear-fragment map per origin.

Each page shows the origin sequence as a bar with the 3'-appended RB-TnSeq
barcode cassette (U1 | N20 | U2) highlighted at the tail, plus the metadata
needed for ordering and BarSeq amplification.
"""
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).parent
IN_CSV = ROOT / "final_order.csv"
OUT_PDF = ROOT / "final_order_maps.pdf"

SPACER = "A"
U1 = "GATGTCCACGAGGTCTCT"
U2 = "CGTACGCTGCAGGTCGAC"

METHOD_COLOR = {
    "twist_synthesis": "#2E7D32",
    "genomic_pcr_amplification": "#C62828",
}
TIER_COLOR = {
    "repABC_primary": "#1565C0",
    "repABC_alt_cloning": "#C62828",
    "diversity_PLSDB": "#6A1B9A",
}
U1_COLOR = "#FFB74D"
BC_COLOR = "#FFF176"
U2_COLOR = "#4FC3F7"


def draw_fragment(ax, row):
    origin_len = len(row["sequence"])
    cassette_len = len(SPACER) + len(U1) + len(row["barcode_N20"]) + len(U2)
    total = origin_len + cassette_len
    height = 1.0
    y0 = 0.0

    # Origin body
    tier = row["source_tier"]
    ax.add_patch(mpatches.Rectangle((0, y0), origin_len, height,
                                     facecolor=TIER_COLOR.get(tier, "#555"),
                                     edgecolor="black", lw=0.8, alpha=0.85))
    ax.text(origin_len / 2, y0 + height / 2,
            f"{row['id']}  ·  {origin_len:,} bp",
            ha="center", va="center", fontsize=8, color="white", weight="bold")

    # Cassette: spacer | U1 | N20 | U2
    x = origin_len
    ax.add_patch(mpatches.Rectangle((x, y0), len(SPACER), height,
                                     facecolor="#E0E0E0", edgecolor="black", lw=0.6))
    x += len(SPACER)
    ax.add_patch(mpatches.Rectangle((x, y0), len(U1), height,
                                     facecolor=U1_COLOR, edgecolor="black", lw=0.6))
    ax.text(x + len(U1) / 2, y0 + height + 0.15, "U1", ha="center", fontsize=6, color="#6D4C41")
    x += len(U1)
    ax.add_patch(mpatches.Rectangle((x, y0), len(row["barcode_N20"]), height,
                                     facecolor=BC_COLOR, edgecolor="black", lw=0.6))
    ax.text(x + len(row["barcode_N20"]) / 2, y0 + height + 0.15,
            f"N20", ha="center", fontsize=6, color="#5D4037")
    x += len(row["barcode_N20"])
    ax.add_patch(mpatches.Rectangle((x, y0), len(U2), height,
                                     facecolor=U2_COLOR, edgecolor="black", lw=0.6))
    ax.text(x + len(U2) / 2, y0 + height + 0.15, "U2", ha="center", fontsize=6, color="#01579B")

    # Axis / scale
    ax.set_xlim(-total * 0.01, total * 1.01)
    ax.set_ylim(-1.2, 1.8)
    ax.set_aspect("auto")
    ax.spines[:].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([0, origin_len, total])
    ax.set_xticklabels([0, f"{origin_len:,}", f"{total:,}"], fontsize=7)
    ax.set_xlabel("bp (5'→3', top strand)", fontsize=7)

    # Metadata block
    flags = []
    if row.get("oversize_flag"):
        flags.append("OVERSIZE (>5 kb after cassette — Ian to review)")
    if row.get("junction_homopolymer_flag"):
        flags.append(f"JUNCTION_HOMOPOLYMER: {row['junction_homopolymer_flag']}")
    flag_line = "  ".join(flags) if flags else "none"

    meta = (
        f"Production method:   {row['production_method']}\n"
        f"Pool:                {tier}\n"
        f"Tier / rank:         {row['priority_or_category']}\n"
        f"Rep family:          {row['rep_class_or_family']}\n"
        f"Species / strain:    {row['species_or_strain']}\n"
        f"Taxonomy:            {row['taxonomy_class']} / {row['taxonomy_family']}\n"
        f"GC content:          {row['GC_content']}\n"
        f"\n"
        f"Barcode (N20):       {row['barcode_N20']}\n"
        f"Cassette (3' tail):  {SPACER} | {U1} | N20 | {U2}\n"
        f"                     (spacer | BarSeq_P2 site | barcode | BarSeq_P1 site)\n"
        f"Ordered fragment:    {len(row['barcoded_sequence']):,} bp (origin {origin_len:,} + cassette 57)\n"
        f"Flags:               {flag_line}"
    )
    ax.text(0, -0.25, meta, ha="left", va="top", transform=ax.transData,
            fontsize=7, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F5F5F5", edgecolor="#BDBDBD"))


def main():
    rows = list(csv.DictReader(IN_CSV.open()))
    n = len(rows)
    print(f"Rendering {n} origins -> {OUT_PDF.name}")

    origin_bp = sum(len(r["sequence"]) for r in rows)
    cassette_bp = sum(len(r["barcoded_sequence"]) - len(r["sequence"]) for r in rows)
    order_bp = origin_bp + cassette_bp
    by_method = Counter(r["production_method"] for r in rows)
    by_pool = Counter(r["source_tier"] for r in rows)
    by_tier = Counter(r["priority_or_category"] for r in rows if r["source_tier"] == "diversity_PLSDB")

    with PdfPages(str(OUT_PDF)) as pdf:
        # --- Title / summary page ---
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.text(0.5, 0.90, "Final Synthesis / Cloning Order",
                ha="center", fontsize=22, weight="bold", transform=ax.transAxes)
        ax.text(0.5, 0.84, f"{n} origins  |  {order_bp:,} bp total  |  250 kb budget",
                ha="center", fontsize=13, transform=ax.transAxes)

        summary = [
            f"Origin bp:                  {origin_bp:>7,}",
            f"Cassette bp (56 bp each):   {cassette_bp:>7,}",
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
            "Barcode cassette, top-strand 5'->3' (appended to 3' end of each origin):",
            f"  1 bp spacer '{SPACER}'     (breaks any G-homopolymer extension into U1)",
            f"  U1 ({len(U1)} bp) = {U1}   (BarSeq_P2 priming site)",
            f"  N20 = 20 bp unique per-origin barcode (Hamming >= 3, GC 40-60%, no homopolymer >= 4)",
            f"  U2 ({len(U2)} bp) = {U2}   (BarSeq_P1 priming site)",
            "",
            "Priming sites verified by direct match against pKMW7 (Wetmore et al. 2015, mBio).",
        ]
        ax.text(0.08, 0.74, "\n".join(summary), ha="left", va="top",
                transform=ax.transAxes, fontsize=9, family="monospace")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- One page per origin ---
        order_sort = {"repABC_primary": 0, "repABC_alt_cloning": 1, "diversity_PLSDB": 2}
        sorted_rows = sorted(rows, key=lambda r: (order_sort.get(r["source_tier"], 9), r["id"]))
        for i, row in enumerate(sorted_rows):
            fig, ax = plt.subplots(figsize=(11, 6))
            draw_fragment(ax, row)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            if (i + 1) % 20 == 0:
                print(f"  rendered {i + 1}/{n}")

    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
