"""Build `final_order_maps.pdf` by extracting the already-rendered per-origin
pages from the two upstream pipelines and concatenating them in final-order
sequence. No re-rendering, no barcode overlay — these are the exact maps
produced by `analyze_repABC_origins.py` (repABC) and `plasmid_origin_pipeline.py`
(diversity) for the 73 origins in `final_order.csv`.

Page layout:
  1.  Title / summary page
  2...  One page per origin, in this order:
        - repABC primary  (26 origins, Twist)
        - repABC alt_cloning  (7 origins, genomic PCR)
        - diversity  (40 origins, Twist), sorted by selection tier

Source PDFs (must already exist):
  - results_repABC/primary_repABC_maps.pdf  (title + 33 repABC primary maps)
  - results_plsdb/origin_plasmid_maps.pdf  (title + 112 diversity maps)

Barcodes + cassette stay in the CSV / FASTA only — they are not drawn on
these structural maps.
"""
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "final_order"
OUT_DIR.mkdir(exist_ok=True)
IN_CSV = OUT_DIR / "final_order.csv"
OUT_PDF = OUT_DIR / "final_order_maps.pdf"

REPABC_CSV = ROOT / "results_repABC" / "primary_synthesis_origins.csv"
REPABC_PDF = ROOT / "results_repABC" / "primary_repABC_maps.pdf"
DIVERSITY_META = ROOT / "results_plsdb" / "diversity_library_metadata.csv"
DIVERSITY_PDF = ROOT / "results_plsdb" / "origin_plasmid_maps.pdf"


def make_title_pdf(rows, path):
    """Render a single-page summary PDF and return its path."""
    origin_bp = sum(len(r["sequence"]) for r in rows)
    cassette_bp = sum(len(r["barcoded_sequence"]) - len(r["sequence"]) for r in rows)
    order_bp = origin_bp + cassette_bp
    by_method = Counter(r["production_method"] for r in rows)
    by_pool = Counter(r["source_tier"] for r in rows)
    by_tier = Counter(r["priority_or_category"] for r in rows if r["source_tier"] == "diversity_PLSDB")

    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.text(0.5, 0.92, "Final Synthesis / Cloning Order",
            ha="center", fontsize=22, weight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.86,
            f"{len(rows)} origins  |  {origin_bp:,} bp of origin + {cassette_bp:,} bp of barcode cassette  |  {order_bp:,} bp total",
            ha="center", fontsize=12, transform=ax.transAxes)

    summary_lines = [
        f"Origin bp (plotted on the following pages): {origin_bp:>8,}",
        f"Cassette bp (57 bp per origin, NOT shown):  {cassette_bp:>8,}",
        f"Total order bp:                             {order_bp:>8,}",
        "",
        "Pool breakdown:",
        *[f"  {k:<30}{v:>3} origins" for k, v in by_pool.items()],
        "",
        "Production method:",
        *[f"  {k:<30}{v:>3} origins" for k, v in by_method.items()],
        "",
        "Diversity-tier mix:",
        *[f"  {k:<34}{v:>3}" for k, v in by_tier.items()],
        "",
        "Barcode cassette (top-strand 5'→3', appended to 3' end of each origin):",
        "  A  |  GATGTCCACGAGGTCTCT  |  [N20]  |  CGTACGCTGCAGGTCGAC",
        "       (spacer)        (U1)         (barcode)        (U2)",
        "",
        "Barcodes are in final_order.csv (barcode_N20, barcoded_sequence) and",
        "final_order_barcoded.fasta. See barcode_diversity_report.{txt,pdf}",
        "for the Hamming/Levenshtein/positional QC audit.",
        "",
        "The maps on the following pages are the structural annotations from:",
        "  - analyze_repABC_origins.py  (repABC, from ARC_repABC_loci_fna_gbk)",
        "  - plasmid_origin_pipeline.py (diversity, with Pfam-annotated ORFs)",
    ]
    ax.text(0.06, 0.78, "\n".join(summary_lines), ha="left", va="top",
            transform=ax.transAxes, fontsize=9, family="monospace")

    with PdfPages(str(path)) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def build_page_maps():
    """Return {locus_id → PageObject} for both sources."""
    # repABC: page 0 = title, pages 1..N match primary_synthesis_origins.csv row order
    repabc_rows = list(csv.DictReader(REPABC_CSV.open()))
    repabc_pdf = PdfReader(str(REPABC_PDF))
    assert len(repabc_pdf.pages) == len(repabc_rows) + 1, \
        f"repABC PDF pages {len(repabc_pdf.pages)} != {len(repabc_rows)} + 1"
    repabc_map = {r["locus_id"]: repabc_pdf.pages[i + 1] for i, r in enumerate(repabc_rows)}

    # Diversity: page 0 = title, pages 1..N match diversity_library_metadata.csv row order
    div_rows = list(csv.DictReader(DIVERSITY_META.open()))
    div_pdf = PdfReader(str(DIVERSITY_PDF))
    assert len(div_pdf.pages) == len(div_rows) + 1, \
        f"diversity PDF pages {len(div_pdf.pages)} != {len(div_rows)} + 1"
    div_map = {r["NUCCORE_UID"]: div_pdf.pages[i + 1] for i, r in enumerate(div_rows)}

    return repabc_map, div_map


TIER_ORDER = {
    "Tier 3: Cross-phylum BHR": 0,
    "Tier 2: Alphaproteobacteria": 1,
    "Tier 4: Phylum-exclusive": 2,
    "Tier 5: Diversity fill": 3,
}
POOL_ORDER = {"repABC_primary": 0, "repABC_alt_cloning": 1, "diversity_PLSDB": 2}


def main():
    rows = list(csv.DictReader(IN_CSV.open()))
    repabc_map, div_map = build_page_maps()

    # Order final_order origins for deterministic page layout
    sorted_rows = sorted(
        rows,
        key=lambda r: (
            POOL_ORDER.get(r["source_tier"], 9),
            TIER_ORDER.get(r["priority_or_category"], 9),
            r["id"],
        ),
    )

    # Render the summary title page to a temp PDF
    title_pdf_path = OUT_DIR / "_final_order_title.pdf"
    make_title_pdf(rows, title_pdf_path)

    # Concatenate: title + per-origin pages
    writer = PdfWriter()
    title_reader = PdfReader(str(title_pdf_path))
    writer.add_page(title_reader.pages[0])

    missing = []
    added = 0
    for r in sorted_rows:
        src_page = None
        if r["source_tier"] in ("repABC_primary", "repABC_alt_cloning"):
            src_page = repabc_map.get(r["id"])
        elif r["source_tier"] == "diversity_PLSDB":
            src_page = div_map.get(r["id"])
        if src_page is None:
            missing.append((r["id"], r["source_tier"]))
            continue
        writer.add_page(src_page)
        added += 1

    with OUT_PDF.open("wb") as fh:
        writer.write(fh)
    title_pdf_path.unlink(missing_ok=True)

    print(f"Wrote {OUT_PDF.name}: {added + 1} pages ({added} origins + 1 title)")
    if missing:
        print(f"WARNING: {len(missing)} origin(s) had no page in the source PDFs:")
        for oid, pool in missing:
            print(f"  {oid}  ({pool})")


if __name__ == "__main__":
    main()
