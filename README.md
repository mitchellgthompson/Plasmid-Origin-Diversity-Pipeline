# Plasmid Origin Diversity Pipeline

End-to-end pipeline for identifying, characterising, and annotating plasmid
replication origins (oriV) from the **PLSDB 2025** database. Designed to
produce a small, diverse, broad-host-range library of origins suitable for
synthetic biology applications.

## What it does

1. **Loads PLSDB 2025** metadata (72,556 plasmids) and filters to circular,
   bacterial, non-`repABC`, non-narrow-host-range *E. coli* plasmids ≤ 6 kb.
2. **Streams** the PLSDB sequence FASTA in 2,000-sequence batches:
   - Predicts ORFs with `pyrodigal` (parallelised across CPU cores).
   - Identifies Rep proteins by `hmmsearch` against a curated RIP HMM set.
   - Scores intergenic sequences flanking each Rep gene with an embedded
     OriVFinder algorithm (iteron detection, Z'-curve AT-rich regions,
     DnaA box / Fis / IHF pattern matching).
3. **Defines origin windows** that contain the Rep gene, the best-scoring
   IGS, and a ≥ 200 bp buffer around any ORF inside the window.
4. **Deduplicates** by Rep HMM type (keeps up to 5 representatives per type).
5. **Selects a diversity library** in tiers (Agrobacterium-native →
   Alphaproteobacteria → cross-phylum BHR → phylum-exclusive → fill).
6. **Validates** each origin (Rep + Par + structural element scoring).
7. **Annotates ORFs** against `RIP.hmm` and `Pfam-A`, then **trims** each
   origin to its functional core (removing phage / metabolism / AMR ORFs
   that were padding the 6 kb window).
8. **Writes outputs** as a FASTA, CSV, plasmid-map PDF, overview PNG, and
   one annotated GenBank file per origin.

## Repository contents

| File / dir | Purpose |
|------------|---------|
| `plasmid_origin_pipeline.py` | Main PLSDB pipeline (runs steps S1–S12) |
| `origins_to_genbank.py` | Pfam annotation, origin trimming, GenBank export |
| `analyze_repABC_origins.py` | **repABC origin analysis, synthesis ranking & PDF maps** |
| `build_combined_synthesis_csv.py` | Merge repABC + diversity outputs into one reference CSV |
| `build_final_order.py` | **Build the final ≤250 kb barcoded order** (constraint filter → dedup → family exclusions → priority fit → RB-TnSeq cassette) |
| `make_final_order_pdf.py` | Render `final_order.csv` as a linear-fragment PDF (one page per origin) |
| `combined_synthesis_origins.csv` | Pre-filter reference CSV (all 230 repABC + diversity origins) |
| `final_order/` | **All final deliverables** (CSV, PDFs, barcode diversity report) |
| `final_order/final_order.csv` | Final order — 73 barcoded origins (≈249 kb) with `production_method`, `barcode_N20`, `barcoded_sequence` columns |
| `final_order/final_order_maps.pdf` | One-page-per-origin linear map of the final order, cassette highlighted |
| `final_order/barcode_diversity_report.txt` | Barcode diversity audit (Hamming, Levenshtein, positional balance) |
| `final_order/barcode_diversity_report.pdf` | Histograms + positional-composition heatmap |
| `analyze_barcode_diversity.py` | Script that produces the barcode-diversity audit outputs |
| `ARC_repABC_loci_fna_gbk/` | 118 curated repABC loci (.fna + .gbk) from Agrobacterium/Rhizobium |
| `ARC_repABC_loci_fna_gbk/parS.table.xlsx` | parS partition-site annotations for repABC loci |
| `results_repABC/` | repABC pipeline outputs (metadata, FASTAs, PDFs, summary) |
| `RIPs/RIP.hmm` | Rep-protein HMM profile (35 families) |
| `RIPs/RIP.fasta` | Source sequences for the RIP HMM |
| `RIPs/Collected_RIPs.xlsx` | Curated metadata for the RIP families |
| `results_plsdb/` | Example output from the PLSDB pipeline |
| `results_plsdb/genbanks/` | One annotated GenBank per origin |
| `results_plsdb/JBx_109896.gb` | pKMW7 reference (RB-TnSeq transposon vector, Wetmore 2015) used to verify U1/U2 priming sites |
| `requirements.txt` | Python dependencies |

## External data (not in repo)

These files are too large to commit and must be downloaded separately:

- **PLSDB 2025**: <https://ccb-microbe.cs.uni-saarland.de/plsdb2025/>
  Place the relational tables and `sequences.fasta` in `plsdb2025/`.
- **Pfam-A HMM**: <https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz>
  Decompress and run `hmmpress Pfam-A.hmm` once to build the binary index.

## Dependencies

```bash
pip install -r requirements.txt
```

You also need **HMMER 3** on `PATH` (`hmmsearch`, `hmmscan`, `hmmpress`):

```bash
conda install -c bioconda hmmer
# or: brew install hmmer
```

## Quick start

```bash
# 1. Download PLSDB 2025 (one-off)
mkdir -p plsdb2025
cd plsdb2025
wget https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download/plsdb_2025_03_28.zip
unzip plsdb_2025_03_28.zip
bunzip2 sequences.fasta.bz2
cd ..

# 2. Download and index Pfam-A (one-off, ~2 GB)
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

# 3. Run the pipeline (S1–S12, ~2 h on a laptop)
python plasmid_origin_pipeline.py \
    --plsdb-dir   plsdb2025/ \
    --rip-hmm     RIPs/RIP.hmm \
    --email       you@example.com \
    --output-dir  results_plsdb/ \
    --min-origins 75 \
    --target-bp   250000

# 4. Annotate ORFs with Pfam-A and trim non-origin flanking ORFs
python origins_to_genbank.py \
    --fasta       results_plsdb/origin_sequences.fasta \
    --csv         results_plsdb/diversity_library_metadata.csv \
    --rip-hmm     RIPs/RIP.hmm \
    --pfam-hmm    Pfam-A.hmm \
    --plsdb-fasta plsdb2025/sequences.fasta \
    --outdir      results_plsdb/genbanks \
    --buffer      200
```

## repABC Origin Analysis Pipeline

A separate pipeline for characterising **repABC-family replication origins**
from curated Agrobacterium/Rhizobium loci. It parses annotated GenBank files,
integrates parS partition-site data, runs Twist synthesis feasibility
assessment, computes diversity statistics, and prioritises origins into
**primary** (~125 kb, synthesisable and diverse) and **secondary** (backup)
synthesis tiers.

### Quick start

```bash
python3 analyze_repABC_origins.py
```

No arguments needed — input data lives in `ARC_repABC_loci_fna_gbk/`.

### repABC pipeline steps

1. **Parse** 118 `.fna`/`.gbk` file pairs — extract repA/B/C genes, oriV,
   ctRNA, S-element annotations, UniRef50 cross-references
2. **Integrate parS** data from `parS.table.xlsx` (247 parS hits across 106 loci)
3. **Assign taxonomy** (species, family, class) from strain identifiers
4. **Compute sequence complexity** — homopolymers, direct repeats, local GC
   variation, dinucleotide entropy
5. **Twist synthesis assessment** — GC category, strategy, fragment count,
   estimated cost, ease score (0–100), homopolymer flag
6. **Diversity grouping** by RepC UniRef50 cluster, pTi/pRi type, and
   replicon type (40 unique rep classes)
7. **Composite ranking** (functional completeness 35%, synthesis ease 30%,
   sequence complexity 20%, compactness 15%) and tier assignment:
   - **Primary**: best synthesisable representative per class, budget-capped
     at ~125 kb, prioritising pTi → pRi → named replicons → other classes
   - **Secondary**: all remaining origins as backup
8. **Generate annotated linear plasmid maps** (one PDF page per origin) showing
   repA/B/C arrows, oriV, parS sites, ctRNA, S-element, AT-rich regions,
   GC% track, and metadata box

### repABC outputs

After running, `results_repABC/` contains:

| File | Contents |
|------|----------|
| `repABC_origins_metadata.csv` | Full metadata for all 118 origins |
| `primary_synthesis_origins.csv` | Primary tier metadata |
| `secondary_synthesis_origins.csv` | Secondary tier metadata |
| `primary_synthesis_origins.fasta` | Primary tier sequences |
| `secondary_synthesis_origins.fasta` | Secondary tier sequences |
| `primary_repABC_maps.pdf` | Annotated linear maps for primary tier |
| `secondary_repABC_maps.pdf` | Annotated linear maps for secondary tier |
| `all_repABC_maps.pdf` | Maps for all 118 origins |
| `diversity_summary.txt` | Human-readable ranking and statistics |

---

## Final order assembly (`build_final_order.py`)

Once the upstream repABC and PLSDB pipelines have produced their outputs,
`build_final_order.py` merges them into a single **barcoded, synthesis-ready
order** under a 250 kb total-order budget. The recipe:

1. Start from `combined_synthesis_origins.csv` (33 repABC primary + 112
   diversity origins, 145 total).
2. **Drop Twist constraint-fails** from Ian Blaby's review (2026-04-20):
   7 repABC (all pTi-family) + 33 diversity origins. The 7 repABC do **not**
   leave the library — they are tagged `production_method=genomic_pcr_amplification`
   and will be amplified from genomic/plasmid templates instead.
3. **Exact-sequence deduplicate** the diversity set (drops 4 redundant copies).
4. **Exclude** `MobT`, `Phg_2220_C`, and `Gemini_AL1__RCR` Rep-family origins.
5. **Priority-ordered greedy fit** across diversity tiers
   (Tier 3 cross-phylum BHR > Tier 2 Alphaproteobacteria > Tier 4
   phylum-exclusive > Tier 5 fill) using *effective* fragment size
   (origin + 57 bp cassette) so the total order lands at ≤ 250 kb.
6. **Append an RB-TnSeq barcode cassette** to the 3′ end of every origin,
   top-strand 5′→3′:
   ```
   origin ... — A — GATGTCCACGAGGTCTCT — [N20] — CGTACGCTGCAGGTCGAC
                |   └─ U1 (BarSeq_P2) ─┘   └─ barcode ─┘   └─ U2 (BarSeq_P1) ─┘
                └─ 1 bp spacer (prevents G-homopolymer extension into U1)
   ```
   Priming sites are the canonical Wetmore *et al.* 2015 (mBio) BarSeq U1/U2
   sequences, verified by direct match against the **pKMW7** transposon
   vector (`results_plsdb/JBx_109896.gb`). Each 20 bp barcode satisfies:
   GC 40–60 %, no homopolymer ≥ 4, pairwise Hamming ≥ 3, no overlap with
   U1/U2 or their reverse complements, and absent (both strands) from every
   target origin. Barcodes are regenerated deterministically from a fixed
   RNG seed, so reruns produce identical assignments.
7. **Flag post-cassette issues** for human review:
   - `oversize_flag=YES` if the barcoded fragment exceeds 5,000 bp (Twist
     clonal-gene standard tier ceiling).
   - `junction_homopolymer_flag` if the cassette extends an origin-terminal
     homopolymer across the junction into a ≥ 6 bp G/C or ≥ 10 bp A/T run.

Current result: **73 origins, 248,974 bp total order** — 66 twist_synthesis
(26 repABC + 40 diversity) + 7 genomic_pcr_amplification; 6 origins carry
`oversize_flag=YES` for Ian to confirm with Twist's extended clonal-gene tier.

All final deliverables live in [`final_order/`](final_order):

```
final_order/
├─ final_order.csv                  # 73-origin order, one row per origin, includes
│                                     barcode_N20 and barcoded_sequence columns
├─ final_order_maps.pdf             # one page per origin: linear fragment +
│                                     highlighted cassette + metadata
├─ barcode_diversity_report.txt     # Hamming / Levenshtein / positional audit
└─ barcode_diversity_report.pdf     # histograms + positional-composition heatmap
```

Rebuild end-to-end from the reference CSV:

```bash
python3 build_final_order.py          # -> final_order/final_order.csv
python3 make_final_order_pdf.py       # -> final_order/final_order_maps.pdf
python3 analyze_barcode_diversity.py  # -> final_order/barcode_diversity_report.{txt,pdf}
```

### Barcode diversity audit

`analyze_barcode_diversity.py` audits the 73 N20 barcodes for BarSeq
demultiplexing robustness. The current set passes every check:

| Check | Target | Observed |
|-------|--------|----------|
| min pairwise Hamming distance | ≥ 3 | **7** |
| mean pairwise Hamming vs random baseline | ≈ baseline | 14.97 vs 14.93 |
| min pairwise Levenshtein (edit) distance | ≥ 3 | **6** |
| positional Shannon entropy (bits, max 2.00) | ≥ 1.7 | median 1.970 |
| positional max-base frequency | < 40 % | max 39.7 % |
| GC content | [40 %, 60 %] | all within range, mean 49.9 % |
| revcomp / palindrome collisions | 0 | 0 |

---

## PLSDB Pipeline Outputs

After both scripts complete, `results_plsdb/` contains:

- `diversity_library_metadata.csv` — full metadata for every selected origin
- `origin_sequences.fasta` — trimmed origin sequences
- `origin_plasmid_maps.pdf` — one annotated map per origin (ORFs coloured
  by function)
- `diversity_library_overview.png` — summary charts (classification, GC,
  phylum coverage, tier composition, synthesis ease)
- `genbanks/<Accession>.gb` — one fully-annotated GenBank file per origin,
  ready for opening in SnapGene / Benchling / ApE / etc.

Each GenBank CDS feature carries:
- `/product=` Pfam description (or RIP HMM name for Rep proteins)
- `/function=` category (`rep`, `partition`, `mobilization`, `recombinase`,
  `amr`, `structural`, `hypothetical`, `other`)
- `/db_xref="Pfam:..."` accession when applicable
- `/note="origin_relevant=yes"` if the ORF is plausibly origin-related
- `/ApEinfo_fwdcolor`, `/ApEinfo_revcolor` ApE/SnapGene-compatible colours

## Pipeline parameters

| Flag | Default | Purpose |
|------|---------|---------|
| `--plsdb-dir` | (required) | PLSDB 2025 directory containing `nuccore.csv`, `taxonomy.csv`, `typing.csv`, `sequences.fasta` |
| `--rip-hmm` | `RIPs/RIP.hmm` | Rep-protein HMM profile |
| `--email` | (required) | Contact email (NCBI requirement; only used if NCBI fallback is triggered) |
| `--output-dir` | `results_plsdb/` | Output directory |
| `--min-origins` | 100 | Minimum origins in the diversity library |
| `--target-bp` | 500_000 | Target total bp (the pipeline overshoots ~1.5× to compensate for downstream trimming) |
| `--hmm-evalue` | 1e-5 | HMM hit E-value cutoff |

## Pipeline architecture

```
┌──────────────────────────────────────────────────────────┐
│  plasmid_origin_pipeline.py                              │
│  ┌────────────────────────────────────────────────────┐  │
│  │ S1  Load PLSDB metadata (72k → 47k)                │  │
│  │ S2  Filter (circular, size, not repABC, not NHR)   │  │
│  │ S3  Stream FASTA, predict ORFs (parallel)          │  │
│  │ S4  hmmsearch vs RIP.hmm                           │  │
│  │ S5  OriVFinder IGS scoring                         │  │
│  │ S6  Origin window (≥200 bp ORF buffer, ≤6 kb)      │  │
│  │ S7  Local ORF-based trimming                       │  │
│  │ S8  Dedup by RIP type (5/type)                     │  │
│  │ S9  Tiered diversity selection (1.5× overshoot)    │  │
│  │ S10 Functional OriV validation + ORF labels        │  │
│  │ S11 Twist synthesis assessment                     │  │
│  │ S12 Outputs (CSV, FASTA, PDF, PNG)                 │  │
│  └────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────┐
│  origins_to_genbank.py                                   │
│  ┌────────────────────────────────────────────────────┐  │
│  │ Predict ORFs with pyrodigal                        │  │
│  │ Run hmmsearch vs RIP.hmm  (Rep proteins)           │  │
│  │ Run hmmscan  vs Pfam-A    (everything else)        │  │
│  │ Trim origins to origin-relevant core               │  │
│  │ Iteratively extend to enforce 200 bp ORF buffer    │  │
│  │ Re-predict & re-classify on extended windows       │  │
│  │ Write per-origin .gb files + regenerate PDF/FASTA  │  │
│  └────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────┘
```

## Citing this work

Please cite the upstream tools and databases:

- **PLSDB 2025**: <https://ccb-microbe.cs.uni-saarland.de/plsdb2025/>
- **OriVFinder**: doi:10.1093/nar/gkaf341
- **Pfam**: <https://pfam.xfam.org/>
- **HMMER**: <http://hmmer.org/>
- **pyrodigal**: doi:10.21105/joss.04296
- **Biopython**: doi:10.1093/bioinformatics/btp163
