# Plasmid Origin Diversity Pipeline

End-to-end pipeline that produces a synthesis-ready, barcoded library of
bacterial plasmid replication origins (oriV) for parallel profiling via
RB-TnSeq BarSeq. Origins come from two sources — a curated repABC set from
Agrobacterium/Rhizobium, and a broad-host-range diversity library mined from
**PLSDB 2025** — and are merged, constraint-filtered, deduplicated, and
tagged with unique 20 bp barcodes flanked by the canonical Deutschbauer-lab
BarSeq priming sites (Wetmore et al. 2015).

## What it does

The pipeline runs in four stages, each with its own script:

**A. Diversity library discovery from PLSDB 2025** — `plasmid_origin_pipeline.py`
   → `origins_to_genbank.py`

1. Load PLSDB 2025 metadata (72k plasmids) and filter to circular, bacterial,
   non-`repABC`, non-narrow-host-range *E. coli* plasmids ≤ 6 kb.
2. Stream the PLSDB FASTA in batches: predict ORFs with `pyrodigal`, call Rep
   proteins with `hmmsearch` vs `RIPs/RIP.hmm`, and score intergenic
   sequences with the embedded OriVFinder algorithm (iteron detection,
   Z'-curve AT-rich regions, DnaA box / Fis / IHF pattern matching).
3. Define origin windows that contain the Rep gene + best-scoring IGS + ORF
   buffer; deduplicate by Rep HMM type; select a tiered diversity library
   (Alphaproteobacteria → cross-phylum BHR → phylum-exclusive → fill).
4. Annotate ORFs against `Pfam-A`, trim each origin to its functional core,
   and emit FASTA, CSV, PDF maps, overview PNG, and one GenBank per origin.

**B. repABC origin analysis from curated Agrobacterium/Rhizobium loci** —
   `analyze_repABC_origins.py`

5. Parse 118 `.fna` + `.gbk` pairs in `ARC_repABC_loci_fna_gbk/`; extract
   repA/B/C genes, oriV, ctRNA, S-element, UniRef50 cross-references;
   integrate parS partition-site annotations from `parS.table.xlsx`.
6. Score each locus for Twist synthesizability and functional completeness;
   rank into a **primary** tier (33 diverse, synthesizable representatives)
   and a **secondary** tier (85 backups).
7. Write metadata CSVs, FASTAs, and per-origin plasmid-map PDFs.

**C. Merge + constraint-filter + budget fit + RB-TnSeq barcoding** —
   `build_combined_synthesis_csv.py` → `build_final_order.py`

8. Merge repABC primary + PLSDB diversity into `combined_synthesis_origins.csv`
   (145 origins).
9. Apply Ian Blaby's Twist constraint-fail list (2026-04-20): drop 33
   diversity fails; re-tag 7 repABC fails (all pTi-family) as
   `production_method=genomic_pcr_amplification` — still in the library,
   just produced from genomic templates instead of Twist.
10. Exact-sequence dedup; exclude `MobT`, `Phg_2220_C`, `Gemini_AL1__RCR`
    Rep families; priority-greedy-fit to a **250 kb total-order budget**
    that accounts for the 57 bp barcode cassette appended to every origin.
11. Append a cassette to the 3′ end of each origin (top-strand 5′→3′):
    **`A — U1(GATGTCCACGAGGTCTCT) — [N20] — U2(CGTACGCTGCAGGTCGAC)`** where
    U1/U2 are the canonical Wetmore 2015 BarSeq priming sites
    (verified by direct match against pKMW7 at `results_plsdb/JBx_109896.gb`)
    and each N20 is deterministic, Hamming ≥ 3, GC 40–60 %, no homopolymer
    ≥ 4, absent both strands in every target origin.
12. Render the barcoded library as a one-page-per-origin linear PDF
    (`make_final_order_pdf.py`).

**D. Barcode diversity audit** — `analyze_barcode_diversity.py`

13. Audit the N20 set for BarSeq demultiplexing robustness: pairwise Hamming
    and Levenshtein distance distributions, positional base-composition
    entropy, GC content, revcomp collisions, shared-*k*-mer overlap, and
    comparison against a random 20-mer baseline.

Final deliverables for the collaborator live in **[`final_order/`](final_order)**.

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
| `pLANT_master_ORI_template.gb` | **Destination binary vector** (9,056 bp) into which every origin in the final order will be cloned — Gibson-assembly slot flanked by FseI / AscI cut sites |
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

# 5. Analyse the curated repABC locus set (no arguments — inputs in ARC_repABC_loci_fna_gbk/)
python3 analyze_repABC_origins.py

# 6. Merge repABC + PLSDB diversity into one reference CSV
python3 build_combined_synthesis_csv.py

# 7. Build the barcoded, budget-fit, synthesis-ready order
python3 build_final_order.py          # -> final_order/final_order.csv
python3 make_final_order_pdf.py       # -> final_order/final_order_maps.pdf
python3 analyze_barcode_diversity.py  # -> final_order/barcode_diversity_report.{txt,pdf}
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

### Destination vector

Every origin in the final order will be cloned into [`pLANT_master_ORI_template.gb`](pLANT_master_ORI_template.gb)
— a 9,056 bp circular binary vector ("Master Binary Vector"). The vector
carries the plant-transformation payload (T-DNA LB/RB, Rhodo GFP, plant
pNOS::eGFP, NOS terminator, C58 Overdrive, dPCR QC targets, Kanamycin
selection, ColE1 ori for *E. coli* propagation, conjugal-transfer origin
from pEX18-Gm) and a single dropout slot for the replication origin:

```
… Gibson RB Super Site 2 (30 bp) | FseI | ORI dropout (20 bp stuffer) | AscI | Gibson ORI (30 bp) …
                              ^cut^                                 ^cut^
   5' homology arm                                                 3' homology arm
```

Linearise the vector with FseI + AscI to release the 20 bp stuffer, then
Gibson-assemble each barcoded origin in. The 30 bp flanking sequences are
not included in the ordered fragments — they are part of the vector
backbone and provide the Gibson homology directly.

### Final deliverables

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

`plasmid_origin_pipeline.py` flags:

| Flag | Default | Purpose |
|------|---------|---------|
| `--plsdb-dir` | (required) | PLSDB 2025 directory containing `nuccore.csv`, `taxonomy.csv`, `typing.csv`, `sequences.fasta` |
| `--rip-hmm` | `RIPs/RIP.hmm` | Rep-protein HMM profile |
| `--email` | (required) | Contact email (NCBI requirement; only used if NCBI fallback is triggered) |
| `--output-dir` | `results_plsdb/` | Output directory |
| `--min-origins` | 100 | Minimum origins in the diversity library |
| `--target-bp` | 500_000 | Target total bp (the pipeline overshoots ~1.5× to compensate for downstream trimming) |
| `--hmm-evalue` | 1e-5 | HMM hit E-value cutoff |

`analyze_repABC_origins.py`, `build_combined_synthesis_csv.py`,
`build_final_order.py`, `make_final_order_pdf.py`, and
`analyze_barcode_diversity.py` take **no CLI arguments** — configurable
knobs (budget, excluded Rep families, constraint-fail lists, RNG seed,
cassette structure) are pinned as constants at the top of each script so
the library build is reproducible.

## Pipeline architecture

```
 Stage A — PLSDB diversity discovery
┌────────────────────────────────────────────────────────────────────┐
│  plasmid_origin_pipeline.py                                        │
│    S1  Load PLSDB metadata (72k → 47k filtered)                    │
│    S2  Filter (circular, ≤6 kb, not repABC, not E. coli NHR)       │
│    S3  Stream FASTA, predict ORFs with pyrodigal (parallel)        │
│    S4  hmmsearch vs RIPs/RIP.hmm                                   │
│    S5  OriVFinder IGS scoring (iterons, Z'-curve, DnaA/Fis/IHF)    │
│    S6  Origin windowing (≥200 bp ORF buffer, ≤6 kb)                │
│    S7  Local ORF-based trimming                                    │
│    S8  Deduplicate by RIP type (≤5 per type)                       │
│    S9  Tiered diversity selection (Alpha / BHR / phylum / fill)    │
│    S10 Functional validation + ORF labels                          │
│    S11 Twist synthesis assessment                                  │
│    S12 Write CSV, FASTA, PDF maps, overview PNG                    │
│                         │                                          │
│                         ▼                                          │
│  origins_to_genbank.py                                             │
│    Pfam-A / RIP annotation → per-origin GenBank + origin trimming  │
└─────────────────────────────┬──────────────────────────────────────┘
                              │  results_plsdb/diversity_library_metadata.csv
                              │
 Stage B — repABC analysis    │
┌─────────────────────────────┼──────────────────────────────────────┐
│  analyze_repABC_origins.py  │                                      │
│    Parse 118 Agro/Rhizo .fna+.gbk pairs                            │
│    Integrate parS partition sites (parS.table.xlsx)                │
│    Rank + tier (primary / secondary) for Twist synthesis           │
│    Write CSV, FASTA, primary/secondary/all PDF maps                │
└─────────────────────────────┬──────────────────────────────────────┘
                              │  results_repABC/primary_synthesis_origins.{csv,fasta}
                              │
 Stage C — merge + barcoded order
┌─────────────────────────────┼──────────────────────────────────────┐
│  build_combined_synthesis_csv.py                                   │
│    Merge repABC primary (33) + PLSDB diversity (112) -> reference  │
│                         │                                          │
│                         ▼   combined_synthesis_origins.csv         │
│  build_final_order.py                                              │
│    Drop Ian's Twist constraint-fails (re-tag 7 repABC as alt_clone)│
│    Exact-sequence dedup                                            │
│    Exclude MobT / Phg_2220_C / Gemini_AL1 Rep families             │
│    Priority greedy fit to 250 kb (effective size = origin + 57 bp) │
│    Append cassette:  A | U1 | N20 | U2  (Wetmore 2015 BarSeq)      │
│    Flag post-cassette oversize / junction-homopolymer issues       │
│                         │                                          │
│                         ▼                                          │
│  make_final_order_pdf.py  → one page per origin, cassette shown    │
└─────────────────────────────┬──────────────────────────────────────┘
                              │  final_order/final_order.csv
                              │
 Stage D — barcode diversity audit
┌─────────────────────────────┼──────────────────────────────────────┐
│  analyze_barcode_diversity.py                                      │
│    Pairwise Hamming + Levenshtein distance distributions           │
│    Positional base composition + Shannon entropy                   │
│    GC distribution, revcomp collisions, shared k-mer overlap       │
│    PASS/FAIL verdict vs random-20-mer baseline                     │
└────────────────────────────────────────────────────────────────────┘
                              │  final_order/barcode_diversity_report.{txt,pdf}
                              ▼
                    Deliver `final_order/` to collaborator
```

## Citing this work

Please cite the upstream tools and databases:

- **PLSDB 2025**: <https://ccb-microbe.cs.uni-saarland.de/plsdb2025/>
- **OriVFinder**: doi:10.1093/nar/gkaf341
- **Pfam**: <https://pfam.xfam.org/>
- **HMMER**: <http://hmmer.org/>
- **pyrodigal**: doi:10.21105/joss.04296
- **Biopython**: doi:10.1093/bioinformatics/btp163
