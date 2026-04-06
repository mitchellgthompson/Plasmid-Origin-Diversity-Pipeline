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
| `plasmid_origin_pipeline.py` | Main pipeline (runs steps S1–S12) |
| `origins_to_genbank.py` | Pfam annotation, origin trimming, GenBank export |
| `RIPs/RIP.hmm` | Rep-protein HMM profile (35 families) |
| `RIPs/RIP.fasta` | Source sequences for the RIP HMM |
| `RIPs/Collected_RIPs.xlsx` | Curated metadata for the RIP families |
| `results_plsdb/` | Example output from a recent run |
| `results_plsdb/genbanks/` | One annotated GenBank per origin |
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

## Outputs

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
