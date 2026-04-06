#!/usr/bin/env python3
"""
plasmid_origin_pipeline.py

Builds a diversity library of functional plasmid origins of replication
starting directly from the full PLSDB (or a PLSDB-format metadata file),
with integrated OriV finding, filtering, origin trimming, and validation.

PIPELINE STAGES
───────────────
  1. Load PLSDB metadata and normalise columns
  2. Metadata filters: circular/complete/bacterial, not repABC, not NHR E. coli
  3. OriV finding:
       batch ORF prediction (pyrodigal) + Rep-protein HMM search (hmmsearch)
       on the PLSDB sequence FASTA (or NCBI-fetched sequences as fallback)
  4. Origin window definition per plasmid
       · centre on the detected Rep gene
       · enforce ≥200 bp buffer between Rep gene edges and window edges
       · hard cap at 6 000 bp
  5. Window trimming: remove non-origin CDS (mobilisation, AMR, transposons)
       from each edge while preserving the 200 bp Rep-gene buffer zone
  6. Deduplicate by RIP type and build filtered origins table
  7. Tiered diversity selection (~100 origins, ~500 kb)
  8. Fetch GenBank annotations for selected origins (NCBI Entrez)
  9. Functional OriV validation
       · Rep protein  — HMM + GenBank CDS
       · Par system   — GenBank CDS + Walker-A motif scan
       · Structure    — iterons (spacing-aware), AT-rich region, DnaA boxes
 10. Twist synthesis assessment (real-sequence homopolymer check included)
 11. Outputs: diversity_library_metadata.csv, origin_sequences.fasta,
              origin_plasmid_maps.pdf, diversity_library_overview.png

USAGE
─────
  # Recommended – provide PLSDB FASTA for fast local OriV finding
  python plasmid_origin_pipeline.py \\
      --plsdb-meta  plsdb.tsv \\
      --plsdb-fasta plsdb.fna \\
      --rip-hmm     RIPs/RIP.hmm \\
      --email       your@email.com

  # Slower fallback – batch-fetch sequences from NCBI
  python plasmid_origin_pipeline.py \\
      --plsdb-meta  plsdb.tsv \\
      --rip-hmm     RIPs/RIP.hmm \\
      --email       your@email.com

  # Optional: provide NCBI API key to raise rate limit to 10 req/s
      --ncbi-api-key <key>

DEPENDENCIES
────────────
  pip install biopython pyrodigal matplotlib pandas numpy
  hmmsearch must be on PATH  (HMMER 3.x, from http://hmmer.org)
"""

import argparse
import concurrent.futures
import csv
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

from Bio import Entrez, SeqIO
from Bio.Seq import Seq

import pyrodigal


def _pyrodigal_worker(args_tuple):
    """Module-level worker for parallel ORF prediction (must be picklable)."""
    acc, seq_str = args_tuple
    finder = pyrodigal.GeneFinder(meta=True)
    genes = finder.find_genes(seq_str.encode())
    orfs = []
    for k, g in enumerate(genes):
        aa = g.translate().rstrip("*")
        if len(aa) >= 30:
            orfs.append({"id": f"{acc}_orf{k+1}", "start": g.begin,
                         "end": g.end,
                         "strand": "+" if g.strand == 1 else "-", "aa": aa})
    return acc, orfs


try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import argrelextrema
    _SCIPY_OK = True
except ImportError:
    _SCIPY_OK = False

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS AND KEYWORD SETS
# ─────────────────────────────────────────────────────────────────────────────

AGRO_GC         = 0.59
ORIGIN_MAX_BP   = 6_000
REP_BUFFER_BP   = 200       # minimum bp between Rep gene edge and window edge
NCBI_SLEEP      = 0.4       # seconds between NCBI requests (no API key)
NCBI_SLEEP_KEY  = 0.12      # seconds with API key

# ── repABC detection ─────────────────────────────────────────────────────────
# repABC is a tripartite replication system common in Alphaproteobacteria.
# Detected by: (1) PLSDB/MOB-suite rep_type field, (2) CDS annotation keywords,
# (3) presence of three consecutive Rep HMM hits (tripartite system).
REPABC_KW = frozenset({
    "repabc", "rep_abc", "repa_c", "repb_c", "repc_c",
    "trfa",                          # TrfA is RepA in repABC-type systems
    "plasmid replication protein a", # RepA helicase
    "plasmid replication protein b", # RepB centromere-binding
    "plasmid replication protein c", # RepC initiator
    "rep_cluster_1",                 # MOB-suite clusters for repABC family
    "rep_cluster_3",
    "rep_cluster_5",
})

# ── NHR E. coli detection ────────────────────────────────────────────────────
# Narrow-host-range E. coli plasmids (ColE1-type, IncF, IncI, etc.)
NHR_ECOLI_GENERA  = frozenset({"Escherichia", "Shigella"})
NHR_INC_GROUPS    = frozenset({"IncF", "IncFIA", "IncFIB", "IncFIC", "IncFII",
                                "IncI1", "IncI2", "IncB/O/K/Z"})
COLEI_ORIGIN_KW   = frozenset({"colei", "cole1", "rnai/rnaii", "rnai_rnaii",
                                "p15a", "pbr322", "pmb1", "puc", "col_e"})

# ── CDS classification keywords ──────────────────────────────────────────────
REP_KW  = frozenset({"replic", "initiator", "trfa", "repa", "repb", "repc",
                      "helix-turn-helix", "dna-binding", "rep protein",
                      "plasmid replication"})
PAR_KW  = frozenset({"partiti", "parb", "para", "spo0j", "sopa", "sopb",
                      "segregation", "plasmid stabiliz", "noc", "parg",
                      "ribbon-helix"})
MOB_KW  = frozenset({"mobili", "relax", "transfer", "conjugal", "trai", "traa",
                      "trbb", "trbc", "mob protein", "type iv"})
REC_KW  = frozenset({"transposase", "integrase", "recombinas", "resolvase",
                      "is element", "insertion sequence"})
AMR_KW  = frozenset({"resistance", "beta-lactam", "aminoglycoside", "tetracycline",
                      "sulfonamide", "chloramphenicol", "fluoroquinolone",
                      "macrolide", "rifampin", "colistin", "vancomycin"})

# Non-origin CDS: categories whose genes should be trimmed from window edges
NON_ORIGIN_TRIM_CATS = {"mobilization", "recombinase", "amr"}

# ── Motifs ───────────────────────────────────────────────────────────────────
WALKER_A_RE   = re.compile(r"G[ST].{3,5}GK[TS]")
DNAA_CONSENSUS = "TTATCCACA"
DNAA_RC        = str(Seq(DNAA_CONSENSUS).reverse_complement())

# ── Visualisation ────────────────────────────────────────────────────────────
TIER_COLORS = {
    "Tier 1: Agrobacterium-native": "#D32F2F",
    "Tier 2: Alphaproteobacteria":  "#F57C00",
    "Tier 3: Cross-phylum BHR":     "#388E3C",
    "Tier 4: Phylum-exclusive":     "#1976D2",
    "Tier 5: Diversity fill":       "#7B1FA2",
}
FEATURE_COLORS = {
    "rep":          "#E53935",
    "partition":    "#1565C0",
    "mobilization": "#2E7D32",
    "recombinase":  "#7B1FA2",
    "amr":          "#FF6F00",
    "structural":   "#FF8F00",
    "hypothetical": "#9E9E9E",
    "other":        "#78909C",
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1  PLSDB METADATA LOADING AND COLUMN NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────

# Map from PLSDB raw column names to internal standard names.
# Handles PLSDB v2.0+ TSV and the user's pre-processed CSV format.
_PLSDB_COL_MAP = {
    # PLSDB raw format
    "NUCCORE_ACC":             "Accession",
    "NUCCORE_Length":          "Length_bp",
    "NUCCORE_GC_Content":      "GC_content",
    "NUCCORE_Topology":        "Topology",
    "NUCCORE_Completeness":    "Completeness",
    "NUCCORE_Create_Date":     "Create_date",
    "NUCCORE_Definition":      "Description",
    "phylum_name":             "Phylum",
    "class_name":              "Class",
    "order_name":              "Order",
    "family_name":             "Family",
    "genus_name":              "Genus",
    "species_name":            "Species",
    "mobility":                "Predicted_mobility",
    "relaxase_type":           "Relaxase_types",
    "mpf_type":                "MPF_type",
    "rep_type":                "Replicon_types",
    "inc_type":                "Incompatibility_groups",
    "orit_type":               "OriT_types",
    # Some PLSDB versions
    "BIOSAMPLE_Latitude":      "Latitude",
    "BIOSAMPLE_Longitude":     "Longitude",
    "BIOSAMPLE_SampleType":    "Ecosystem",
}


def load_plsdb(path: str) -> pd.DataFrame:
    """
    Load a PLSDB metadata file (TSV or CSV) and normalise column names.
    Accepts raw PLSDB format OR the user's pre-processed CSV format.
    """
    sep = "\t" if path.endswith(".tsv") or path.endswith(".tab") else ","
    df  = pd.read_csv(path, sep=sep, low_memory=False)

    # Rename only columns that need renaming
    rename = {k: v for k, v in _PLSDB_COL_MAP.items() if k in df.columns}
    df = df.rename(columns=rename)

    # Ensure mandatory columns exist
    if "Accession" not in df.columns:
        raise ValueError(
            "Cannot find accession column. Expected 'NUCCORE_ACC' (raw PLSDB) "
            "or 'Accession' (pre-processed format)."
        )

    # Normalise GC content to fraction if stored as percentage
    if "GC_content" in df.columns and df["GC_content"].dropna().max() > 1.0:
        df["GC_content"] = df["GC_content"] / 100.0

    print(f"Loaded PLSDB: {len(df):,} entries, {df.columns.tolist()[:6]} ...")
    return df


def load_plsdb_2025(plsdb_dir: str) -> pd.DataFrame:
    """
    Load and join the PLSDB 2025 relational tables (nuccore + taxonomy + typing)
    and normalise to the internal column format.

    Expected files in plsdb_dir:
        nuccore.csv    – core plasmid records
        taxonomy.csv   – NCBI taxonomy per plasmid
        typing.csv     – MOB-suite replicon/relaxase typing
    """
    d = Path(plsdb_dir)
    nuccore  = pd.read_csv(str(d / "nuccore.csv"),  low_memory=False)
    taxonomy = pd.read_csv(str(d / "taxonomy.csv"), low_memory=False)
    typing   = pd.read_csv(str(d / "typing.csv"),   low_memory=False)

    # Join taxonomy on TAXONOMY_UID
    df = nuccore.merge(taxonomy, on="TAXONOMY_UID", how="left")
    # Join typing on NUCCORE_ACC
    df = df.merge(typing, on="NUCCORE_ACC", how="left")

    # Normalise to internal column names
    rename = {
        "NUCCORE_ACC":             "Accession",
        "NUCCORE_Length":          "Length_bp",
        "NUCCORE_GC":              "GC_content",
        "NUCCORE_Topology":        "Topology",
        "NUCCORE_Completeness":    "Completeness",
        "NUCCORE_CreateDate":      "Create_date",
        "NUCCORE_Description":     "Description",
        "TAXONOMY_phylum":         "Phylum",
        "TAXONOMY_class":          "Class",
        "TAXONOMY_order":          "Order",
        "TAXONOMY_family":         "Family",
        "TAXONOMY_genus":          "Genus",
        "TAXONOMY_species":        "Species",
        "predicted_mobility":      "Predicted_mobility",
        "relaxase_type(s)":        "Relaxase_types",
        "mpf_type":                "MPF_type",
        "rep_type(s)":             "Replicon_types",
        "predicted_host_range_overall_name": "Incompatibility_groups",
        "orit_type(s)":            "OriT_types",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    if "GC_content" in df.columns and df["GC_content"].dropna().max() > 1.0:
        df["GC_content"] = df["GC_content"] / 100.0

    print(f"Loaded PLSDB 2025: {len(df):,} entries (nuccore+taxonomy+typing joined)")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1b  OriVFinder CLASSES  (embedded from OriVFinder v1.02)
# Iteron, ZCurveAnalyzer, PatternFinder — pure-Python scoring of intergenic
# sequences adjacent to Rep (RIP) genes.
# Reference: doi.org/10.1093/nar/gkaf341
# ─────────────────────────────────────────────────────────────────────────────

class _Iteron:
    """Detect and score direct repeats (iterons) in a DNA sequence."""

    def __init__(self, min_iteron: int = 10, max_iteron: int = 22):
        self.min_iteron = min_iteron
        self.max_iteron = max_iteron

    @staticmethod
    def _dna_score(s1: str, s2: str) -> float:
        if len(s1) != len(s2): return 0.0
        return sum(2 if a == b else 0 for a, b in zip(s1, s2))

    @staticmethod
    def _entropy_score(seqs: list[str]) -> float:
        if not seqs: return 0.0
        positions = list(zip(*seqs))
        L = len(positions)
        total_e = 0.0
        for pos in positions:
            from collections import Counter as _C
            cnt = _C(pos)
            tot = sum(cnt.values())
            probs = [c / tot for c in cnt.values()]
            total_e += sum(1 - (-p * np.log2(p)) for p in probs if p > 0)
        return max(0, (total_e / L) * 100)

    @staticmethod
    def _cluster_index(distances: list[int]) -> float:
        if len(distances) <= 1: return 1.0
        span = sum(distances)
        return 1.0 / (1.0 + np.log(1.0 + span / len(distances)))

    def analyze(self, seq: str) -> tuple[float, list]:
        """
        Fast numpy-based iteron detector.

        Strategy: for each kmer length kl in [min_iteron, max_iteron], encode
        the DNA into a rolling hash using numpy views to count exact-match
        kmer occurrences in O(L) time.  Any kmer that appears ≥ 3 times with
        inter-copy distances ≤ 300 bp is reported as an iteron candidate.
        Score = kl × cluster_index × copies × (entropy/100), matching the
        original formula but computed ~100× faster.
        """
        seq = seq.upper()
        L   = len(seq)
        if L < self.min_iteron * 3:
            return 0.0, []

        # Encode A→0, C→1, G→2, T→3, other→4 using a lookup table
        _enc = np.frombuffer(
            bytes(seq.encode("ascii")), dtype=np.uint8
        )
        _lut = np.full(256, 4, dtype=np.uint8)
        for ch, val in zip(b"ACGT", range(4)):
            _lut[ch] = val
        enc = _lut[_enc]          # shape (L,), values 0-4

        MAX_GAP  = 300            # maximum inter-copy gap for a valid iteron
        MIN_COPIES = 3            # minimum copies (including first)
        all_results = []

        for kl in range(self.min_iteron, self.max_iteron + 1):
            if L < kl * MIN_COPIES:
                continue
            # Build all kmers as rows of shape (L-kl+1, kl)
            idx = np.arange(kl)[None, :] + np.arange(L - kl + 1)[:, None]
            kmers = enc[idx]                  # (N_pos, kl)
            # Skip any kmer containing an ambiguous base (value 4)
            valid = ~np.any(kmers == 4, axis=1)   # (N_pos,)
            if not valid.any():
                continue

            # Hash each kmer to a single integer (base-5)
            powers = (5 ** np.arange(kl, dtype=np.int64))
            hashes = (kmers * powers).sum(axis=1)   # (N_pos,)
            hashes[~valid] = -1                      # mark invalid

            # Find kmers that appear ≥ MIN_COPIES times
            from collections import Counter as _C
            h_valid  = hashes[valid]
            pos_valid = np.where(valid)[0]
            cnt = _C(h_valid.tolist())
            for h, n in cnt.items():
                if n < MIN_COPIES:
                    continue
                positions = sorted(pos_valid[hashes[pos_valid] == h].tolist())
                # Check that copies are within MAX_GAP of each other
                dists = [positions[i+1] - positions[i] - kl
                         for i in range(len(positions) - 1)]
                if any(d > MAX_GAP for d in dists):
                    continue
                ci  = self._cluster_index(dists)
                # Reconstruct the kmer string for entropy
                kmer_str = seq[positions[0]: positions[0] + kl]
                alns = [seq[p: p+kl] for p in positions]
                eb   = self._entropy_score(alns)
                final_score = kl * ci * n * (eb / 100)
                if final_score > 0:
                    all_results.append({
                        "kmer": kmer_str, "copies": n, "length": kl,
                        "clustering_index": ci, "entropy": eb,
                        "final_score": final_score,
                        "positions": [(p, p + kl) for p in positions],
                    })

        if not all_results:
            return 0.0, []
        best = max(all_results, key=lambda x: x["final_score"])
        return best["final_score"], all_results


class _ZCurveAnalyzer:
    """Z'-curve analysis to identify AT-rich replication initiation regions."""

    def __init__(self, sigma: float = 5.0):
        self.sigma = sigma

    @staticmethod
    def _regression_slope(x: np.ndarray, y: np.ndarray) -> float:
        xm, ym = x.mean(), y.mean()
        Sxx = ((x - xm)**2).sum()
        if Sxx == 0: return 0.0
        Sxy = ((x - xm) * (y - ym)).sum()
        return Sxy / Sxx

    def analyze(self, seq: str) -> tuple[float, dict | None]:
        if not _SCIPY_OK or len(seq) < 10:
            return 0.0, None
        seq = seq.upper()
        L = len(seq)
        Zn = np.zeros(L)
        sa = sc = sg = st = 0
        for i, b in enumerate(seq):
            if b == 'A': sa += 1
            elif b == 'C': sc += 1
            elif b == 'G': sg += 1
            elif b == 'T': st += 1
            Zn[i] = (sa + st) - (sg + sc)

        x = np.arange(1, L+1, dtype=float)
        k = self._regression_slope(x, Zn)
        Zn_prime = Zn - k * x
        Zn_smooth = gaussian_filter1d(Zn_prime, self.sigma)

        peaks   = argrelextrema(Zn_smooth, np.greater)[0]
        valleys = argrelextrema(-Zn_smooth, np.greater)[0]
        splits  = np.sort(np.concatenate(([0], peaks, valleys, [L-1])))

        max_at_score, best = -np.inf, None
        for i in range(len(splits)-1):
            s, e = int(splits[i]), int(splits[i+1])
            if e - s <= 0: continue
            sub = seq[s:e]
            at = (sub.count('A') + sub.count('T')) / len(sub)
            slope = self._regression_slope(x[s:e], Zn_smooth[s:e])
            sc = (e - s) * slope * at
            if sc > max_at_score:
                max_at_score = sc
                best = {"start": s, "end": e, "at_content": at, "slope": slope, "at_score": sc}

        return (max_at_score if best else 0.0), best


class _PatternFinder:
    """Find DnaA boxes, Fis/IHF binding sites, and other oriV motifs."""

    _IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]',
              'W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]',
              'H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
    _IUPAC_RC = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','S':'S','W':'W',
                 'K':'M','M':'K','B':'V','D':'H','H':'D','V':'B','N':'N'}

    _PATTERNS = {
        "DnaA_box":         {"seq": "TTATCCACA", "mm": 1,  "score": 20},
        "ctra":             {"seq": r"TTAA.{7}TTAA", "regex": True, "score": 5},
        "Fis_binding":      {"seq": "GNNNAWWWWWTNNNC", "iupac": True, "score": 5},
        "IHF_binding":      {"seq": "WATCAANNNNTTR",   "iupac": True, "score": 5},
    }

    @staticmethod
    def _mm(s1: str, s2: str) -> int:
        return sum(a != b for a, b in zip(s1, s2))

    def _iupac_re(self, p: str) -> str:
        return ''.join(self._IUPAC[b] for b in p)

    def _iupac_rc(self, p: str) -> str:
        return ''.join(self._IUPAC_RC[b] for b in reversed(p))

    def find_patterns(self, seq: str) -> tuple[float, list]:
        seq = seq.upper()
        rc_seq = str(Seq(seq).reverse_complement())
        total_score = 0.0
        hits = []
        seen = set()

        for name, info in self._PATTERNS.items():
            score = info["score"]
            if info.get("regex"):
                pat = info["seq"]
                for m in re.finditer(pat, seq):
                    k = (m.start(), m.end())
                    if k not in seen:
                        hits.append({"name": name, "start": m.start(), "end": m.end(),
                                     "strand": "+", "mismatches": 0})
                        total_score += score; seen.add(k)
                rc_pat = str(Seq(re.sub(r'\.\{.*?\}', lambda m: m.group(0).replace('.','.'), pat)).reverse_complement())
                for m in re.finditer(pat, rc_seq):
                    k = (len(seq)-m.end(), len(seq)-m.start())
                    if k not in seen:
                        hits.append({"name": name, "start": k[0], "end": k[1],
                                     "strand": "-", "mismatches": 0})
                        total_score += score; seen.add(k)
            elif info.get("iupac"):
                re_fwd = self._iupac_re(info["seq"])
                re_rev = self._iupac_re(self._iupac_rc(info["seq"]))
                for strand, target, re_p in [("+", seq, re_fwd), ("-", rc_seq, re_rev)]:
                    for m in re.finditer(re_p, target):
                        pos = (m.start(), m.end()) if strand == "+" else (len(seq)-m.end(), len(seq)-m.start())
                        if pos not in seen:
                            hits.append({"name": name, "start": pos[0], "end": pos[1],
                                         "strand": strand, "mismatches": 0})
                            total_score += score; seen.add(pos)
            else:
                pat_seq = info["seq"]; mm_max = info["mm"]
                plen = len(pat_seq)
                rc_pat = str(Seq(pat_seq).reverse_complement())
                for strand, target, p in [("+", seq, pat_seq), ("-", rc_seq, rc_pat)]:
                    for i in range(len(target) - plen + 1):
                        if self._mm(target[i:i+plen], p) <= mm_max:
                            pos = (i, i+plen) if strand == "+" else (len(seq)-i-plen, len(seq)-i)
                            if pos not in seen:
                                hits.append({"name": name, "start": pos[0], "end": pos[1],
                                             "strand": strand, "mismatches": self._mm(target[i:i+plen], p)})
                                mm = self._mm(target[i:i+plen], p)
                                total_score += score if mm == 0 else score // 2
                                seen.add(pos)
        return total_score, hits


def score_igs(seq: str,
              iteron_min: int = 10, iteron_max: int = 22,
              zcurve_sigma: float = 5.0) -> dict:
    """
    Score an intergenic sequence using OriVFinder's composite algorithm.

    Returns a dict with keys:
        iteron_score, zcurve_score, pattern_score, total_score,
        iteron_details, zcurve_details, pattern_hits
    """
    it_score, it_details = _Iteron(iteron_min, iteron_max).analyze(seq)
    zc_score, zc_details = _ZCurveAnalyzer(zcurve_sigma).analyze(seq)
    pt_score, pt_hits    = _PatternFinder().find_patterns(seq)
    total = it_score + zc_score + abs(pt_score)
    return {
        "iteron_score":  it_score,
        "zcurve_score":  zc_score,
        "pattern_score": pt_score,
        "total_score":   total,
        "iteron_details": it_details,
        "zcurve_details": zc_details,
        "pattern_hits":   pt_hits,
    }


def extract_igs_near_rips(plasmid_seq: str, orfs: list[dict],
                          rep_orf_ids: set[str],
                          min_igs: int = 50) -> list[dict]:
    """
    Extract intergenic sequences adjacent to Rep (RIP) genes.

    Strategy (mirrors OriVFinder AnnotationAnalyzer):
      1. Sort all ORFs by start position.
      2. Compute gaps between consecutive ORFs — these are IGSs.
      3. Also include the circular wrap-around gap (last ORF end → first ORF start).
      4. Flag each IGS as "near_rip" if it is between two genes where at
         least one of the flanking genes is a Rep ORF, or is adjacent to
         the Rep gene within ±2 gene positions.

    Returns list of dicts:
        igs_seq, igs_start, igs_end, near_rip, flanking_rep_ids
    """
    L = len(plasmid_seq)
    if not orfs:
        return [{"igs_seq": plasmid_seq, "igs_start": 0, "igs_end": L,
                 "near_rip": False, "flanking_rep_ids": []}]

    sorted_orfs = sorted(orfs, key=lambda o: o["start"])
    n = len(sorted_orfs)
    igs_list = []

    # IGSs between consecutive ORFs
    for i in range(n - 1):
        gap_s = sorted_orfs[i]["end"]
        gap_e = sorted_orfs[i+1]["start"]
        if gap_e - gap_s < min_igs:
            continue
        near = (sorted_orfs[i]["id"] in rep_orf_ids or
                sorted_orfs[i+1]["id"] in rep_orf_ids)
        # also ±2 gene neighbours
        reps_nearby = [sorted_orfs[j]["id"]
                       for j in range(max(0, i-1), min(n, i+3))
                       if sorted_orfs[j]["id"] in rep_orf_ids]
        near = near or bool(reps_nearby)
        igs_list.append({
            "igs_seq":          plasmid_seq[gap_s:gap_e],
            "igs_start":        gap_s,
            "igs_end":          gap_e,
            "near_rip":         near,
            "flanking_rep_ids": reps_nearby,
        })

    # Circular wrap-around gap
    wrap_s = sorted_orfs[-1]["end"]
    wrap_e = sorted_orfs[0]["start"]
    if wrap_s != wrap_e:
        circ_seq = (plasmid_seq[wrap_s:] + plasmid_seq[:wrap_e]) if wrap_e < wrap_s \
                   else plasmid_seq[wrap_s:wrap_e]
        if len(circ_seq) >= min_igs:
            near = (sorted_orfs[-1]["id"] in rep_orf_ids or
                    sorted_orfs[0]["id"]  in rep_orf_ids)
            igs_list.append({
                "igs_seq":          circ_seq,
                "igs_start":        wrap_s,
                "igs_end":          wrap_e,
                "near_rip":         near,
                "flanking_rep_ids": [o["id"] for o in [sorted_orfs[-1], sorted_orfs[0]]
                                     if o["id"] in rep_orf_ids],
                "circular":         True,
            })

    return igs_list


def find_best_oriv(acc: str, plasmid_seq: str,
                   orfs: list[dict], hmm_hits: dict,
                   rep_buffer: int = REP_BUFFER_BP,
                   max_window: int = ORIGIN_MAX_BP) -> dict | None:
    """
    Identify the best origin window for a plasmid using OriVFinder scoring.

    1. Find Rep ORF from HMM hits.
    2. Extract all IGSs.
    3. Score each IGS with OriVFinder composite score.
    4. Select the best-scoring IGS near the Rep gene.
    5. Compute origin window = best_IGS + Rep_gene + rep_buffer on each side.

    Returns a dict with origin coordinates, scores, and the best IGS sequence.
    Returns None if no Rep ORF is found.
    """
    hit_info = hmm_hits.get(acc, {})
    if hit_info.get("hmm_tier", "none") == "none":
        return None

    # Find the best Rep ORF
    best_hit = None
    if hit_info.get("hits"):
        best_hit = min(hit_info["hits"], key=lambda h: h["evalue"])
    if not best_hit:
        return None

    rep_orf = next((o for o in orfs if o["id"] == best_hit["orf_id"]), None)
    if not rep_orf:
        return None

    rep_start = rep_orf["start"]
    rep_end   = rep_orf["end"]
    rep_ids   = {o["id"] for o in orfs
                 if abs(o["start"] - rep_start) < 5 or o["id"] == rep_orf["id"]}

    # Extract and score all IGSs
    igs_list = extract_igs_near_rips(plasmid_seq, orfs, rep_ids)

    if not igs_list:
        # Fallback: define window from Rep gene position
        ws, we = define_origin_window(rep_start, rep_end, len(plasmid_seq),
                                       max_window, rep_buffer)
        return {
            "Accession": acc, "Rep_orf_start": rep_start, "Rep_orf_end": rep_end,
            "Rep_hmm_name": best_hit["name"], "Rep_hmm_evalue": best_hit["evalue"],
            "Origin_start_bp": ws, "Origin_end_bp": we,
            "Origin_span_bp": we - ws,
            "OriVFinder_evidence": "Rep-only",
            "OriVFinder_total_score": 0.0,
            "OriVFinder_iteron_score": 0.0,
            "OriVFinder_zcurve_score": 0.0,
            "OriVFinder_pattern_score": 0.0,
        }

    # Score only near-Rep IGSs first (avoids scoring all ~80+ IGSs per plasmid)
    near_igs  = [g for g in igs_list if g.get("near_rip")]
    score_set = near_igs if near_igs else igs_list  # fallback: score all
    scored = []
    for igs in score_set:
        s = score_igs(igs["igs_seq"])
        scored.append({**igs, **s})

    best_igs = max(scored, key=lambda x: x["total_score"])

    # Build origin window: encompass Rep gene + best IGS + buffer
    plen  = len(plasmid_seq)
    igs_s = best_igs["igs_start"]
    igs_e = best_igs["igs_end"]

    # Window spans: min(rep_start, igs_s) - buffer  to  max(rep_end, igs_e) + buffer
    raw_s = min(rep_start, igs_s)
    raw_e = max(rep_end,   igs_e)
    ws    = max(0,    raw_s - rep_buffer)
    we    = min(plen, raw_e + rep_buffer)

    # Clamp to max_window (centred on the Rep gene)
    if (we - ws) > max_window:
        ws, we = define_origin_window(rep_start, rep_end, plen, max_window, rep_buffer)

    # Re-check buffer around Rep gene
    if rep_start - ws < rep_buffer:
        ws = max(0, rep_start - rep_buffer)
    if we - rep_end < rep_buffer:
        we = min(plen, rep_end + rep_buffer)

    # Ensure 200 bp buffer between window edge and ANY origin-critical ORF
    # (i.e. any ORF whose body falls inside the window).  If an ORF starts
    # within 200 bp of the left edge, extend left; likewise for right edge.
    for orf in orfs:
        os_, oe_ = orf["start"], orf["end"]
        # Only care about ORFs that overlap or are inside the window
        if oe_ <= ws or os_ >= we:
            continue
        if os_ - ws < rep_buffer and os_ > ws:
            ws = max(0, os_ - rep_buffer)
        if we - oe_ < rep_buffer and oe_ < we:
            we = min(plen, oe_ + rep_buffer)

    # Clamp to max_window after buffer expansion
    if (we - ws) > max_window:
        ws, we = define_origin_window(rep_start, rep_end, plen, max_window, rep_buffer)

    evidence = "near-RIP" if near_igs else "best-IGS"
    return {
        "Accession":               acc,
        "Rep_orf_start":           rep_start,
        "Rep_orf_end":             rep_end,
        "Rep_hmm_name":            best_hit["name"],
        "Rep_hmm_evalue":          best_hit["evalue"],
        "Origin_start_bp":         ws,
        "Origin_end_bp":           we,
        "Origin_span_bp":          we - ws,
        "OriVFinder_evidence":     evidence,
        "OriVFinder_total_score":  best_igs["total_score"],
        "OriVFinder_iteron_score": best_igs["iteron_score"],
        "OriVFinder_zcurve_score": best_igs["zcurve_score"],
        "OriVFinder_pattern_score":best_igs["pattern_score"],
        "OriVFinder_igs_start":    igs_s,
        "OriVFinder_igs_end":      igs_e,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2  METADATA-LEVEL FILTERING
# ─────────────────────────────────────────────────────────────────────────────

def _is_repabc(row: pd.Series) -> bool:
    """Return True if this plasmid appears to use a repABC replication system."""
    targets = " ".join(str(row.get(c, "") or "")
                       for c in ("Replicon_types", "Origin_RIP_type",
                                 "Origin_elements", "Description")).lower()
    return any(k in targets for k in REPABC_KW)


def _is_nhr_ecoli(row: pd.Series) -> bool:
    """
    Return True if this is a narrow-host-range E. coli plasmid.

    Criteria:
      (a) Host is Escherichia or Shigella AND incompatibility group is
          an IncF/IncI-family group (NHR-associated), OR
      (b) Replicon/RIP type contains ColE1/RNAI-origin keywords (NHR regardless
          of current host — these replicate only in enteric bacteria).
    """
    genus = str(row.get("Genus", "") or "").strip()
    inc   = str(row.get("Incompatibility_groups", "") or "").lower()
    rtype = str(row.get("Replicon_types", "") or "").lower()
    rip   = str(row.get("Origin_RIP_type", "") or "").lower()

    ecoli_host = genus in NHR_ECOLI_GENERA
    nhr_inc    = any(g.lower() in inc for g in NHR_INC_GROUPS)
    colei_ori  = any(k in rtype + " " + rip for k in COLEI_ORIGIN_KW)

    return (ecoli_host and nhr_inc) or colei_ori


def metadata_filter(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply metadata-level filters before OriV finding.

    Keeps entries that are:
      • Circular topology (circular plasmids only)
      • Bacterial kingdom
      • NOT repABC family
      • NOT narrow-host-range E. coli
      • Length ≤ 500 kb  (very large plasmids are genomic islands, not useful)
      • Length ≥ 1 kb    (tiny sequences unreliable)
    """
    n0 = len(df)
    report = {}

    # Circular topology
    if "Topology" in df.columns:
        df = df[df["Topology"].str.lower().str.contains("circ", na=False)]
    report["circular"] = len(df)

    # Size bounds (plasmid-level)
    if "Length_bp" in df.columns:
        df = df[(df["Length_bp"] >= 1_000) & (df["Length_bp"] <= 500_000)]
    report["size_bounds"] = len(df)

    # Not repABC
    before = len(df)
    df = df[~df.apply(_is_repabc, axis=1)]
    report["not_repABC"]  = len(df)
    report["_repABC_removed"] = before - len(df)

    # Not NHR E. coli
    before = len(df)
    df = df[~df.apply(_is_nhr_ecoli, axis=1)]
    report["not_NHR_Ecoli"]  = len(df)
    report["_NHR_removed"]   = before - len(df)

    print(f"\nMetadata filtering ({n0:,} → {len(df):,}):")
    print(f"  After circular filter:      {report['circular']:,}")
    print(f"  After size filter:          {report['size_bounds']:,}")
    print(f"  After repABC removal:       {report['not_repABC']:,}  "
          f"(removed {report['_repABC_removed']:,})")
    print(f"  After NHR E. coli removal:  {report['not_NHR_Ecoli']:,}  "
          f"(removed {report['_NHR_removed']:,})")

    return df.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3  OriV FINDING
# ─────────────────────────────────────────────────────────────────────────────
# Strategy:
#   1. For each plasmid sequence, predict ORFs with pyrodigal (meta mode).
#   2. Write all ORF protein translations to a temporary FASTA.
#   3. Run hmmsearch (HMMER) once against RIP.hmm for the entire batch.
#   4. Parse the table output to identify which plasmids have a Rep protein hit.
#   5. Record the genomic coordinates of each Rep-protein ORF (needed for
#      origin window definition in Section 4).

def _run_hmmsearch(proteins_fasta: str, hmm_path: str,
                   evalue: float = 1e-5) -> dict[str, list[dict]]:
    """
    Run hmmsearch and return hits as {accession: [{orf_id, evalue, score, name}, ...]}.
    Requires hmmsearch on PATH.
    """
    if not shutil.which("hmmsearch"):
        raise RuntimeError("hmmsearch not found on PATH. Install HMMER 3.x.")

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
            orf_id   = parts[0]          # e.g. ACC_orf3
            hmm_name = parts[2]
            evalue_v = float(parts[4])
            score_v  = float(parts[5])
            # Recover accession: last token after splitting off the orf suffix
            acc = orf_id.rsplit("_orf", 1)[0]
            hits.setdefault(acc, []).append({
                "orf_id":   orf_id,
                "name":     hmm_name,
                "evalue":   evalue_v,
                "score":    score_v,
            })
    return hits


def find_oriv_in_plsdb(meta_df: pd.DataFrame,
                       plsdb_fasta: str | None,
                       hmm_path: str,
                       ncbi_email: str,
                       ncbi_api_key: str | None = None,
                       batch_size: int = 500,
                       hmm_evalue: float = 1e-5) -> pd.DataFrame:
    """
    Identify origins of replication in the filtered PLSDB candidates.

    If plsdb_fasta is provided: process sequences locally (fast).
    Otherwise: batch-fetch from NCBI (slow, rate-limited).

    Returns an expanded DataFrame with per-accession Rep protein location columns:
        Rep_orf_start, Rep_orf_end, Rep_orf_strand, Rep_hmm_name, Rep_hmm_evalue
    """
    Entrez.email = ncbi_email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    accessions = meta_df["Accession"].tolist()
    n = len(accessions)
    print(f"\nOriV finding across {n:,} candidates ...")

    # Collect ORF info by accession: {acc: [(start,end,strand,aa), ...]}
    all_orfs: dict[str, list[dict]] = {}
    # Collect plasmid sequences for window extraction later
    plasmid_seqs: dict[str, str] = {}

    # ── Local FASTA path ─────────────────────────────────────────────────────
    if plsdb_fasta:
        print(f"  Streaming {plsdb_fasta} ...")
        accs_needed = set(accessions)
        finder = pyrodigal.GeneFinder(meta=True)

        processed = 0
        with open(plsdb_fasta) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                # Match by accession — PLSDB FASTA IDs may have version suffixes
                acc = record.id.split()[0]
                # Strip version if needed (e.g. NZ_CP000001.1 → NZ_CP000001)
                # but keep the versioned form to match metadata
                if acc not in accs_needed:
                    # Try without version suffix
                    acc_no_ver = acc.rsplit(".", 1)[0]
                    if acc_no_ver in accs_needed:
                        acc = acc_no_ver
                    else:
                        continue

                seq_str = str(record.seq).upper()
                plasmid_seqs[acc] = seq_str

                genes = finder.find_genes(seq_str.encode())
                orfs  = []
                for j, g in enumerate(genes):
                    aa = g.translate().rstrip("*")
                    if len(aa) < 30:
                        continue
                    orfs.append({
                        "id":     f"{acc}_orf{j+1}",
                        "start":  g.begin,
                        "end":    g.end,
                        "strand": "+" if g.strand == 1 else "-",
                        "aa":     aa,
                    })
                all_orfs[acc] = orfs
                processed += 1
                if processed % 5000 == 0:
                    print(f"    ... {processed:,}/{n:,} sequences processed")

        print(f"  Sequences processed: {processed:,}")

    # ── NCBI batch-fetch fallback ─────────────────────────────────────────────
    else:
        print("  No PLSDB FASTA provided — batch-fetching from NCBI. "
              "This may take several hours for large datasets.")
        finder = pyrodigal.GeneFinder(meta=True)
        sleep  = NCBI_SLEEP_KEY if ncbi_api_key else NCBI_SLEEP

        for i in range(0, n, batch_size):
            batch_accs = accessions[i: i + batch_size]
            ids_str    = ",".join(batch_accs)
            try:
                handle = Entrez.efetch(db="nucleotide", id=ids_str,
                                       rettype="fasta", retmode="text")
                for rec in SeqIO.parse(StringIO(handle.read()), "fasta"):
                    acc     = rec.id.split()[0].rsplit(".", 1)[0]
                    seq_str = str(rec.seq).upper()
                    plasmid_seqs[acc] = seq_str
                    genes = finder.find_genes(seq_str.encode())
                    orfs  = []
                    for j, g in enumerate(genes):
                        aa = g.translate().rstrip("*")
                        if len(aa) < 30:
                            continue
                        orfs.append({
                            "id":     f"{acc}_orf{j+1}",
                            "start":  g.begin,
                            "end":    g.end,
                            "strand": "+" if g.strand == 1 else "-",
                            "aa":     aa,
                        })
                        all_orfs[acc] = orfs
                handle.close()
            except Exception as e:
                print(f"    Batch {i//batch_size+1} error: {e}")
            time.sleep(sleep)
            if (i // batch_size + 1) % 10 == 0:
                print(f"    ... {min(i+batch_size, n):,}/{n:,}")

    # ── Write all ORF proteins to temp FASTA and run hmmsearch ───────────────
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                     delete=False, prefix="oriv_orfs_") as tmp:
        prot_fasta = tmp.name
        total_prots = 0
        for acc, orfs in all_orfs.items():
            for orf in orfs:
                tmp.write(f">{orf['id']}\n{orf['aa']}\n")
                total_prots += 1

    print(f"  Total ORFs written: {total_prots:,}")
    print(f"  Running hmmsearch (E ≤ {hmm_evalue}) ...")
    hmm_hits = _run_hmmsearch(prot_fasta, hmm_path, hmm_evalue)
    os.unlink(prot_fasta)
    if os.path.exists(prot_fasta + ".tbl"):
        os.unlink(prot_fasta + ".tbl")

    # ── Filter for repABC by HMM hit name (tripartite detection) ─────────────
    # If all top-3 hits for a plasmid match repABC-family names, discard.
    def _hits_are_repabc(hits: list[dict]) -> bool:
        names_lc = " ".join(h["name"].lower() for h in hits[:5])
        return any(k in names_lc for k in REPABC_KW)

    repabc_by_hmm: set[str] = set()
    for acc, hits in hmm_hits.items():
        if _hits_are_repabc(hits):
            repabc_by_hmm.add(acc)
    if repabc_by_hmm:
        print(f"  Additional repABC plasmids flagged by HMM: {len(repabc_by_hmm):,}")

    # ── Build per-accession Rep protein location table ────────────────────────
    # For each accession with HMM hits, find the best-scoring ORF and its coords.
    orf_lookup: dict[str, dict] = {}   # orf_id → orf dict
    for acc, orfs in all_orfs.items():
        for orf in orfs:
            orf_lookup[orf["id"]] = orf

    rep_locs: list[dict] = []
    for acc in accessions:
        if acc in repabc_by_hmm:
            continue
        hits = hmm_hits.get(acc, [])
        if not hits:
            continue
        best_hit = min(hits, key=lambda h: h["evalue"])
        orf      = orf_lookup.get(best_hit["orf_id"])
        if orf is None:
            continue
        meta_row = meta_df[meta_df["Accession"] == acc]
        if meta_row.empty:
            continue
        rep_locs.append({
            **meta_row.iloc[0].to_dict(),
            "Rep_orf_start":  orf["start"],
            "Rep_orf_end":    orf["end"],
            "Rep_orf_strand": orf["strand"],
            "Rep_hmm_name":   best_hit["name"],
            "Rep_hmm_evalue": best_hit["evalue"],
            "Plasmid_seq":    plasmid_seqs.get(acc, ""),
        })

    result_df = pd.DataFrame(rep_locs)
    print(f"  Origins with Rep protein detected: {len(result_df):,}")
    return result_df


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4  ORIGIN WINDOW DEFINITION WITH 200 bp REP-BUFFER
# ─────────────────────────────────────────────────────────────────────────────

def define_origin_window(rep_start: int, rep_end: int,
                         plasmid_len: int,
                         max_window: int = ORIGIN_MAX_BP,
                         buffer: int = REP_BUFFER_BP) -> tuple[int, int]:
    """
    Define a ≤ max_window bp origin window that:
      (a) contains the Rep gene entirely,
      (b) has ≥ buffer bp of sequence upstream AND downstream of the Rep gene.

    Centering strategy:
      • Compute available space: max_window - rep_gene_length.
      • Split that equally between upstream and downstream flanks,
        each clamped to ≥ buffer bp.
      • If the Rep gene alone + 2×buffer > max_window, the window exceeds
        max_window — accept this (the origin is still ≤ 6 kb in practice
        because Rep genes are < 3 kb).

    Returns (window_start, window_end) in 0-based half-open coordinates.
    """
    rep_len  = rep_end - rep_start
    avail    = max(max_window - rep_len, 2 * buffer)
    flank    = avail // 2
    # Ensure flank ≥ buffer
    flank_up = max(flank, buffer)
    flank_dn = max(avail - flank_up, buffer)

    w_start = max(0,           rep_start - flank_up)
    w_end   = min(plasmid_len, rep_end   + flank_dn)

    # Second pass: if one side was truncated by plasmid boundary, extend the other
    if w_start == 0 and (w_end - w_start) < max_window:
        w_end = min(plasmid_len, max_window)
    if w_end == plasmid_len and (w_end - w_start) < max_window:
        w_start = max(0, plasmid_len - max_window)

    # Final buffer sanity checks
    if rep_start - w_start < buffer:
        w_start = max(0, rep_start - buffer)
    if w_end - rep_end < buffer:
        w_end = min(plasmid_len, rep_end + buffer)

    return w_start, w_end


def add_origin_windows(oriv_df: pd.DataFrame) -> pd.DataFrame:
    """Apply define_origin_window to all rows and keep only ≤ ORIGIN_MAX_BP windows."""
    rows = []
    for _, row in oriv_df.iterrows():
        plen = int(row["Length_bp"])
        rs   = int(row["Rep_orf_start"])
        re_  = int(row["Rep_orf_end"])
        ws, we = define_origin_window(rs, re_, plen)
        span = we - ws
        if span > ORIGIN_MAX_BP:
            continue    # skip if Rep gene alone pushes beyond cap (extremely rare)
        seq  = row.get("Plasmid_seq", "")
        orig_seq = seq[ws:we] if seq else ""
        row = row.to_dict()
        row.update({
            "Origin_start_bp":        ws,
            "Origin_end_bp":          we,
            "Origin_span_bp":         span,
            "Origin_GenBank_location": f"{ws}..{we}",
            "Origin_sequence":         orig_seq,
        })
        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"\nOrigins ≤ {ORIGIN_MAX_BP} bp with 200 bp Rep buffer: {len(df):,}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5  NON-ORIGIN CDS TRIMMING
# ─────────────────────────────────────────────────────────────────────────────

def _cds_category(product: str) -> str:
    pl = product.lower()
    if any(k in pl for k in REP_KW):                     return "rep"
    if any(k in pl for k in PAR_KW):                     return "partition"
    if any(k in pl for k in MOB_KW):                     return "mobilization"
    if any(k in pl for k in REC_KW):                     return "recombinase"
    if any(k in pl for k in AMR_KW):                     return "amr"
    if "hypothetical" in pl or "unknown" in pl:          return "hypothetical"
    return "other"


def _is_non_origin_cds(product: str) -> bool:
    cat = _cds_category(product)
    return cat in NON_ORIGIN_TRIM_CATS


def trim_origin_window(w_start: int, w_end: int,
                       gb_features: list[dict],
                       rep_start: int, rep_end: int,
                       buffer: int = REP_BUFFER_BP) -> tuple[int, int]:
    """
    Trim non-origin CDS from the edges of an origin window.

    Rules:
      • Only trim CDS in categories: mobilization, recombinase, amr.
      • Trimming proceeds from each edge inward until the first
        non-trimmable CDS (rep, partition, other).
      • NEVER trim inside the Rep-gene protection zone:
          [rep_start − buffer, rep_end + buffer]
      • Trimming at the left edge sets new_start = trimmed_CDS.end
      • Trimming at the right edge sets new_end  = trimmed_CDS.start

    Returns (trimmed_start, trimmed_end).
    """
    safe_l = max(w_start, rep_start - buffer)
    safe_r = min(w_end,   rep_end   + buffer)

    feats = sorted(gb_features, key=lambda f: f["start"])
    new_l, new_r = w_start, w_end

    # Trim left edge
    for feat in feats:
        # Feature must be entirely left of the Rep protection zone
        if feat["end"] > safe_l:
            break
        if _is_non_origin_cds(feat.get("product", "")):
            candidate = min(feat["end"], safe_l)
            if candidate > new_l:
                new_l = candidate
        else:
            break   # non-trimmable CDS found — stop left trimming

    # Trim right edge
    for feat in reversed(feats):
        # Feature must be entirely right of the Rep protection zone
        if feat["start"] < safe_r:
            break
        if _is_non_origin_cds(feat.get("product", "")):
            candidate = max(feat["start"], safe_r)
            if candidate < new_r:
                new_r = candidate
        else:
            break   # non-trimmable CDS found — stop right trimming

    return new_l, new_r


def apply_trimming_from_annotations(oriv_df: pd.DataFrame,
                                    ncbi_email: str,
                                    ncbi_api_key: str | None = None) -> pd.DataFrame:
    """
    Fetch GenBank annotations for each origin, then apply CDS trimming.
    Only fetches the origin window + small flanking to minimise data transfer.

    Updates Origin_start_bp, Origin_end_bp, Origin_span_bp, Origin_sequence.
    Also adds an Origin_elements column describing what CDS are in the window.
    """
    Entrez.email = ncbi_email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    sleep  = NCBI_SLEEP_KEY if ncbi_api_key else NCBI_SLEEP
    n      = len(oriv_df)
    print(f"\nFetching annotations + trimming windows for {n:,} candidates ...")

    updated_rows = []
    for i, (_, row) in enumerate(oriv_df.iterrows()):
        acc   = row["Accession"]
        ws    = int(row["Origin_start_bp"])
        we    = int(row["Origin_end_bp"])
        rs    = int(row["Rep_orf_start"])
        re_   = int(row["Rep_orf_end"])
        plen  = int(row["Length_bp"])
        seq   = row.get("Origin_sequence", "")

        gb_features: list[dict] = []

        # Fetch annotation for the origin window region
        try:
            fetch_s = max(1, ws - 200 + 1)      # 1-based for NCBI, with small extra flank
            fetch_e = min(plen, we + 200)
            handle  = Entrez.efetch(db="nucleotide", id=acc, rettype="gb",
                                    retmode="text", seq_start=fetch_s, seq_stop=fetch_e)
            record  = SeqIO.read(StringIO(handle.read()), "genbank")
            handle.close()

            # Offset: feature coords in the fetched sub-record are relative to fetch_s
            offset = fetch_s - 1  # convert back to 0-based absolute coords
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                product    = feat.qualifiers.get("product", ["unknown"])[0]
                protein_id = feat.qualifiers.get("protein_id", [""])[0]
                feat_start = int(feat.location.start) + offset
                feat_end   = int(feat.location.end)   + offset
                gb_features.append({
                    "product":    product,
                    "protein_id": protein_id,
                    "start":      feat_start,
                    "end":        feat_end,
                    "strand":     feat.location.strand,
                    "category":   _cds_category(product),
                })
        except Exception:
            pass   # If fetch fails, use untrimmed window

        # Apply CDS trimming
        new_s, new_e = trim_origin_window(ws, we, gb_features, rs, re_)

        # Extract trimmed sequence
        plasmid_seq = row.get("Plasmid_seq", "")
        if plasmid_seq:
            origin_seq = plasmid_seq[new_s:new_e]
        elif seq:
            # Adjust relative to existing window
            rel_s = new_s - ws
            rel_e = new_e - ws
            origin_seq = seq[max(0, rel_s): rel_e]
        else:
            origin_seq = ""

        # Build origin elements string
        rep_cds = [f["product"] for f in gb_features if f["category"] == "rep"]
        par_cds = [f["product"] for f in gb_features if f["category"] == "partition"]
        all_cds_in_window = [f["product"] for f in gb_features
                             if new_s <= f["start"] < new_e]
        elements_str = "; ".join(all_cds_in_window[:6]) if all_cds_in_window else ""

        row_d = row.to_dict()
        row_d.update({
            "Origin_start_bp":        new_s,
            "Origin_end_bp":          new_e,
            "Origin_span_bp":         new_e - new_s,
            "Origin_GenBank_location": f"{new_s}..{new_e}",
            "Origin_sequence":         origin_seq,
            "Origin_elements":         elements_str,
            "n_rep_elements":          len(rep_cds),
            "n_par_elements":          len(par_cds),
            "has_identified_rip":      len(rep_cds) > 0,
            "has_par":                 len(par_cds) > 0,
            "_gb_features":            gb_features,   # kept internally, dropped before CSV
        })
        updated_rows.append(row_d)

        if (i + 1) % 100 == 0 or i == 0:
            print(f"  [{i+1}/{n}] {acc}: window {ws}→{new_s} .. {we}→{new_e}"
                  f"  ({new_e-new_s} bp, {len(gb_features)} CDS)")
        time.sleep(sleep)

    df = pd.DataFrame(updated_rows)
    # Final size filter after trimming
    df = df[df["Origin_span_bp"].between(500, ORIGIN_MAX_BP)].reset_index(drop=True)
    print(f"After trimming: {len(df):,} origins ≤ {ORIGIN_MAX_BP} bp")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6  DEDUPLICATION AND FILTERING
# ─────────────────────────────────────────────────────────────────────────────

def deduplicate_and_build_metadata(oriv_df: pd.DataFrame,
                                    max_per_type: int = 1) -> pd.DataFrame:
    """
    Deduplicate by RIP type: keep the top N representatives per RIP type,
    ranked by diversity (GC distance from Agrobacterium, origin span).
    Also assigns Origin_RIP_type from the Rep HMM hit name.

    max_per_type: how many representatives to keep per RIP type.
        1 = strict dedup (original behaviour)
        >1 = keep up to N diverse representatives per type
    """
    # Derive RIP type from HMM hit name
    if "Origin_RIP_type" not in oriv_df.columns:
        oriv_df["Origin_RIP_type"] = oriv_df.get("Rep_hmm_name", pd.Series("", index=oriv_df.index))

    # Keep only confirmed Rep-containing origins
    has_rip = oriv_df["Origin_RIP_type"].notna() & (oriv_df["Origin_RIP_type"] != "")
    oriv_df = oriv_df[has_rip].copy()

    # Rank within each RIP type: minimise GC distance from Agrobacterium
    def _top_reps(group):
        g = group.copy()
        g["_gc_dist"] = abs(g["GC_content"] - AGRO_GC)
        mx = g["Origin_span_bp"].max() or 1
        g["_sel"] = g["_gc_dist"] - 0.1 * (g["Origin_span_bp"] / mx)
        return g.nsmallest(max_per_type, "_sel")

    best = oriv_df.groupby("Origin_RIP_type", group_keys=False).apply(_top_reps).reset_index(drop=True)
    print(f"\nAfter deduplication: {len(best):,} origins ({best['Origin_RIP_type'].nunique()} RIP types, "
          f"up to {max_per_type} per type)")
    return best


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7  TIERED DIVERSITY SELECTION
# ─────────────────────────────────────────────────────────────────────────────

def compute_rip_stats(df: pd.DataFrame) -> pd.DataFrame:
    stats = df.groupby("Origin_RIP_type").agg(
        count     = ("Accession", "count"),
        n_genera  = ("Genus",     "nunique"),
        n_phyla   = ("Phylum",    "nunique"),
        n_species = ("Species",   "nunique"),
        has_alpha = ("Class",     lambda x: (x == "Alphaproteobacteria").any()),
        has_agro  = ("Genus",     lambda x: (x == "Agrobacterium").any()),
        mob_types = ("Predicted_mobility", "nunique"),
    ).reset_index()
    stats["diversity_score"] = (
        np.log1p(stats["n_genera"]) * 2.0 +
        np.log1p(stats["n_phyla"])  * 3.0 +
        np.log1p(stats["count"])    * 0.5 +
        stats["has_alpha"].astype(int) * 3.0 +
        stats["has_agro"].astype(int)  * 5.0 +
        stats["mob_types"]             * 0.5
    )
    return stats.sort_values("diversity_score", ascending=False).reset_index(drop=True)


def select_diversity_library(df: pd.DataFrame,
                              rip_stats: pd.DataFrame,
                              target_bp:    int = 500_000,
                              min_origins:  int = 100) -> pd.DataFrame:
    """
    Tiered greedy selection to maximise diversity.

    Tier 1  – RIP types found in Agrobacterium
    Tier 2  – Other Alphaproteobacteria RIP types
    Tier 3  – Cross-phylum broad-host-range (≥ 3 phyla)
    Tier 4  – One representative per underrepresented phylum
    Tier 5  – Remainder by diversity score until targets met
    """
    phylum_exclusive: dict[str, set] = {}
    for phylum in df["Phylum"].dropna().unique():
        ph  = set(df[df["Phylum"] == phylum]["Origin_RIP_type"])
        oth = set(df[df["Phylum"] != phylum]["Origin_RIP_type"])
        excl = ph - oth
        if excl:
            phylum_exclusive[phylum] = excl

    tier1 = set(rip_stats[rip_stats["has_agro"]]["Origin_RIP_type"])
    tier2 = set(rip_stats[rip_stats["has_alpha"]]["Origin_RIP_type"]) - tier1
    tier3 = set(rip_stats[rip_stats["n_phyla"] >= 3]["Origin_RIP_type"]) - tier1 - tier2

    tier4: set = set()
    for phylum, excl in phylum_exclusive.items():
        if phylum == "Pseudomonadota":
            continue
        cands = rip_stats[rip_stats["Origin_RIP_type"].isin(excl)]
        if len(cands):
            tier4.add(cands.iloc[0]["Origin_RIP_type"])
    tier4 -= tier1 | tier2 | tier3

    rip_order = rip_stats["Origin_RIP_type"].tolist()
    rows, labels, seen_acc, total_bp = [], [], set(), 0

    def _add(rips, label, budget=False):
        nonlocal total_bp
        for rip in rips:
            if budget and total_bp >= target_bp and len(rows) >= min_origins:
                return
            sub = df[df["Origin_RIP_type"] == rip]
            if sub.empty:
                continue
            for _, rep in sub.iterrows():
                acc = rep["Accession"]
                if acc in seen_acc:
                    continue
                if budget and total_bp >= target_bp and len(rows) >= min_origins:
                    return
                rows.append(rep); labels.append(label)
                seen_acc.add(acc); total_bp += int(rep["Origin_span_bp"])

    _add(sorted(tier1), "Tier 1: Agrobacterium-native")
    _add([r for r in rip_order if r in tier2], "Tier 2: Alphaproteobacteria")
    _add([r for r in rip_order if r in tier3], "Tier 3: Cross-phylum BHR")
    _add(sorted(tier4), "Tier 4: Phylum-exclusive")
    seen_rips = {r["Origin_RIP_type"] for r in rows}
    _add([r for r in rip_order if r not in seen_rips], "Tier 5: Diversity fill", budget=True)

    lib = pd.DataFrame(rows).reset_index(drop=True)
    lib["Selection_tier"] = labels
    print(f"\nDiversity library: {len(lib)} origins, {lib['Origin_span_bp'].sum():,} bp")
    for lbl in sorted(lib["Selection_tier"].unique()):
        sub = lib[lib["Selection_tier"] == lbl]
        print(f"  {lbl}: {len(sub)} origins, {sub['Origin_span_bp'].sum():,} bp")
    return lib


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8  FETCH GENBANK ANNOTATIONS FOR SELECTED ORIGINS
# ─────────────────────────────────────────────────────────────────────────────

def fetch_annotations_for_library(lib: pd.DataFrame,
                                   ncbi_email: str,
                                   ncbi_api_key: str | None = None) -> tuple[dict, dict]:
    """
    Fetch GenBank records for the final selected origins.
    Returns (sequences_dict, records_dict).
    If sequences are already embedded in the DataFrame (Origin_sequence column),
    they are used directly; only the annotation record is fetched.
    """
    Entrez.email = ncbi_email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    sleep   = NCBI_SLEEP_KEY if ncbi_api_key else NCBI_SLEEP
    seqs    = {}
    records = {}

    for i, (_, row) in enumerate(lib.iterrows()):
        acc   = row["Accession"]
        ws    = int(row["Origin_start_bp"])
        we    = int(row["Origin_end_bp"])
        plen  = int(row["Length_bp"])
        loc   = str(row.get("Origin_GenBank_location", ""))

        # Use embedded sequence if available
        embedded = row.get("Origin_sequence", "")
        if isinstance(embedded, str) and len(embedded) >= 100:
            seqs[acc] = embedded
            # Still need GenBank record for annotation
            try:
                fs, fe = max(1, ws - 1000 + 1), min(plen, we + 1000)
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb",
                                       retmode="text", seq_start=fs, seq_stop=fe)
                records[acc] = SeqIO.read(StringIO(handle.read()), "genbank")
                handle.close()
            except Exception:
                pass
        else:
            # Fetch sequence
            try:
                if "join" in loc:
                    handle = Entrez.efetch(db="nucleotide", id=acc,
                                           rettype="gb", retmode="text")
                    rec = SeqIO.read(StringIO(handle.read()), "genbank")
                    handle.close()
                    full = str(rec.seq)
                    seqs[acc]    = full[ws:] + full[:we]
                    records[acc] = rec
                else:
                    fs = max(1, ws - 1000 + 1)
                    fe = min(plen, we + 1000)
                    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb",
                                           retmode="text", seq_start=fs, seq_stop=fe)
                    rec = SeqIO.read(StringIO(handle.read()), "genbank")
                    handle.close()
                    off = ws - (fs - 1)
                    seqs[acc]    = str(rec.seq[off: off + (we - ws)])
                    records[acc] = rec
            except Exception as e:
                print(f"  WARNING: {acc}: {e}")

        if (i + 1) % 10 == 0 or i == 0:
            print(f"  [{i+1}/{len(lib)}] {acc}: {len(seqs.get(acc,''))} bp")
        time.sleep(sleep)

    return seqs, records


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9  FUNCTIONAL OriV VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def parse_genbank_features(records: dict) -> dict[str, list[dict]]:
    gb: dict = {}
    for acc, rec in records.items():
        feats = []
        recs = [rec] if not isinstance(rec, tuple) else list(rec)
        for r in recs:
            if not hasattr(r, "features"):
                continue
            for f in r.features:
                if f.type != "CDS":
                    continue
                product = f.qualifiers.get("product", ["unknown"])[0]
                feats.append({
                    "product":    product,
                    "protein_id": f.qualifiers.get("protein_id", [""])[0],
                    "start":      int(f.location.start),
                    "end":        int(f.location.end),
                    "strand":     f.location.strand,
                    "category":   _cds_category(product),
                })
        gb[acc] = feats
    return gb


def predict_orfs_for_library(seqs: dict) -> dict[str, list[dict]]:
    finder = pyrodigal.GeneFinder(meta=True)
    all_orfs = {}
    for acc, seq in seqs.items():
        genes = finder.find_genes(seq.encode())
        orfs  = []
        for j, g in enumerate(genes):
            aa = g.translate().rstrip("*")
            orfs.append({"id": f"{acc}_orf{j+1}", "start": g.begin, "end": g.end,
                          "strand": "+" if g.strand == 1 else "-", "aa": aa})
        all_orfs[acc] = orfs
    return all_orfs


def run_rep_hmm_on_library(all_orfs: dict, hmm_path: str,
                            evalue_strict: float = 1e-5,
                            evalue_relaxed: float = 0.1) -> dict[str, dict]:
    """Two-pass hmmsearch (strict then relaxed) on the final selected library."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp:
        p = tmp.name
        for acc, orfs in all_orfs.items():
            for orf in orfs:
                if len(orf["aa"]) >= 30:
                    tmp.write(f">{orf['id']}\n{orf['aa']}\n")

    strict  = _run_hmmsearch(p, hmm_path, evalue_strict)
    no_hit  = {acc: orfs for acc, orfs in all_orfs.items() if acc not in strict}
    # Write only no-hit ORFs for relaxed pass
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp2:
        p2 = tmp2.name
        for acc, orfs in no_hit.items():
            for orf in orfs:
                if len(orf["aa"]) >= 30:
                    tmp2.write(f">{orf['id']}\n{orf['aa']}\n")
    relaxed = _run_hmmsearch(p2, hmm_path, evalue_relaxed) if no_hit else {}
    os.unlink(p); os.unlink(p2)
    for f in [p+".tbl", p2+".tbl"]:
        if os.path.exists(f): os.unlink(f)

    results = {}
    for acc, orfs in all_orfs.items():
        if acc in strict:
            hits = strict[acc]; tier = "strict"
        elif acc in relaxed:
            hits = relaxed[acc]; tier = "relaxed"
        else:
            hits = []; tier = "none"
        best = min(hits, key=lambda h: h["evalue"]) if hits else None
        results[acc] = {"hits": hits, "best_evalue": best["evalue"] if best else None,
                        "best_name": best["name"] if best else None, "hmm_tier": tier}
    n_s = sum(1 for v in results.values() if v["hmm_tier"] == "strict")
    n_r = sum(1 for v in results.values() if v["hmm_tier"] == "relaxed")
    print(f"  Validation HMM: {n_s} strict, {n_r} relaxed, {len(results)-n_s-n_r} no-hit")
    return results


def detect_par(acc: str, gb_feats: list, orfs: list) -> dict:
    par_prods = [f["product"] for f in gb_feats if f["category"] == "partition"]
    par_types: set = set()
    for p in par_prods:
        pl = p.lower()
        if any(k in pl for k in {"para", "sopa", "walker", "segreag"}): par_types.add("ParA")
        if any(k in pl for k in {"parb", "sopb", "spo0j", "noc"}):      par_types.add("ParB")
        if any(k in pl for k in {"parg", "ribbon", "rhh", "copg"}):     par_types.add("RHH")
        if not par_types: par_types.add("Par-unknown")
    has_walker = any(len(o["aa"]) >= 100 and WALKER_A_RE.search(o["aa"]) for o in orfs)
    if has_walker and "ParA" not in par_types: par_types.add("ParA-motif")
    return {"par_detected": bool(par_prods) or has_walker,
            "has_par_genbank": bool(par_prods),
            "par_types": par_types,
            "par_products": par_prods,
            "has_para_walker": has_walker}


def detect_iterons(seq: str, min_len: int = 10, max_len: int = 22,
                   min_copies: int = 2, max_spacing: int = 200) -> list:
    seq_u = seq.upper()
    found = []
    for kl in range(min_len, min(max_len + 1, len(seq_u) // 2)):
        pos_map: dict[str, list] = defaultdict(list)
        for i in range(len(seq_u) - kl + 1):
            k = seq_u[i: i + kl]
            if "N" not in k:
                pos_map[k].append(i)
        for kmer, positions in pos_map.items():
            if len(positions) < min_copies:
                continue
            spacings  = [positions[j+1] - positions[j] - kl for j in range(len(positions)-1)]
            clustered = all(s <= max_spacing for s in spacings)
            found.append({"seq": kmer, "copies": len(positions), "length": kl,
                          "positions": positions, "clustered": clustered})
    if not found:
        return []
    found.sort(key=lambda x: (not x["clustered"], -x["length"], -x["copies"]))
    kept = []
    for c in found[:50]:
        if not any(c["seq"] in k["seq"] or k["seq"] in c["seq"] for k in kept):
            kept.append(c)
        if len(kept) >= 5: break
    return kept


def detect_at_rich(seq: str, window: int = 100,
                   threshold: float = 0.65, step: int = 10) -> list:
    seq_u = seq.upper()
    hits  = []
    for i in range(0, len(seq_u) - window + 1, step):
        sub = seq_u[i: i + window]
        at  = (sub.count("A") + sub.count("T")) / len(sub)
        if at >= threshold:
            hits.append((i, i + window, at))
    if not hits: return []
    merged = [list(hits[0])]
    for s, e, at in hits[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e); merged[-1][2] = max(merged[-1][2], at)
        else:
            merged.append([s, e, at])
    return [tuple(r) for r in merged]


def detect_dnaa_boxes(seq: str) -> list:
    seq_u = seq.upper()
    hits  = []
    n     = len(DNAA_CONSENSUS)
    def _h1(a, b): return sum(x != y for x, y in zip(a, b)) <= 1
    for i in range(len(seq_u) - n + 1):
        sub = seq_u[i: i + n]
        if _h1(sub, DNAA_CONSENSUS): hits.append((i, i + n, "fwd"))
        if _h1(sub, DNAA_RC):        hits.append((i, i + n, "rev"))
    return hits


def score_oriv(acc: str, row: pd.Series, seq: str,
               hmm_res: dict, gb_feats: list, orfs: list,
               par_res: dict, iterons: list, at_rich: list, dnaa: list) -> dict:
    """
    Score functional OriV content and classify the origin.

    Scoring (max ~16 pts)
    ─────────────────────
    Rep protein
      +4  HMM strict (E < 1e-5)
      +3  HMM relaxed
      +2  GenBank CDS annotation
      +1  MMseqs2 pre-confirmed (from upstream metadata)

    Par system
      +3  Both ParA + ParB detected (GenBank or motif)
      +2  One Par component (GenBank)
      +1  Walker-A motif only

    OriV structural elements
      +3  Clustered iterons (≥ 2 copies, spacing ≤ 200 bp)
      +1  Non-clustered repeats only
      +2  AT-rich melting region (≥ 65% AT in ≥ 100 bp)
      +1  DnaA box(es)

    Classification
    ──────────────
    COMPLETE    Rep + Par + structural  score ≥ 9
    FUNCTIONAL  Rep + structural        score ≥ 6
    MINIMAL     Rep only               score ≥ 4
    INCOMPLETE  insufficient evidence   score < 4
    """
    sc = 0; ev = []; has_rep = False

    if hmm_res["hmm_tier"] == "strict":  sc += 4; ev.append("HMM-Rep(strict)");  has_rep = True
    elif hmm_res["hmm_tier"] == "relaxed": sc += 3; ev.append("HMM-Rep(relaxed)"); has_rep = True

    rep_cds = [f for f in gb_feats if f["category"] == "rep"]
    if rep_cds: sc += 2; ev.append(f"GenBank-Rep({len(rep_cds)})"); has_rep = True

    if pd.notna(row.get("MMseqs2_evalue")): sc += 1; ev.append("MMseqs2"); has_rep = True
    if pd.notna(row.get("Rep_hmm_evalue")):
        if not has_rep: sc += 2; ev.append("OriVFinder-Rep"); has_rep = True

    # Par
    n_types = len({t for t in par_res["par_types"] if t not in {"Par-unknown", "ParA-motif"}})
    if par_res["has_par_genbank"]:
        if n_types >= 2: sc += 3; ev.append("ParA+ParB-GB")
        else: sc += 2; ev.append(f"Par-GB({','.join(sorted(par_res['par_types']))})")
    elif par_res["has_para_walker"]: sc += 1; ev.append("ParA-WalkerA")

    # Structural
    has_struct = False
    cl_it = [it for it in iterons if it["clustered"]]
    if cl_it:
        best = max(cl_it, key=lambda x: x["copies"])
        sc += 3; ev.append(f"Iterons(clustered,{best['copies']}x{best['length']}bp)"); has_struct = True
    elif iterons:
        sc += 1; ev.append("Repeats(non-clustered)")

    if at_rich:
        best_at = max(at_rich, key=lambda x: x[2])
        sc += 2; ev.append(f"AT-rich({best_at[2]*100:.0f}%)"); has_struct = True

    if dnaa: sc += 1; ev.append(f"DnaA-box({len(dnaa)})")

    if   has_rep and par_res["par_detected"] and has_struct and sc >= 9: cls = "COMPLETE"
    elif has_rep and has_struct and sc >= 6:                              cls = "FUNCTIONAL"
    elif has_rep and sc >= 4:                                             cls = "MINIMAL"
    else:                                                                 cls = "INCOMPLETE"

    return {
        "oriv_score":          sc,
        "oriv_class":          cls,
        "oriv_evidence":       "; ".join(ev),
        "has_rep":             has_rep,
        "has_par":             par_res["par_detected"],
        "par_types":           ", ".join(sorted(par_res["par_types"])),
        "has_par_genbank":     par_res["has_par_genbank"],
        "has_para_walker":     par_res["has_para_walker"],
        "par_products":        "; ".join(par_res["par_products"][:3]),
        "n_iterons":           sum(it["copies"] for it in iterons),
        "n_clustered_iterons": len(cl_it),
        "best_iteron_len":     max((it["length"] for it in iterons), default=0),
        "n_at_rich":           len(at_rich),
        "max_at_content":      max((r[2] for r in at_rich), default=0.0),
        "n_dnaa_boxes":        len(dnaa),
        "hmm_rep_evalue":      hmm_res["best_evalue"],
        "hmm_rep_name":        hmm_res["best_name"],
        "hmm_tier":            hmm_res["hmm_tier"],
        "n_rep_cds_genbank":   len(rep_cds),
        "n_orfs":              len(orfs),
    }


def validate_library(lib: pd.DataFrame, seqs: dict,
                     records: dict, hmm_path: str) -> pd.DataFrame:
    print("\nStep 9a: ORF prediction ...")
    all_orfs = predict_orfs_for_library(seqs)
    print("Step 9b: Rep-protein HMM (library validation) ...")
    hmm_res  = run_rep_hmm_on_library(all_orfs, hmm_path)
    print("Step 9c: Parsing GenBank CDS features ...")
    gb_map   = parse_genbank_features(records)
    print("Step 9d: Scoring Rep, Par, and structural OriV elements ...")

    rows = []
    for _, row in lib.iterrows():
        acc  = row["Accession"]
        seq  = seqs.get(acc, "")
        if not seq:
            rows.append({"Accession": acc, "oriv_class": "INCOMPLETE",
                         "oriv_score": 0, "oriv_evidence": "no-seq"}); continue
        orfs    = all_orfs.get(acc, [])
        gb_f    = gb_map.get(acc, [])
        hr      = hmm_res.get(acc, {"hmm_tier":"none","best_evalue":None,
                                    "best_name":None,"hits":[]})
        par     = detect_par(acc, gb_f, orfs)
        itr     = detect_iterons(seq)
        at      = detect_at_rich(seq)
        dnaa    = detect_dnaa_boxes(seq)
        vrow    = score_oriv(acc, row, seq, hr, gb_f, orfs, par, itr, at, dnaa)
        vrow["Accession"] = acc
        rows.append(vrow)

    val_df = pd.DataFrame(rows)
    for cls in ["COMPLETE","FUNCTIONAL","MINIMAL","INCOMPLETE"]:
        print(f"  {cls:12s}: {(val_df['oriv_class']==cls).sum()}/{len(val_df)}")
    return val_df


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10  SYNTHESIS ASSESSMENT
# ─────────────────────────────────────────────────────────────────────────────

def assess_synthesis(lib: pd.DataFrame, seqs: dict | None = None) -> pd.DataFrame:
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
        gp = gc*100
        gs = 50 if 40<=gp<=60 else 40 if 30<=gp<=65 else 25 if 25<=gp<=70 else 10 if 20<=gp<=75 else 0
        ls = 50 if bp<=1800 else 40 if bp<=3000 else 30 if bp<=5000 else 20 if bp<=6000 else 10
        return gs + ls
    def _ease_cat(s):
        return "Easy" if s>=80 else "Moderate" if s>=60 else "Challenging" if s>=40 else "Difficult"

    lib = lib.copy()
    lib["Twist_GC_cat"]      = lib["GC_content"].apply(_gc_cat)
    lib["Twist_strategy"]    = lib["Origin_span_bp"].apply(_strat)
    lib["Twist_n_fragments"] = lib["Origin_span_bp"].apply(_frags)
    lib["Twist_est_cost"]    = lib.apply(lambda r: _cost(r["Origin_span_bp"], r["GC_content"]), axis=1)
    lib["Twist_ease_score"]  = lib.apply(lambda r: _ease(r["GC_content"], r["Origin_span_bp"]), axis=1)
    lib["Twist_ease_cat"]    = lib["Twist_ease_score"].apply(_ease_cat)
    if seqs:
        lib["Twist_has_homopolymer"] = lib["Accession"].map(
            lambda a: bool(re.search(r"([ACGT])\1{10,}", seqs.get(a,"").upper())))
    return lib


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11  OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

def write_fasta(seqs: dict, path: str):
    with open(path, "w") as fh:
        for acc, seq in seqs.items():
            fh.write(f">{acc}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    print(f"FASTA → {path}  ({len(seqs)} sequences)")


def draw_origin_map(ax, acc, row, seq, gb_feats, orfs, val, iterons, at_rich):
    span = len(seq)
    ax.set_xlim(-span*0.05, span*1.15); ax.set_ylim(-3.5, 5.5); ax.axis("off")
    ax.plot([0,span],[0,0],"k-",lw=3,solid_capstyle="butt")
    tick = max(500, (span//10//500)*500)
    for p in range(0, span+1, tick):
        ax.plot([p,p],[-0.2,0.2],"k-",lw=0.8)
        ax.text(p,-0.5,f"{p/1000:.1f}kb",ha="center",fontsize=6,color="gray")
    for s,e,at in at_rich:
        ax.axvspan(s,e,ymin=0.42,ymax=0.58,color="#FFF176",alpha=0.8,zorder=1)
        ax.text((s+e)/2,0.28,f"AT\n{at*100:.0f}%",ha="center",fontsize=5,color="#F9A825")

    y_fwd, y_rev, ah = 0.9, -1.2, 0.55
    def _arrow(st, en, strand, color, label):
        y = y_fwd if strand >= 0 else y_rev
        hd = min(span*0.025, abs(en-st)*0.25, 300)
        if strand >= 0:
            body = max(st, en-hd)
            ax.fill_between([st,body],y-ah/2,y+ah/2,color=color,alpha=0.75,lw=0.5,edgecolor="k")
            ax.fill([body,en,body],[y+ah/2,y,y-ah/2],color=color,alpha=0.75,lw=0.5,edgecolor="k")
            lx,ly,va=(st+en)/2,y+ah/2+0.18,"bottom"
        else:
            body=min(en,st+hd)
            ax.fill_between([body,en],y-ah/2,y+ah/2,color=color,alpha=0.75,lw=0.5,edgecolor="k")
            ax.fill([body,st,body],[y+ah/2,y,y-ah/2],color=color,alpha=0.75,lw=0.5,edgecolor="k")
            lx,ly,va=(st+en)/2,y-ah/2-0.15,"top"
        short = label[:28]+"…" if len(label)>28 else label
        ax.text(lx,ly,short,ha="center",va=va,fontsize=5.5,style="italic",
                bbox=dict(boxstyle="round,pad=0.1",fc="white",alpha=0.7,ec="none"),clip_on=True)

    drawn = False
    if gb_feats:
        for f in gb_feats:
            _arrow(f["start"],f["end"],f["strand"] or 1,
                   FEATURE_COLORS.get(f["category"],FEATURE_COLORS["other"]),f["product"])
        drawn = True
    if not drawn:
        for orf in orfs:
            _arrow(orf["start"],orf["end"],1 if orf["strand"]=="+" else -1,
                   FEATURE_COLORS["other"],f"ORF({len(orf['aa'])}aa)")

    win = max(50, span//100); gc_p, gc_v = [],[]
    for i in range(0, span-win+1, win//2):
        sub = seq[i:i+win].upper()
        gc_v.append((sub.count("G")+sub.count("C"))/len(sub)); gc_p.append(i+win/2)
    if gc_v:
        ga = np.array(gc_v); rng = ga.max()-ga.min()+1e-6
        yg = -2.5+(ga-ga.min())/rng
        ax.fill_between(gc_p,-2.5,yg,color="#90CAF9",alpha=0.45)
        ax.plot(gc_p,yg,color="#1565C0",lw=0.7)
        ax.text(-span*0.04,-2.2,"GC",ha="right",fontsize=6,color="#1565C0")

    cls = val.get("oriv_class","?")
    cc  = {"COMPLETE":"#2E7D32","FUNCTIONAL":"#1565C0","MINIMAL":"#F9A825","INCOMPLETE":"#C62828"}.get(cls,"gray")
    ax.set_title(f"{acc}  |  {row.get('Origin_RIP_type','?')}  |  {row.get('Species','?')}",
                 fontsize=9,fontweight="bold",pad=4)
    par_t = val.get("par_types","") or "—"
    meta  = (f"Span: {span:,} bp  |  GC: {row.get('GC_content',0)*100:.1f}%  |  "
             f"Tier: {str(row.get('Selection_tier','?')).split(':')[-1].strip()}\n"
             f"OriV: {cls} (score={val.get('oriv_score','?')})  |  "
             f"Par: {'Yes' if val.get('has_par') else 'No'} ({par_t})  |  "
             f"Evidence: {val.get('oriv_evidence','?')}")
    ax.text(span*0.5,4.9,meta,ha="center",va="top",fontsize=6,color=cc,
            bbox=dict(boxstyle="round",fc="#F1F8E9",alpha=0.85,ec=cc,lw=1.2))
    hdls = [mpatches.Patch(color=c,alpha=0.75,label=k.replace("_"," ").title())
            for k,c in FEATURE_COLORS.items()]
    hdls += [mpatches.Patch(color="#FFF176",alpha=0.8,label="AT-rich")]
    ax.legend(handles=hdls,loc="upper right",fontsize=5,ncol=2,framealpha=0.8,edgecolor="gray")


def generate_pdf(lib, seqs, records, all_orfs, gb_map, val_df, path):
    val_d = val_df.set_index("Accession").to_dict("index")
    n     = len(lib)
    with PdfPages(path) as pdf:
        # Title page
        fig,ax = plt.subplots(figsize=(11,8.5)); ax.axis("off")
        ax.text(0.5,0.60,"Plasmid Origin Diversity Library",ha="center",fontsize=26,
                fontweight="bold",transform=ax.transAxes)
        ax.text(0.5,0.45,f"{n} origins  |  {lib['Origin_span_bp'].sum():,} bp total",
                ha="center",fontsize=16,transform=ax.transAxes)
        c_stat = " | ".join(f"{cls}: {(val_df['oriv_class']==cls).sum()}"
                            for cls in ["COMPLETE","FUNCTIONAL","MINIMAL","INCOMPLETE"])
        ax.text(0.5,0.30,c_stat,ha="center",fontsize=12,transform=ax.transAxes,color="#1B5E20")
        pdf.savefig(fig,bbox_inches="tight"); plt.close(fig)

        for idx,(_, row) in enumerate(lib.iterrows()):
            acc = row["Accession"]
            seq = seqs.get(acc,"")
            if not seq: continue
            try:
                fig,ax = plt.subplots(figsize=(14,5))
                draw_origin_map(ax,acc,row,seq,
                                gb_map.get(acc,[]),all_orfs.get(acc,[]),
                                val_d.get(acc,{}),detect_iterons(seq),detect_at_rich(seq))
                plt.tight_layout(); pdf.savefig(fig,bbox_inches="tight")
            except Exception as e:
                fig,ax = plt.subplots(figsize=(14,5))
                ax.text(0.5,0.5,f"{acc}: {e}",ha="center",va="center",transform=ax.transAxes,color="red")
                pdf.savefig(fig,bbox_inches="tight")
            finally:
                plt.close(fig)
            if (idx+1) % 20 == 0: print(f"  Maps: {idx+1}/{n}")
    print(f"PDF → {path}")


def generate_overview(lib, val_df, path):
    fig  = plt.figure(figsize=(22,24))
    gs_  = gridspec.GridSpec(4,3,hspace=0.38,wspace=0.32,
                             left=0.06,right=0.96,top=0.95,bottom=0.04)
    tp   = [t for t in TIER_COLORS if t in lib["Selection_tier"].values]
    tc   = [TIER_COLORS[t] for t in tp]

    def _ax(pos): return fig.add_subplot(gs_[pos])

    ax = _ax((0,0))
    cls_o = ["COMPLETE","FUNCTIONAL","MINIMAL","INCOMPLETE"]
    cls_c = ["#2E7D32","#1565C0","#F9A825","#C62828"]
    cnts  = [(val_df["oriv_class"]==c).sum() for c in cls_o]
    bars = ax.bar(cls_o,cnts,color=cls_c,edgecolor="white",lw=1.5)
    for b,c in zip(bars,cnts): ax.text(b.get_x()+b.get_width()/2,b.get_height()+0.3,str(c),ha="center",fontweight="bold")
    ax.set_title("A. OriV Classification",fontweight="bold",fontsize=12)
    ax.set_ylim(0,max(cnts)*1.25)

    ax = _ax((0,1))
    py = val_df["has_par"].sum(); pn = len(val_df)-py
    ax.pie([py,pn],labels=[f"Par\n({py})",f"No Par\n({pn})"],colors=["#1565C0","#EF9A9A"],
           autopct="%1.0f%%",startangle=90,textprops={"fontsize":10})
    ax.set_title("B. Par System",fontweight="bold",fontsize=12)

    ax = _ax((0,2))
    ax.hist(val_df["oriv_score"],bins=range(0,int(val_df["oriv_score"].max())+2),
            color="#1976D2",edgecolor="white",alpha=0.85)
    for th,lbl,col in [(4,"MINIMAL","#F9A825"),(6,"FUNCTIONAL","#1565C0"),(9,"COMPLETE","#2E7D32")]:
        ax.axvline(th,color=col,ls="--",lw=1.5,label=lbl)
    ax.legend(fontsize=7); ax.set_xlabel("OriV Score"); ax.set_title("C. Score Distribution",fontweight="bold",fontsize=12)

    ax = _ax((1,0))
    en = ["Iterons\n(clustered)","AT-rich","DnaA box","Rep HMM\n(strict)","Par\n(GenBank)"]
    ec = [(val_df["n_clustered_iterons"]>0).sum(),(val_df["n_at_rich"]>0).sum(),
          (val_df["n_dnaa_boxes"]>0).sum(),(val_df["hmm_tier"]=="strict").sum(),
          val_df["has_par_genbank"].sum()]
    ax.barh(en,ec,color=["#FF8F00","#FDD835","#26A69A","#E53935","#1565C0"],edgecolor="white")
    for i,c in enumerate(ec): ax.text(c+0.3,i,f"{c}/{len(val_df)}",va="center",fontsize=9)
    ax.set_xlim(0,len(val_df)*1.25)
    ax.set_title("D. Structural Elements",fontweight="bold",fontsize=12)

    ax = _ax((1,1))
    tc2   = lib["Selection_tier"].value_counts()
    tbp   = lib.groupby("Selection_tier")["Origin_span_bp"].sum()
    bars  = ax.barh(range(len(tp)),[tc2.get(t,0) for t in tp],color=tc,edgecolor="white",height=0.6)
    for b,t in zip(bars,tp):
        ax.text(b.get_width()+0.3,b.get_y()+b.get_height()/2,
                f"{tc2.get(t,0)}  ({tbp.get(t,0)/1000:.0f} kb)",va="center",fontsize=8)
    ax.set_yticks(range(len(tp))); ax.set_yticklabels([t.split(":")[-1].strip() for t in tp],fontsize=9)
    ax.invert_yaxis(); ax.set_title("E. Tier Composition",fontweight="bold",fontsize=12)

    ax = _ax((1,2))
    for t in tp:
        sub = lib[lib["Selection_tier"]==t]
        ax.scatter(sub["GC_content"],np.random.uniform(0.3,1.0,len(sub)),
                   color=TIER_COLORS[t],s=30,alpha=0.8,label=t.split(":")[0])
    ax.axvline(AGRO_GC,color="red",ls="--",lw=1.5,label=f"Agrobacterium ({AGRO_GC*100:.0f}%)")
    ax.axvspan(0.25,0.65,alpha=0.05,color="green"); ax.legend(fontsize=7)
    ax.set_xlabel("GC Content"); ax.set_title("F. GC Content",fontweight="bold",fontsize=12)

    ax = _ax((2,0))
    sd = [lib[lib["Selection_tier"]==t]["Origin_span_bp"].values for t in tp]
    bp_= ax.boxplot(sd,vert=True,patch_artist=True,labels=[t.split(":")[-1].strip() for t in tp])
    for p,c in zip(bp_["boxes"],tc): p.set_facecolor(c); p.set_alpha(0.7)
    ax.set_ylabel("bp"); ax.tick_params(axis="x",rotation=25)
    ax.set_title("G. Origin Size",fontweight="bold",fontsize=12)

    ax = _ax((2,1))
    pc = lib["Phylum"].value_counts()
    ax.barh(range(len(pc)),pc.values,color="#2196F3",edgecolor="white")
    ax.set_yticks(range(len(pc))); ax.set_yticklabels(pc.index,fontsize=8)
    for i,v in enumerate(pc.values): ax.text(v+0.1,i,str(v),va="center",fontsize=8)
    ax.invert_yaxis(); ax.set_title("H. Phylum Coverage",fontweight="bold",fontsize=12)

    ax = _ax((2,2))
    eo = ["Easy","Moderate","Challenging","Difficult"]
    ec_ = ["#4CAF50","#FFC107","#FF9800","#F44336"]
    ecs = [(lib["Twist_ease_cat"]==e).sum() for e in eo]
    bs  = ax.bar(eo,ecs,color=ec_,edgecolor="white")
    for b,c in zip(bs,ecs):
        if c: ax.text(b.get_x()+b.get_width()/2,b.get_height()+0.3,str(c),ha="center",fontweight="bold")
    ax.set_title("I. Twist Synthesis Ease",fontweight="bold",fontsize=12)

    ax = _ax((3,slice(None)))
    ax.axis("off")
    summary = (
        f"LIBRARY:  {len(lib)} origins  |  {lib['Origin_span_bp'].sum():,} bp  |  "
        f"{lib['Origin_RIP_type'].nunique()} RIP types  |  "
        f"{lib['Phylum'].nunique()} phyla  |  {lib['Genus'].nunique()} genera\n"
        f"OriV:  Complete {(val_df['oriv_class']=='COMPLETE').sum()}  "
        f"Functional {(val_df['oriv_class']=='FUNCTIONAL').sum()}  "
        f"Minimal {(val_df['oriv_class']=='MINIMAL').sum()}  "
        f"Incomplete {(val_df['oriv_class']=='INCOMPLETE').sum()}  |  "
        f"Par system: {val_df['has_par'].sum()}/{len(val_df)}\n"
        f"SYNTHESIS:  Est. total ${lib['Twist_est_cost'].sum():,.0f}  |  "
        f"{lib['Twist_n_fragments'].sum()} fragments  |  "
        f"{(lib['Twist_ease_cat']=='Easy').sum()} easy  "
        f"{(lib['Twist_ease_cat']=='Moderate').sum()} moderate  "
        f"{(lib['Twist_ease_cat'].isin(['Challenging','Difficult'])).sum()} challenging"
    )
    ax.text(0.5,0.5,summary,ha="center",va="center",fontsize=11,transform=ax.transAxes,
            fontfamily="monospace",bbox=dict(boxstyle="round",fc="#E3F2FD",alpha=0.85))

    fig.suptitle("Plasmid Origin Diversity Library — Summary",
                 fontsize=18,fontweight="bold",y=0.98)
    fig.savefig(path,dpi=150,bbox_inches="tight"); plt.close(fig)
    print(f"Overview → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# PLSDB FULL MODE  (PLSDB 2025 + OriVFinder scoring)
# ─────────────────────────────────────────────────────────────────────────────

def run_plsdb_full(args, out: Path):
    """
    Full PLSDB 2025 pipeline with integrated OriVFinder IGS scoring.

    Stages:
      S1. Load and join PLSDB 2025 metadata (nuccore + taxonomy + typing)
      S2. Metadata filter (circular, size, not repABC, not NHR E. coli)
      S3. Stream filtered sequences from sequences.fasta
      S4. Batch ORF prediction (pyrodigal) + Rep-protein HMM (hmmsearch)
      S5. OriVFinder IGS scoring: extract IGSs near each Rep gene,
          score with Iteron + ZCurve + PatternFinder, pick best IGS per plasmid
      S6. Define origin windows with ≥200 bp Rep-buffer (≤6 kb)
      S7. Fetch GenBank annotations + trim non-origin CDS from edges
      S8. Deduplicate by RIP type
      S9. Tiered diversity selection
      S10. Functional OriV validation
      S11. Synthesis assessment
      S12. Outputs
    """
    sep = "="*70
    Entrez.email = args.email
    if args.ncbi_api_key:
        Entrez.api_key = args.ncbi_api_key

    # S1 — Load PLSDB 2025
    print(f"\n{sep}\nS1: Load PLSDB 2025 metadata\n{sep}", flush=True)
    plsdb_dir = Path(args.plsdb_dir)
    df = load_plsdb_2025(str(plsdb_dir))

    # S2 — Metadata filter
    print(f"\n{sep}\nS2: Metadata filter\n{sep}", flush=True)
    df = metadata_filter(df)
    print(f"After metadata filter: {len(df):,} plasmids", flush=True)

    # S3 — Load sequences from PLSDB FASTA
    print(f"\n{sep}\nS3: Load plasmid sequences from PLSDB FASTA\n{sep}", flush=True)
    fasta_path = plsdb_dir / "sequences.fasta"
    if not fasta_path.exists():
        # try bz2
        bz2_path = plsdb_dir / "sequences.fasta.bz2"
        if bz2_path.exists():
            print("  Decompressing sequences.fasta.bz2 ...", flush=True)
            import bz2 as _bz2
            with _bz2.open(str(bz2_path), "rb") as fi, \
                 open(str(fasta_path), "wb") as fo:
                fo.write(fi.read())
        else:
            raise FileNotFoundError(
                f"PLSDB FASTA not found. Expected {fasta_path} or {bz2_path}")

    # ── Streaming batch pipeline: ORF prediction + hmmsearch + OriVFinder scoring ──
    # Process sequences in batches of BATCH_SIZE to keep RAM < 1GB.
    # For each batch:
    #   1. Predict ORFs (pyrodigal)
    #   2. Write proteins to temp FASTA
    #   3. Run hmmsearch against RIP.hmm
    #   4. For hits: run OriVFinder IGS scoring
    #   5. Accumulate results; discard sequences/ORFs from memory
    # ─────────────────────────────────────────────────────────────────────────────
    STREAM_BATCH = 2000   # sequences per batch — tunes RAM vs. hmmsearch overhead
    N_WORKERS    = max(1, (os.cpu_count() or 4) - 1)  # leave one core for OS
    acc_needed   = set(df["Accession"].tolist())
    meta_lookup  = df.set_index("Accession")
    oriv_rows: list[dict] = []

    print(f"\n{sep}\nS3+S4+S5: Streaming ORF prediction + HMM + OriVFinder (batch={STREAM_BATCH})\n{sep}",
          flush=True)
    if not _SCIPY_OK:
        print("  WARNING: scipy not available — ZCurve scoring disabled", flush=True)

    def _hits_are_repabc(hits: list) -> bool:
        names = " ".join(h["name"].lower() for h in hits[:5])
        return any(k in names for k in REPABC_KW)

    batch_seqs: dict[str, str] = {}
    total_loaded = 0
    total_hits   = 0

    # ThreadPoolExecutor — pyrodigal releases the GIL during find_genes(),
    # so threads get true parallelism with zero pickling/IPC overhead.
    _executor = concurrent.futures.ThreadPoolExecutor(max_workers=N_WORKERS)

    def _flush_batch(b_seqs: dict) -> None:
        nonlocal total_hits
        if not b_seqs:
            return
        t0 = time.time()
        # Parallel ORF prediction via threads (no pickling, no IPC)
        b_orfs = dict(_executor.map(_pyrodigal_worker, b_seqs.items()))
        t1 = time.time()
        total_orfs = sum(len(v) for v in b_orfs.values())
        print(f"    ORF prediction: {t1-t0:.1f}s  ({total_orfs:,} ORFs total)", flush=True)
        # Write proteins for this batch
        with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                         delete=False, prefix="plsdb_batch_") as tmp:
            prot_path = tmp.name
            for acc2, orfs2 in b_orfs.items():
                for orf2 in orfs2:
                    tmp.write(f">{orf2['id']}\n{orf2['aa']}\n")
        t2 = time.time()
        print(f"    Write proteins: {t2-t1:.1f}s  ({os.path.getsize(prot_path)//1024:,} KB)", flush=True)
        raw = _run_hmmsearch(prot_path, args.rip_hmm, args.hmm_evalue)
        t3 = time.time()
        print(f"    hmmsearch: {t3-t2:.1f}s  ({len(raw):,} seqs with hits)", flush=True)
        os.unlink(prot_path)
        tbl = prot_path + ".tbl"
        if os.path.exists(tbl): os.unlink(tbl)

        n_scored = 0
        for acc2, orfs2 in b_orfs.items():
            hits = raw.get(acc2, [])
            if not hits or _hits_are_repabc(hits):
                continue
            best = min(hits, key=lambda h: h["evalue"])
            hmm_info = {"hits": hits, "best_evalue": best["evalue"],
                        "best_name": best["name"], "hmm_tier": "strict",
                        "orf_id": best["orf_id"]}
            seq2 = b_seqs.get(acc2, "")
            if not seq2:
                continue
            n_scored += 1
            result = find_best_oriv(acc2, seq2, orfs2, {acc2: hmm_info})
            if result is None:
                continue
            if acc2 not in meta_lookup.index:
                continue
            rd = meta_lookup.loc[acc2].to_dict()
            rd.update(result)
            rd.setdefault("Origin_RIP_type", result.get("Rep_hmm_name", ""))
            rd["Plasmid_seq"] = seq2
            oriv_rows.append(rd)
            total_hits += 1
        t4 = time.time()
        print(f"    OriVFinder scoring: {t4-t3:.1f}s  ({n_scored} scored → {total_hits} total hits)",
              flush=True)

    print(f"  Using {N_WORKERS} parallel workers for ORF prediction (fork pool)", flush=True)
    fasta_size = fasta_path.stat().st_size
    fasta_read  = 0
    with open(str(fasta_path)) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            fasta_read += len(rec.seq) + len(rec.id) + 2   # rough byte estimate
            acc = rec.id.split()[0]
            if acc not in acc_needed:
                acc_nv = acc.rsplit(".", 1)[0]
                if acc_nv in acc_needed:
                    acc = acc_nv
                else:
                    continue
            batch_seqs[acc] = str(rec.seq).upper()
            total_loaded += 1

            if total_loaded % 500 == 0:
                pct = 100 * fasta_read / fasta_size
                print(f"  reading: {total_loaded:,} sequences accepted ({pct:.1f}% of FASTA)", flush=True)

            if len(batch_seqs) >= STREAM_BATCH:
                print(f"  flushing batch at {total_loaded:,} seqs …", flush=True)
                _flush_batch(batch_seqs)
                batch_seqs.clear()
                print(f"  [{total_loaded:,} seqs / {total_hits:,} origins so far]", flush=True)

    # Flush remaining
    if batch_seqs:
        _flush_batch(batch_seqs)
        batch_seqs.clear()

    _executor.shutdown(wait=False)

    print(f"  Loaded {total_loaded:,} sequences; identified {total_hits:,} origin candidates",
          flush=True)

    oriv_df = pd.DataFrame(oriv_rows)
    print(f"  Origins identified: {len(oriv_df):,}", flush=True)

    # S6 — Origin windows (200 bp buffer already enforced in find_best_oriv)
    print(f"\n{sep}\nS6: Origin windows (≥{REP_BUFFER_BP} bp Rep-buffer, ≤{ORIGIN_MAX_BP} bp)\n{sep}",
          flush=True)
    oriv_df = oriv_df[oriv_df["Origin_span_bp"].between(200, ORIGIN_MAX_BP)].copy()
    oriv_df["Origin_sequence"] = oriv_df.apply(
        lambda r: r["Plasmid_seq"][int(r["Origin_start_bp"]): int(r["Origin_end_bp"])]
                  if r.get("Plasmid_seq") else "", axis=1)
    oriv_df["Origin_GenBank_location"] = (
        oriv_df["Origin_start_bp"].astype(str) + ".." + oriv_df["Origin_end_bp"].astype(str))
    print(f"  After size filter: {len(oriv_df):,} origins", flush=True)

    # S7 — Local ORF-based CDS trimming (no NCBI fetch — uses pyrodigal ORFs)
    # We already have full plasmid sequences and Rep ORF coordinates from S3–S5.
    # Re-predict ORFs on the origin window region and trim non-Rep ORFs from edges.
    print(f"\n{sep}\nS7: Local ORF-based origin trimming (no GenBank fetch)\n{sep}", flush=True)
    trimmed_rows = []
    s7_finder = pyrodigal.GeneFinder(meta=True)

    for i, (_, row) in enumerate(oriv_df.iterrows()):
        acc  = row["Accession"]
        ws   = int(row["Origin_start_bp"])
        we   = int(row["Origin_end_bp"])
        rs   = int(row.get("Rep_orf_start", ws + REP_BUFFER_BP))
        re_  = int(row.get("Rep_orf_end",   we - REP_BUFFER_BP))
        full_seq = row.get("Plasmid_seq", "")
        seq  = row.get("Origin_sequence", "")

        if (i+1) % 5000 == 0 or i == 0:
            print(f"  [{i+1}/{len(oriv_df)}] {acc}", flush=True)

        # Predict ORFs in a slightly extended window for edge trimming
        ext_s = max(0, ws - 500)
        ext_e = min(len(full_seq), we + 500) if full_seq else we + 500
        ext_seq = full_seq[ext_s:ext_e] if full_seq else ""
        orf_features = []
        if ext_seq and len(ext_seq) >= 60:
            try:
                genes = s7_finder.find_genes(ext_seq.encode())
                for g in genes:
                    aa = g.translate().rstrip("*")
                    if len(aa) < 30:
                        continue
                    orf_features.append({
                        "product":  "predicted_orf",
                        "start":    g.begin + ext_s,
                        "end":      g.end   + ext_s,
                        "strand":   g.strand,
                        "category": "other",
                    })
            except Exception:
                pass

        # Mark ORFs overlapping the Rep gene as "rep"
        for f in orf_features:
            if f["start"] < re_ and f["end"] > rs:
                f["category"] = "rep"
                f["product"]  = row.get("Rep_hmm_name", "rep_protein")

        new_ws, new_we = trim_origin_window(ws, we, orf_features, rs, re_)
        if seq:
            rel_s, rel_e = new_ws - ws, new_we - ws
            new_seq = seq[max(0, rel_s): max(0, rel_e)]
        else:
            new_seq = seq

        rd = row.to_dict()
        rd.update({"Origin_start_bp": new_ws, "Origin_end_bp": new_we,
                   "Origin_span_bp": new_we - new_ws,
                   "Origin_sequence": new_seq,
                   "Origin_elements": row.get("Rep_hmm_name", ""),
                   "n_rep_elements": 1,
                   "n_par_elements": 0,
                   "has_par": False})
        trimmed_rows.append(rd)

    oriv_df = pd.DataFrame(trimmed_rows)
    oriv_df = oriv_df[oriv_df["Origin_span_bp"].between(200, ORIGIN_MAX_BP)].reset_index(drop=True)
    seqs = {r["Accession"]: r.get("Origin_sequence", "") for _, r in oriv_df.iterrows()}
    seqs = {k: v for k, v in seqs.items() if v}
    print(f"After trimming: {len(oriv_df):,} origins", flush=True)

    # S8 — Deduplicate by RIP type (keep up to 5 per type for larger candidate pool)
    print(f"\n{sep}\nS8: Deduplicate by RIP type\n{sep}", flush=True)
    deduped = deduplicate_and_build_metadata(oriv_df, max_per_type=5)

    # S9 — Diversity selection
    # Overshoot target_bp by ~50% to compensate for Pfam-based origin trimming
    # (trimming typically removes ~30-40% of bp from padded 6kb windows)
    overshoot_bp = int(args.target_bp * 1.5)
    overshoot_n  = int(args.min_origins * 1.5)
    print(f"\n{sep}\nS9: Tiered diversity selection (overshooting to ~{overshoot_bp//1000}kb / "
          f"~{overshoot_n} origins to compensate for post-selection trimming)\n{sep}", flush=True)
    rip_stats = compute_rip_stats(deduped)
    lib = select_diversity_library(deduped, rip_stats, overshoot_bp, overshoot_n)

    # S10 — Validation + ORF classification for maps
    print(f"\n{sep}\nS10: Functional OriV validation + ORF classification\n{sep}", flush=True)
    lib_seqs = {acc: seqs[acc] for acc in lib["Accession"] if acc in seqs}
    records  = {}   # no GenBank records in PLSDB mode
    lib_recs = {}
    val_df   = validate_library(lib, lib_seqs, lib_recs, args.rip_hmm)
    all_orfs_lib = predict_orfs_for_library(lib_seqs)

    # Classify ORFs by running hmmsearch on origin proteins to get Rep hits,
    # then label remaining ORFs by size/position heuristics
    print("  Classifying ORFs for map labels ...", flush=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                     delete=False, prefix="lib_orfs_") as tmp:
        lib_prot_path = tmp.name
        for acc, orfs_list in all_orfs_lib.items():
            for orf in orfs_list:
                tmp.write(f">{orf['id']}\n{orf['aa']}\n")
    lib_hmm_raw = _run_hmmsearch(lib_prot_path, args.rip_hmm, args.hmm_evalue)
    os.unlink(lib_prot_path)
    tbl = lib_prot_path + ".tbl"
    if os.path.exists(tbl):
        os.unlink(tbl)

    # Build gb_map_lib from ORF predictions + HMM classification
    # Map each orf_id to its HMM hit name (if any)
    orf_hmm_names: dict[str, str] = {}
    for acc_hits in lib_hmm_raw.values():
        for hit in acc_hits:
            orf_hmm_names[hit["orf_id"]] = hit["name"]

    gb_map_lib: dict[str, list[dict]] = {}
    for acc, orfs_list in all_orfs_lib.items():
        feats = []
        for orf in orfs_list:
            oid = orf["id"]
            if oid in orf_hmm_names:
                hmm_name = orf_hmm_names[oid]
                category = _cds_category(hmm_name)
                if category == "other":
                    category = "rep"  # HMM hits are Rep proteins
                product = hmm_name
            else:
                # Classify by size heuristic
                aa_len = len(orf.get("aa", ""))
                if aa_len < 80:
                    product = f"small ORF ({aa_len} aa)"
                    category = "hypothetical"
                elif aa_len < 200:
                    product = f"ORF ({aa_len} aa)"
                    category = "hypothetical"
                else:
                    product = f"ORF ({aa_len} aa)"
                    category = "other"
            feats.append({
                "product":  product,
                "start":    orf["start"],
                "end":      orf["end"],
                "strand":   1 if orf["strand"] == "+" else -1,
                "category": category,
            })
        gb_map_lib[acc] = feats
    print(f"  Classified {sum(len(v) for v in gb_map_lib.values()):,} ORFs across {len(gb_map_lib)} origins",
          flush=True)

    # S11 — Synthesis
    print(f"\n{sep}\nS11: Synthesis assessment\n{sep}", flush=True)
    lib = assess_synthesis(lib, lib_seqs)

    # S12 — Outputs
    print(f"\n{sep}\nS12: Writing outputs\n{sep}", flush=True)
    final = lib.merge(val_df, on="Accession", how="left")
    drop  = [c for c in final.columns if c.startswith("_") or c == "Plasmid_seq"]
    final = final.drop(columns=drop, errors="ignore")

    csv_path = out / "diversity_library_metadata.csv"
    final.to_csv(str(csv_path), index=False)
    print(f"CSV → {csv_path}  ({final.shape[0]} × {final.shape[1]})", flush=True)

    fasta_path_out = out / "origin_sequences.fasta"
    write_fasta(lib_seqs, str(fasta_path_out))

    generate_pdf(lib, lib_seqs, lib_recs, all_orfs_lib, gb_map_lib, val_df,
                 str(out / "origin_plasmid_maps.pdf"))
    generate_overview(lib, val_df, str(out / "diversity_library_overview.png"))

    print(f"\n{sep}\nDONE\n{sep}", flush=True)
    for cls in ["COMPLETE", "FUNCTIONAL", "MINIMAL", "INCOMPLETE"]:
        print(f"  {cls:12s} {(val_df['oriv_class']==cls).sum():3d}/{len(val_df)}")


def parse_args():
    p = argparse.ArgumentParser(
        description="Plasmid origin diversity library — from full PLSDB to validated library",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ── Input modes (mutually exclusive) ──────────────────────────────────────
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument("--plsdb-dir",
                     help="PLSDB 2025 directory containing nuccore.csv, taxonomy.csv, "
                          "typing.csv, and sequences.fasta (or sequences.fasta.bz2). "
                          "Runs full pipeline with OriVFinder IGS scoring.")
    inp.add_argument("--plsdb-meta",
                     help="Single PLSDB metadata TSV/CSV. "
                          "Runs the original complete pipeline including OriV finding.")
    inp.add_argument("--prefiltered",
                     help="Pre-filtered origins CSV with Origin_start_bp / Origin_end_bp "
                          "already defined. Skips OriV finding; applies buffer checks, "
                          "trimming, validation and outputs only.")

    p.add_argument("--plsdb-fasta",  default=None,
                   help="PLSDB sequence FASTA for local batch OriV finding (plsdb-meta mode). "
                        "If omitted, sequences are fetched from NCBI (slow).")
    p.add_argument("--seqs-fasta",   default=None,
                   help="Pre-fetched origin sequences FASTA (prefiltered mode). "
                        "If omitted, sequences are fetched from NCBI.")
    p.add_argument("--rip-hmm",      required=True,  help="Path to RIP.hmm")
    p.add_argument("--email",        required=True,  help="Email for NCBI Entrez")
    p.add_argument("--output-dir",   default="results", help="Output directory")
    p.add_argument("--target-bp",    type=int, default=500_000)
    p.add_argument("--min-origins",  type=int, default=100)
    p.add_argument("--hmm-evalue",   type=float, default=1e-5)
    p.add_argument("--ncbi-api-key", default=None)
    p.add_argument("--batch-size",   type=int, default=500)
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# PREFILTERED MODE  (origins already have coordinates; apply buffer + trim + validate)
# ─────────────────────────────────────────────────────────────────────────────

def _enforce_rep_buffer(acc: str, ws: int, we: int, plen: int,
                        origin_seq: str, all_orfs: dict,
                        hmm_results: dict,
                        ncbi_email: str, ncbi_api_key: str | None,
                        buffer: int = REP_BUFFER_BP) -> tuple[int, int, str]:
    """
    Check whether the best Rep ORF has ≥ buffer bp from each window edge.
    If not, extend the window by fetching additional flanking sequence from NCBI
    and return the new (w_start, w_end, origin_seq).
    """
    hit  = hmm_results.get(acc, {})
    if hit.get("hmm_tier", "none") == "none":
        return ws, we, origin_seq   # no Rep found, nothing to buffer

    # Find the Rep ORF coordinates within the current origin_seq
    best_orf_id = None
    if hit.get("hits"):
        best_orf_id = min(hit["hits"], key=lambda h: h["evalue"])["orf_id"]

    rep_start_in_win = rep_end_in_win = None
    for orf in all_orfs.get(acc, []):
        if orf["id"] == best_orf_id:
            rep_start_in_win = orf["start"]
            rep_end_in_win   = orf["end"]
            break

    if rep_start_in_win is None:
        return ws, we, origin_seq

    need_left  = max(0, buffer - rep_start_in_win)
    need_right = max(0, buffer - (len(origin_seq) - rep_end_in_win))

    if need_left == 0 and need_right == 0:
        return ws, we, origin_seq   # buffer already satisfied

    # Extend window
    new_ws = max(0,     ws - need_left)
    new_we = min(plen,  we + need_right)

    # Clamp to ORIGIN_MAX_BP
    if (new_we - new_ws) > ORIGIN_MAX_BP:
        excess = (new_we - new_ws) - ORIGIN_MAX_BP
        new_ws += excess // 2
        new_we -= excess - excess // 2

    # Re-fetch from NCBI
    try:
        Entrez.email = ncbi_email
        if ncbi_api_key:
            Entrez.api_key = ncbi_api_key
        fs = max(1, new_ws + 1)
        fe = min(plen, new_we)
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta",
                               retmode="text", seq_start=fs, seq_stop=fe)
        rec = SeqIO.read(StringIO(handle.read()), "fasta")
        handle.close()
        time.sleep(NCBI_SLEEP_KEY if ncbi_api_key else NCBI_SLEEP)
        new_seq = str(rec.seq).upper()
        print(f"    {acc}: extended window {ws}-{we} → {new_ws}-{new_we} "
              f"(Rep buffer was {rep_start_in_win} / {len(origin_seq)-rep_end_in_win} bp)")
        return new_ws, new_we, new_seq
    except Exception as e:
        print(f"    {acc}: buffer-extension fetch failed ({e}), keeping original window")
        return ws, we, origin_seq


def run_prefiltered(args, out: Path):
    """
    Run pipeline stages 3-onwards on pre-filtered origins metadata.

    Stages executed:
      P1. Load pre-filtered metadata and apply repABC / NHR E. coli filters
      P2. Load or fetch origin sequences
      P3. Run ORF prediction + HMM to locate Rep gene within each origin
      P4. Enforce 200 bp Rep-gene buffer (extend window if needed)
      P5. Fetch GenBank annotations + trim non-origin CDS from window edges
      P6. Diversity selection (if more origins than target; otherwise keep all)
      P7. Functional OriV validation
      P8. Synthesis assessment
      P9. Write outputs
    """
    sep = "="*70
    Entrez.email = args.email
    if args.ncbi_api_key:
        Entrez.api_key = args.ncbi_api_key
    sleep = NCBI_SLEEP_KEY if args.ncbi_api_key else NCBI_SLEEP

    # P1 — Load and filter
    print(f"\n{sep}\nP1: Load pre-filtered metadata\n{sep}")
    df = load_plsdb(args.prefiltered)
    df = metadata_filter(df)
    n  = len(df)
    print(f"Entries after filters: {n}")

    # P2 — Load sequences
    print(f"\n{sep}\nP2: Load / fetch origin sequences\n{sep}")
    seqs: dict[str, str] = {}

    if args.seqs_fasta and Path(args.seqs_fasta).exists():
        print(f"  Loading sequences from {args.seqs_fasta}")
        for rec in SeqIO.parse(args.seqs_fasta, "fasta"):
            seqs[rec.id] = str(rec.seq).upper()
        print(f"  Loaded {len(seqs)} sequences")

    # Fetch any missing sequences from NCBI
    missing = [row for _, row in df.iterrows() if row["Accession"] not in seqs]
    if missing:
        print(f"  Fetching {len(missing)} sequences from NCBI ...")
        for i, row in enumerate(missing):
            acc  = row["Accession"]
            ws   = int(row["Origin_start_bp"])
            we   = int(row["Origin_end_bp"])
            plen = int(row["Length_bp"])
            loc  = str(row.get("Origin_GenBank_location", ""))
            try:
                if "join" in loc:
                    handle = Entrez.efetch(db="nucleotide", id=acc,
                                           rettype="fasta", retmode="text")
                    rec    = SeqIO.read(StringIO(handle.read()), "fasta")
                    handle.close()
                    full   = str(rec.seq).upper()
                    seqs[acc] = full[ws:] + full[:we]
                else:
                    fs = max(1, ws + 1); fe = min(plen, we)
                    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta",
                                           retmode="text", seq_start=fs, seq_stop=fe)
                    rec    = SeqIO.read(StringIO(handle.read()), "fasta")
                    handle.close()
                    seqs[acc] = str(rec.seq).upper()
            except Exception as e:
                print(f"    {acc}: fetch error — {e}")
            if (i + 1) % 10 == 0:
                print(f"    ... {i+1}/{len(missing)}")
            time.sleep(sleep)

    # P3 — ORF prediction + HMM to locate Rep within each origin
    print(f"\n{sep}\nP3: ORF prediction + Rep-protein HMM\n{sep}")
    seqs_present = {acc: s for acc, s in seqs.items() if s}
    all_orfs = predict_orfs_for_library(seqs_present)
    hmm_res  = run_rep_hmm_on_library(all_orfs, args.rip_hmm)

    # P4 — Enforce 200 bp Rep buffer
    print(f"\n{sep}\nP4: Enforcing ≥{REP_BUFFER_BP} bp Rep-gene buffer\n{sep}")
    updated_rows = []
    for _, row in df.iterrows():
        acc  = row["Accession"]
        ws   = int(row["Origin_start_bp"])
        we   = int(row["Origin_end_bp"])
        plen = int(row["Length_bp"])
        seq  = seqs.get(acc, "")
        if not seq:
            updated_rows.append(row.to_dict()); continue

        new_ws, new_we, new_seq = _enforce_rep_buffer(
            acc, ws, we, plen, seq, all_orfs, hmm_res, args.email, args.ncbi_api_key)

        seqs[acc] = new_seq
        rd = row.to_dict()
        rd.update({"Origin_start_bp": new_ws, "Origin_end_bp": new_we,
                   "Origin_span_bp": new_we - new_ws,
                   "Origin_GenBank_location": f"{new_ws}..{new_we}",
                   "Origin_sequence": new_seq})
        updated_rows.append(rd)

    df = pd.DataFrame(updated_rows)
    # Re-predict ORFs on any extended sequences
    seqs = {r["Accession"]: r.get("Origin_sequence","") or seqs.get(r["Accession"],"")
            for r in updated_rows}
    seqs = {k: v for k, v in seqs.items() if v}
    all_orfs = predict_orfs_for_library(seqs)

    # P5 — Fetch GenBank annotations + CDS trimming
    print(f"\n{sep}\nP5: Fetch GenBank annotations + trim non-origin CDS\n{sep}", flush=True)
    n_p5 = len(df)
    records: dict = {}
    trimmed_rows = []
    for i_p5, (_, row) in enumerate(df.iterrows()):
        acc  = row["Accession"]
        ws   = int(row["Origin_start_bp"])
        we   = int(row["Origin_end_bp"])
        plen = int(row["Length_bp"])
        rs   = int(row.get("Rep_orf_start", ws + REP_BUFFER_BP))
        re_  = int(row.get("Rep_orf_end",   we - REP_BUFFER_BP))
        seq  = seqs.get(acc, "")

        print(f"  [{i_p5+1}/{n_p5}] {acc}  window={ws}-{we}", flush=True)

        gb_features: list[dict] = []
        try:
            fs = max(1, ws - 500 + 1); fe = min(plen, we + 500)
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb",
                                   retmode="text", seq_start=fs, seq_stop=fe)
            gb_text = handle.read()
            handle.close()
            gbr = SeqIO.read(StringIO(gb_text), "genbank")
            records[acc] = gbr
            offset = fs - 1
            for feat in gbr.features:
                if feat.type != "CDS":
                    continue
                product = feat.qualifiers.get("product", ["unknown"])[0]
                gb_features.append({
                    "product":    product,
                    "protein_id": feat.qualifiers.get("protein_id", [""])[0],
                    "start":      int(feat.location.start) + offset,
                    "end":        int(feat.location.end)   + offset,
                    "strand":     feat.location.strand,
                    "category":   _cds_category(product),
                })
            print(f"    → {len(gb_features)} CDS features", flush=True)
        except Exception as e:
            import traceback
            print(f"    WARNING: GenBank fetch failed for {acc}: {e}", flush=True)
            traceback.print_exc()
        time.sleep(sleep)

        new_ws, new_we = trim_origin_window(ws, we, gb_features, rs, re_)
        if seq:
            rel_s = new_ws - ws; rel_e = new_we - ws
            new_seq = seq[max(0, rel_s): max(0, rel_e)]
        else:
            new_seq = seq
        seqs[acc] = new_seq

        rep_cds = [f["product"] for f in gb_features if f["category"] == "rep"]
        par_cds = [f["product"] for f in gb_features if f["category"] == "partition"]
        all_win = [f["product"] for f in gb_features if new_ws <= f["start"] < new_we]

        rd = row.to_dict()
        rd.update({"Origin_start_bp": new_ws, "Origin_end_bp": new_we,
                   "Origin_span_bp": new_we - new_ws,
                   "Origin_sequence": new_seq,
                   "Origin_elements": "; ".join(all_win[:6]),
                   "n_rep_elements": len(rep_cds),
                   "n_par_elements": len(par_cds),
                   "has_par": len(par_cds) > 0,
                   "_gb_features": gb_features})
        trimmed_rows.append(rd)

    df = pd.DataFrame(trimmed_rows)
    df = df[df["Origin_span_bp"].between(200, ORIGIN_MAX_BP)].reset_index(drop=True)
    seqs = {r["Accession"]: r.get("Origin_sequence","") or seqs.get(r["Accession"],"")
            for _, r in df.iterrows()}
    seqs = {k: v for k, v in seqs.items() if v}
    print(f"After trimming: {len(df)} origins")

    # P6 — Diversity selection (if more than target; else keep all)
    print(f"\n{sep}\nP6: Diversity selection\n{sep}")
    if len(df) > args.min_origins and "Origin_RIP_type" in df.columns:
        rip_stats = compute_rip_stats(df)
        lib = select_diversity_library(df, rip_stats, args.target_bp, args.min_origins)
    else:
        lib = df.copy()
        if "Selection_tier" not in lib.columns:
            lib["Selection_tier"] = "Pre-selected"
        print(f"  Keeping all {len(lib)} origins (below min-origins threshold)")

    # P7 — Validation
    print(f"\n{sep}\nP7: Functional OriV validation\n{sep}")
    lib_seqs = {acc: seqs[acc] for acc in lib["Accession"] if acc in seqs}
    lib_recs = {acc: records[acc] for acc in lib["Accession"] if acc in records}
    val_df   = validate_library(lib, lib_seqs, lib_recs, args.rip_hmm)
    all_orfs_lib = predict_orfs_for_library(lib_seqs)
    gb_map_lib   = parse_genbank_features(lib_recs)

    # P8 — Synthesis
    print(f"\n{sep}\nP8: Synthesis assessment\n{sep}")
    lib = assess_synthesis(lib, lib_seqs)

    # P9 — Outputs
    print(f"\n{sep}\nP9: Writing outputs\n{sep}")
    final = lib.merge(val_df, on="Accession", how="left")
    drop  = [c for c in final.columns if c.startswith("_") or c == "Plasmid_seq"]
    final = final.drop(columns=drop, errors="ignore")

    csv_path = out / "diversity_library_metadata.csv"
    final.to_csv(str(csv_path), index=False)
    print(f"CSV → {csv_path}  ({final.shape[0]} × {final.shape[1]})")

    fasta_path = out / "origin_sequences.fasta"
    write_fasta(lib_seqs, str(fasta_path))

    generate_pdf(lib, lib_seqs, lib_recs, all_orfs_lib, gb_map_lib, val_df,
                 str(out / "origin_plasmid_maps.pdf"))
    generate_overview(lib, val_df, str(out / "diversity_library_overview.png"))

    print(f"\n{sep}\nDONE\n{sep}")
    for cls in ["COMPLETE","FUNCTIONAL","MINIMAL","INCOMPLETE"]:
        print(f"  {cls:12s} {(val_df['oriv_class']==cls).sum():3d}/{len(val_df)}")


def main():
    args = parse_args()
    out  = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    sep = "="*70

    if args.prefiltered:
        run_prefiltered(args, out)
        return

    if args.plsdb_dir:
        run_plsdb_full(args, out)
        return

    # ── Legacy PLSDB-meta mode ───────────────────────────────────────────────
    Entrez.email = args.email

    print(f"\n{sep}\nSTAGE 1: Load PLSDB metadata and apply metadata filters\n{sep}")
    plsdb_df  = load_plsdb(args.plsdb_meta)
    filtered  = metadata_filter(plsdb_df)

    print(f"\n{sep}\nSTAGE 2: OriV finding (ORF prediction + Rep-protein HMM)\n{sep}")
    oriv_df   = find_oriv_in_plsdb(
        filtered, args.plsdb_fasta, args.rip_hmm,
        args.email, args.ncbi_api_key, args.batch_size, args.hmm_evalue)

    print(f"\n{sep}\nSTAGE 3: Origin window definition (≥{REP_BUFFER_BP} bp Rep buffer, ≤{ORIGIN_MAX_BP} bp)\n{sep}")
    windowed  = add_origin_windows(oriv_df)

    print(f"\n{sep}\nSTAGE 4: Fetch annotations + trim non-origin CDS\n{sep}")
    trimmed   = apply_trimming_from_annotations(windowed, args.email, args.ncbi_api_key)

    print(f"\n{sep}\nSTAGE 5: Deduplicate by RIP type\n{sep}")
    deduped   = deduplicate_and_build_metadata(trimmed)

    print(f"\n{sep}\nSTAGE 6: Tiered diversity selection\n{sep}")
    rip_stats = compute_rip_stats(deduped)
    lib       = select_diversity_library(deduped, rip_stats, args.target_bp, args.min_origins)

    print(f"\n{sep}\nSTAGE 7: Fetch GenBank annotations for {len(lib)} selected origins\n{sep}")
    seqs, records = fetch_annotations_for_library(lib, args.email, args.ncbi_api_key)
    lib["Origin_sequence"] = lib["Accession"].map(seqs)

    print(f"\n{sep}\nSTAGE 8: Functional OriV validation\n{sep}")
    val_df       = validate_library(lib, seqs, records, args.rip_hmm)
    all_orfs_lib = predict_orfs_for_library(seqs)
    gb_map_lib   = parse_genbank_features(records)

    print(f"\n{sep}\nSTAGE 9: Twist synthesis assessment\n{sep}")
    lib = assess_synthesis(lib, seqs)

    final = lib.merge(val_df, on="Accession", how="left")
    drop  = [c for c in final.columns if c.startswith("_") or c == "Plasmid_seq"]
    final = final.drop(columns=drop, errors="ignore")
    csv_path = out / "diversity_library_metadata.csv"
    final.to_csv(str(csv_path), index=False)
    print(f"\nCSV → {csv_path}  ({final.shape[0]} × {final.shape[1]})")
    fasta_path = out / "origin_sequences.fasta"
    write_fasta(seqs, str(fasta_path))

    print(f"\n{sep}\nSTAGE 10: Generating outputs\n{sep}")
    generate_pdf(lib, seqs, records, all_orfs_lib, gb_map_lib, val_df,
                 str(out / "origin_plasmid_maps.pdf"))
    generate_overview(lib, val_df, str(out / "diversity_library_overview.png"))

    print(f"\n{sep}\nDONE\n{sep}")
    print(f"  {csv_path}\n  {fasta_path}")
    print(f"  {out/'origin_plasmid_maps.pdf'}\n  {out/'diversity_library_overview.png'}")
    print("\nOriV classification:")
    for cls in ["COMPLETE","FUNCTIONAL","MINIMAL","INCOMPLETE"]:
        print(f"  {cls:12s} {(val_df['oriv_class']==cls).sum():3d}/{len(val_df)}")


if __name__ == "__main__":
    import traceback as _tb
    try:
        main()
    except Exception as _exc:
        print(f"\nFATAL ERROR: {_exc}", flush=True)
        _tb.print_exc()
        sys.exit(1)
