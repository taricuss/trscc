# TRCSS: Transcription-Replication Context Score

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)][pyproject]
[![Zenodo](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXX-blue)][zenodo]

A portable, single-feature engineered covariate for **prime editing efficiency prediction**.

TRCSS captures the joint geometry of **replication fork directionality**,
**transcription strand orientation**, and **TSS distance** — three
genomic-context axes known to shape local DNA accessibility for prime
editing but **not exposed as separate plug-in terms** by standard public
prediction pipelines (PRIDICT2, DeepPE, PE-Designer).

> 💡 **Paper:** [TRCSS: Transcription-Replication Context Score for Improved Prime Editing Efficiency Prediction][paper-link] (submitted to *Nucleic Acids Research — Methods*.

---

## Quick start

```bash
# Install the package
git clone https://github.com/example/trcss.git && cd trcss
pip install -e .

# Compute TRCSS on a CSV with precomputed FD (fork directionality) and TxDir (gene strand)
trcss-compute tests/data/sample_input.csv --out-col trcss -o scored.csv
```

```python
from trcss import compute_trcss

# Head-on geometry (high TRCSS → favorable editing context)
compute_trcss(fd=-0.78, txdir=1)   # -> 0.89
```

See the [usage vignette][vignette] for the full 3-workflow tutorial.

---

## What TRCSS is, and is not

| | |
|---|---|
| **✅ A **biologically grounded drop-in single-feature** that combines FD × TxDir geometry into one [0, 1] number | ❌ A statistically orthogonal signal once TSS distance + TxDir are already in your model |
| **✅ Validated on 146 prime editing sites from Mathis et al. (2025) under grouped locus-disjoint CV | ❌ Claimed to work for strict cross-chromosome generalization on the current 146-site dataset |
| **✅ Produces **Ridge ΔR² = +0.097 drop-in improvement over PRIDICT2 + ePRIDICT (p = 2.45e-10) | ❌ A standalone predictor — use it *with* existing predictors, not instead of them |
| **✅ Computable from **only** FD data + gene strand (no transcript-level TSS annotations required) | ❌ A clinical product, a Lynch-Syndrome diagnostic tool, or a biomarker panel |

**Formal definition:**

$$
\mathrm{TRCSS} = \frac{1 - FD \times TxDir}{2}
$$

where FD ∈ [−1, +1] (fork directionality from OK-seq or Repli-seq),
TxDir ∈ {−1, +1} (+1 for plus-strand genes, −1 for minus-strand genes).

- Head-on geometry (FD × TxDir = −1) → TRCSS = 1
- Co-directional geometry (FD × TxDir = +1) → TRCSS = 0

---

## Headline results (from the manuscript)

Reproducible from the CSVs in `results/`:

| Benchmark | Metric | Baseline median | +TRCSS median | Δ median | p-value (Wilcoxon) |
|---|---|---:|---:|---:|---:|
| PRIDICT2 drop-in | R² | −0.089 | 0.019 | **+0.094** | 2.1×10⁻⁸ |
| PRIDICT2 + ePRIDICT SOTA drop-in (Ridge-primary) | R² | 0.113 | 0.232 | **+0.097** | 2.5×10⁻¹⁰ |
| SOTA + TSS + TxDir (fully controlled head-to-head, Ridge) | R² | 0.266 | 0.222 | −0.004 | 0.713 (ns) |
| SOTA + TRCSS vs DeepPE / PE-Designer | Spearman ρ | — | **0.532** | — | <10⁻⁴ |

Interpretation: TRCSS is a genuine **convenience & portability** feature —
it fills a real gap when existing pipelines don't expose TSS-distance / TxDir as
separate terms. When those mediating covariates are already in your model, the
incremental value is appropriately near-zero, as the formal mediation analysis
(~72% indirect effect via TSS distance + TxDir) predicts.

---

## Repository layout

```
trcss-repo/
├── src/trcss/                      # Installable Python package
│   ├── core.py                  # TRCSS formula + validators
│   └── cli.py                   # trcss-compute CLI entry point
│
├── manuscript/
│   └── TRCSS_manuscript.Rmd       # Full NAR Methods manuscript source
│
├── analysis/
│   ├── python/                 # Reproducibility scripts for every
│   │                           # every analysis in the paper
│   └── r/                      # Diagnostic R scripts
│
├── data/
│   ├── validation/             # 146-site validation set
│   ├── processed/          # Processed ClinVar/processed inputs
│   ├── epridict/            # ePRIDICT helper scripts
│   └── liftOver/            # liftOver chain (hg38→hg19)
│
├── results/                   # All pre-computed CSVs consumed by the Rmd
│   ├── casoffinder/            # Cas-OFFinder off-target tables
│   ├── founder/          # Founder-mutation candidate tables
│   └── figures/
│
├── plant_application/          # PLANTI illustrative downstream
│   ├── data/
│   └── results/
│
├── tests/                    # Unit tests + test data
├── docs/vignette.md        # Step-by-step vignette
├── .github/workflows/ci.yml   # CI: pytest + Rmd rendering on each tag
├── pyproject.toml
├── environment.yml            # Full R + Python pinned env (reproducible
├── requirements.txt
└── LICENSE (MIT)
```

---

## Reproducing the manuscript

### 1. Environment setup

```bash
# Option A: conda (includes R + Python)
conda env create -f environment.yml
conda activate trcss-env

# Option B: pip-only (Python package only
pip install -e ".[analysis]"
```

### 2. Re-run the validation analyses

```bash
# Compute TRCSS on the 146-site arrayed validation panel
python analysis/python/compute_trcss_arrayed.py

# Ridge drop-in CV vs PRIDICT2 + ePRIDICT SOTA (primary headline table)
python analysis/python/analyze_epridict_benchmark_ridge.py

# Fully controlled head-to-head: SOTA + TSS + TxDir vs +TRCSS
python analysis/python/analyze_head_to_head_trcss_vs_covariates.py

# LOCO chromosome-holdout
python analysis/python/analyze_regression_chromosome_holdout.py

# Ablation: FD alone vs TxDir alone vs FD×TxDir vs folded TRCSS
python analysis/python/trcss_ablation_real.py

# Score-component correlation matrix
python analysis/python/compute_score_component_correlation_matrix.py
```

### 3. Knit the manuscript

```bash
Rscript -e 'rmarkdown::render("manuscript/TRCSS_manuscript.Rmd")'
```

Outputs: `manuscript/TRCSS_manuscript.html` and a `.docx` version.

---

## Data availability & external accessions

All external datasets are public and re-downloadable with the exact accessions below. Access
dates are all 2026-07-24 per the manuscript's Data Availability section.

| Dataset | Accession(s) | Source |
|---|---|---|
| RPE-1 OK-seq RFD | GEO **GSE114017** | Petryk et al. 2016 / Chen et al. 2019 |
| K562 Repli-seq S-phase fractions (x ENCFF001GSC/GSD/GSE/GSF/GSG/GSH) | ENCODE (hg19) | ENCODE Project |
| K562 DNase-seq | ENCFF513IJB + ENCFF063ZSM | ENCODE |
| HCT116 H3K27ac/H3K36me3 | ENCFF997CJQ / ENCFF415LBK | ENCODE |
| TCGA MSI-H somatic SNVs | GDC COAD / READ / UCEC | NCI GDC |
| ClinVar MMR Path/LikelyPath variants | ClinVar GRCh38 2024-06 release | NCBI |
| 146-site prime editing efficiencies | GEO **GSE270866** | Mathis et al. (2025), CC BY 4.0 |

All processed derivative 146-site CSV is redistributed here under the CC BY 4.0 terms of
the original publication.

## License

MIT — see [LICENSE](LICENSE). Free for any use, commercial or non-commercial.
Please cite the paper if you use TRCSS in your own work.

---

## Citing

```bibtex
@article{trcss2026,
  title   = {TRCSS: Transcription-Replication Context Score for
             Improved Prime Editing Efficiency Prediction},
  journal = {Nucleic Acids Research},
  year    = {2026},
  doi     = {PENDING},
  note    = {Submitted to NAR Methods}
}

@software{trcss-zenodo,
  author    = {TRCSS Contributors},
  title   = {TRCSS software v1.0.0-submission},
  year    = {2026},
  doi     = {10.5281/zenodo.XXXXX},
  url     = {https://doi.org/10.5281/zenodo.XXXXX}
}
```

[paper-link]: https://doi.org/PENDING
[zenodo]: https://doi.org/10.5281/zenodo.XXXXX
[pyproject]: pyproject.toml
[vignette]: docs/vignette.md
