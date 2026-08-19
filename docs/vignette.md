# TRCSS Usage Vignette

This vignette walks through three common workflows for using the TRCSS
Python package:

1. [Basic: compute TRCSS from a CSV of loci](#1-basic-compute-trcss-from-a-csv)
2. [Scoring: append TRCSS to a PRIDICT2 / ePRIDICT predictions table](#2-append-trcss-to-predictor-outputs)
3. [Advanced: compute TRCSS from BED loci + OK-seq bigWig + GTF annotations](#3-advanced-from-loci-beds-bigwigs-gtf)

---

## 1. Basic: compute TRCSS from a CSV

The input CSV needs at minimum a **fork directionality** column (`fd`, in
[-1, +1]) and a **transcription direction** column (`txdir` ∈ {-1, +1}).

Sample `tests/data/sample_input.csv`:

```
chromosome,position,gene,fd,txdir,hgvs_c
chr3,37055850,MLH1,-0.55,-1,c.1528_1529del
```

### Using the CLI

```bash
pip install -e .   # from repo root

# Default column names (fd, txdir), write to stdout:
trcss-compute tests/data/sample_input.csv

# With explicit columns, save to file:
trcss-compute tests/data/sample_input.csv \
  --fd fd --txdir txdir --out-col trcss \
  -o /tmp/output_with_trcss.csv
```

### Using the Python API

```python
import pandas as pd
from trcss import compute_trcss, compute_trcss_dataframe

# Element-wise
trcss_value = compute_trcss(fd=-0.55, txdir=-1)  # head-on -> ~0.775

# On a DataFrame
df = pd.read_csv("tests/data/sample_input.csv")
compute_trcss_dataframe(df, inplace=True)
print(df[["gene", "fd", "txdir", "trcss"]])
```

---

## 2. Append TRCSS to predictor outputs

A common use case: you have a table of pegRNA predictions from PRIDICT2 or
ePRIDICT with the genomic coordinates of each target. You want to append a
TRCSS column as a drop-in single-feature covariate before re-training a
regression model.

```python
import pandas as pd
from trcss import compute_trcss_dataframe

# 1. Load predictions (typically has chr/pos + precomputed fd + txdir)
preds = pd.read_csv("pridict2_output.csv")

# 2. (If not already present) compute FD and TxDir:
#    See Section 3 below for computing FD from bigWigs and TxDir from GTF.

# 3. Append TRCSS
compute_trcss_dataframe(preds, fd_col="okseq_fd", txdir_col="gene_strand",
                        out_col="trcss", inplace=True)

# 4. Now use the trcss column alongside existing predictors.
#    Ridge drop-in use case:
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupShuffleSplit
import numpy as np

X_base = preds[["PRIDICT2_K562", "ePRIDICT"]].values
X_full = preds[["PRIDICT2_K562", "ePRIDICT", "trcss"]].values
y = preds["observed_efficiency_pct"].values
groups = preds["locus_id"].values  # one group per chr:pos

gss = GroupShuffleSplit(n_splits=100, test_size=0.2, random_state=1)
delta_r2 = []
for train_idx, test_idx in gss.split(X_base, y, groups):
    ridge_b = RidgeCV(alphas=np.logspace(-3, 3, 20)).fit(X_base[train_idx], y[train_idx])
    ridge_f = RidgeCV(alphas=np.logspace(-3, 3, 20)).fit(X_full[train_idx], y[train_idx])
    delta_r2.append(
        ridge_f.score(X_full[test_idx], y[test_idx])
        - ridge_b.score(X_base[test_idx], y[test_idx])
    )

print(f"Median Delta R2 (+TRCSS) = {np.median(delta_r2):+.3f}")
```

---

## 3. Advanced: from loci BED, OK-seq bigWig, GTF

If you only have the target loci (not precomputed FD or TxDir), you can
compute them from source files before running TRCSS.

### Dependencies

```bash
pip install pyBigWig bioframe
```

### Compute TxDir from a GTF

```python
import bioframe as bf
import pandas as pd

def txdir_for_loci(loci_df, gtf_path):
    """Return TxDir (+1 or -1) per locus by nearest protein-coding gene."""
    gtf = bf.read_table(gtf_path, schema="gtf")
    genes = gtf[(gtf["feature"] == "gene") & (gtf["gene_type"] == "protein_coding")]
    genes["start"] = genes["start"].astype(int)
    genes["end"] = genes["end"].astype(int)
    nearest = bf.nearest(loci_df, genes, suffixes=("_locus", "_gene"))
    txdir = nearest["strand_gene"].map({"+": 1.0, "-": -1.0}).fillna(0.0)
    return txdir.values
```

### Extract FD from an OK-seq RFD bigWig

```python
import pyBigWig

def fd_from_bigwig(loci_df, bw_path, window_bp=1000):
    """Extract mean FD value within +/- window_bp of each locus."""
    bw = pyBigWig.open(bw_path)
    out = []
    for _, row in loci_df.iterrows():
        chrom = row["chrom"]
        pos = int(row["start"])
        lo = max(0, pos - window_bp)
        hi = pos + window_bp
        try:
            vals = bw.values(chrom, lo, hi)
            vals = [v for v in vals if v is not None]
            out.append(float(pd.Series(vals).mean()) if vals else float("nan"))
        except Exception:
            out.append(float("nan"))
    bw.close()
    return pd.Series(out).values
```

### Full pipeline

```python
from trcss import compute_trcss

loci = pd.read_csv("my_targets.bed", sep="\t", header=None,
                   names=["chrom", "start", "end", "name"])

loci["txdir"] = txdir_for_loci(loci, "Homo_sapiens.GRCh38.110.gtf.gz")
loci["fd"]    = fd_from_bigwig(loci, "GSE114017_RPE1_RFD_hg38.bw")
loci["trcss"] = compute_trcss(loci["fd"], loci["txdir"])

loci.to_csv("my_targets_with_trcss.csv", index=False)
```

---

## 4. Running tests

```bash
pip install -e ".[dev]"
pytest -q
```

Expected output: all tests pass (currently ~15 unit tests covering input
validation, folding, computation, and DataFrame integration).

---

## 5. Regenerating the manuscript from the repo

The R Markdown source for the NAR Methods submission is in
`manuscript/TRCSS_manuscript.Rmd`. It reads all results from the
`results/` directory in the repo root.

### With conda (recommended)

```bash
conda env create -f environment.yml
conda activate trcss-env
Rscript -e 'rmarkdown::render("manuscript/TRCSS_manuscript.Rmd")'
```

### With R + Python separately

Install the R packages:

```r
install.packages(c("tidyverse", "knitr", "kableExtra", "reshape2", "boot",
                   "rmarkdown"))
```

Then:

```bash
Rscript -e 'rmarkdown::render("manuscript/TRCSS_manuscript.Rmd",
                              output_format="html_document")'
```

To produce the Word submission version:

```bash
Rscript -e 'rmarkdown::render("manuscript/TRCSS_manuscript.Rmd",
                              output_format="word_document")'
```
