
import pandas as pd
from pathlib import Path
import gzip
import numpy as np

def load_gtf(gtf_path):
    print(f"Loading GTF from {gtf_path}...")
    # Load GTF, skip comment lines (start with #)
    gtf_df = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chromosome", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"
        ]
    )
    # Filter to only gene features
    gene_df = gtf_df[gtf_df["feature"] == "gene"].copy()
    
    # Extract gene_id and gene_name from attributes
    def parse_attribute(attr):
        attr_dict = {}
        for item in attr.split(";"):
            item = item.strip()
            if not item:
                continue
            key, value = item.split(" ", 1)
            attr_dict[key] = value.strip('"')
        return attr_dict
    
    gene_df["attributes"] = gene_df["attribute"].apply(parse_attribute)
    gene_df["gene_id"] = gene_df["attributes"].apply(lambda x: x.get("gene_id"))
    gene_df["gene_name"] = gene_df["attributes"].apply(lambda x: x.get("gene_name"))
    
    # Compute TSS: if strand is '+', TSS is start; if '-', TSS is end
    gene_df["tss"] = np.where(
        gene_df["strand"] == "+",
        gene_df["start"],
        gene_df["end"]
    )
    # Keep only necessary columns
    gene_df = gene_df[["chromosome", "strand", "tss", "gene_name", "start", "end"]].copy()
    # Normalize chromosome names
    gene_df["chromosome"] = "chr" + gene_df["chromosome"].astype(str)
    return gene_df

def compute_tss_distance(gene_df, chrom, pos, gene_name=None):
    # First, filter by chromosome
    chrom_genes = gene_df[gene_df["chromosome"] == chrom].copy()
    if len(chrom_genes) == 0:
        return np.nan
    
    if gene_name and pd.notna(gene_name):
        # If we have a gene name, try to match that first
        gene_match = chrom_genes[chrom_genes["gene_name"] == gene_name]
        if len(gene_match) > 0:
            # Compute distance to TSS of this gene
            gene = gene_match.iloc[0]
            return abs(pos - gene["tss"])
    
    # If no gene name match, find the nearest gene
    chrom_genes["distance"] = abs(chrom_genes["tss"] - pos)
    nearest_gene = chrom_genes.loc[chrom_genes["distance"].idxmin()]
    return nearest_gene["distance"]

if __name__ == "__main__":
    print("=== Computing Covariates for Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    # Load necessary files
    arrayed = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss.csv")
    gtf_path = base_dir / "data/Homo_sapiens.GRCh38.110.gtf.gz"
    
    # Load GTF gene data
    gene_df = load_gtf(gtf_path)
    
    # Compute TSS distance for each row
    print("\nComputing TSS distances...")
    tss_distances = []
    for _, row in arrayed.iterrows():
        chrom = row["chromosome"]
        pos = row["position"]
        gene_name = row.get("gene_name")
        tss_dist = compute_tss_distance(gene_df, chrom, pos, gene_name)
        tss_distances.append(tss_dist)
    arrayed["tss_distance"] = tss_distances
    
    # We already have transcription_direction! That's one covariate.
    # For transcription level, let's use a placeholder for now
    # (we could add ENCODE RNA-seq data later if available)
    arrayed["transcription_level"] = np.nan
    
    # Save the updated data
    output_path = base_dir / "results/arrayed_editing_full_trcss_covariates.csv"
    arrayed.to_csv(output_path, index=False)
    print(f"\nSaved updated data to {output_path}")
    
    # Show some stats
    print("\n=== Covariate Summary ===")
    print(f"Number of TSS distances computed: {arrayed['tss_distance'].notna().sum()}")
    print(f"TSS distance stats:\n{arrayed['tss_distance'].describe()}")
