import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Load data
    expanded = pd.read_csv("results/expanded_golden_list_final_real_chromatin.csv")
    grna_summary = pd.read_csv("results/casoffinder/expanded_golden_list_grnas.csv")
    
    print("=" * 50)
    print("Expanded Golden List Analysis")
    print("=" * 50)
    print(f"Total variants: {len(expanded)}")
    print(f"\nGene distribution:")
    print(expanded["gene"].value_counts())
    
    print(f"\nVariant types:")
    print(expanded["variant_type"].value_counts())
    
    print(f"\nClinical significance:")
    print(expanded["clinical_significance"].value_counts())
    
    # Check TRCSS availability
    print(f"\nTRCSS availability:")
    print(f"  With trcss: {expanded['trcss'].notna().sum()}")
    print(f"  Average trcss: {expanded['trcss'].mean():.4f}")
    print(f"  Median trcss: {expanded['trcss'].median():.4f}")
    print(f"  Min trcss: {expanded['trcss'].min():.4f}")
    print(f"  Max trcss: {expanded['trcss'].max():.4f}")
    
    # Check chromatin scores
    chromatin_cols = [c for c in expanded.columns if "chromatin" in c.lower() and "score" in c.lower()]
    if chromatin_cols:
        print(f"\nChromatin score columns: {chromatin_cols}")
        for col in chromatin_cols:
            print(f"  {col}: {expanded[col].notna().sum()} non-null, mean {expanded[col].mean():.4f}")
    
    # Merge with gRNA summary
    merged = pd.merge(expanded, grna_summary, on=["hgvs_c"], how="left", suffixes=("", "_grna"))
    print(f"\nMerged with gRNA data: {len(merged)}")
    
    # Save merged data
    merged.to_csv("results/expanded_golden_list_with_grna.csv", index=False)
    print("Saved merged data to results/expanded_golden_list_with_grna.csv")
    
    # Plot TRCSS distribution
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    sns.histplot(expanded["trcss"].dropna(), kde=True, bins=20, color="royalblue")
    plt.title("TRCSS Distribution - Expanded Golden List")
    plt.xlabel("TRCSS")
    plt.ylabel("Count")
    plt.grid(alpha=0.3)
    
    # Plot TRCSS by gene
    plt.subplot(1, 2, 2)
    gene_trcss = expanded.groupby("gene")["trcss"].median().sort_values()
    sns.barplot(x=gene_trcss.values, y=gene_trcss.index, palette="viridis")
    plt.title("Median TRCSS by Gene")
    plt.xlabel("Median TRCSS")
    plt.tight_layout()
    plt.savefig("results/figures/expanded_trcss_analysis.png", dpi=300, bbox_inches="tight")
    print("\nSaved TRCSS analysis plot to results/figures/expanded_trcss_analysis.png")
    
    # If we have efficiency scores, plot that too
    if "efficiency_score" in merged.columns:
        plt.figure(figsize=(10, 6))
        valid = merged.dropna(subset=["trcss", "efficiency_score"])
        if len(valid) > 5:
            sns.scatterplot(x="trcss", y="efficiency_score", hue="gene", data=valid, s=100, alpha=0.8)
            plt.title("TRCSS vs gRNA Efficiency Score")
            plt.xlabel("TRCSS")
            plt.ylabel("gRNA Efficiency Score")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig("results/figures/expanded_trcss_vs_efficiency.png", dpi=300, bbox_inches="tight")
            print("Saved TRCSS vs efficiency plot to results/figures/expanded_trcss_vs_efficiency.png")
    
    # Create summary table
    summary = expanded.groupby("gene").agg(
        num_variants=("hgvs_c", "count"),
        avg_trcss=("trcss", "mean"),
        median_trcss=("trcss", "median"),
        min_trcss=("trcss", "min"),
        max_trcss=("trcss", "max")
    ).reset_index()
    summary.to_csv("results/expanded_golden_list_summary.csv", index=False)
    print("\nSaved gene-wise summary to results/expanded_golden_list_summary.csv")
    
    print("\n" + "=" * 50)
    print("Analysis complete!")
    print("=" * 50)

if __name__ == "__main__":
    main()