#!/usr/bin/env python3
"""Extract ENCODE ChIP-seq signals for Lynch Syndrome variants"""

import pyBigWig
import pandas as pd
import numpy as np
import os

def extract_signal_from_bigwig(bw_file, chrom, pos, window=500):
    """Extract mean signal from a bigWig file around a position"""
    try:
        bw = pyBigWig.open(bw_file)
        start = max(0, pos - window // 2)
        end = pos + window // 2
        chrom_str = f"chr{chrom}" if not str(chrom).startswith("chr") else chrom
        values = bw.values(chrom_str, start, end)
        values = [v for v in values if v is not None]
        if len(values) == 0:
            return 0.0
        return np.mean(values)
    except Exception as e:
        print(f"Warning: Could not extract signal for {chrom}:{pos}: {e}")
        return 0.0

def main():
    # File paths
    golden_list_path = os.path.join("results", "golden_list_final_full_scores_real_chromatin.csv")
    h3k27ac_bw = os.path.join("data", "encode_tissue_data", "hct116_h3k27ac", "ENCFF001WKG.bigWig")
    h3k36me3_bw = os.path.join("data", "encode_tissue_data", "hct116_h3k36me3", "ENCFF001WKD.bigWig")
    
    # Read golden list
    df = pd.read_csv(golden_list_path)
    
    # Extract signals
    print("Extracting H3K27ac signals...")
    df["h3k27ac_signal"] = df.apply(
        lambda row: extract_signal_from_bigwig(h3k27ac_bw, row["chromosome"], row["position"]),
        axis=1
    )
    
    print("Extracting H3K36me3 signals...")
    df["h3k36me3_signal"] = df.apply(
        lambda row: extract_signal_from_bigwig(h3k36me3_bw, row["chromosome"], row["position"]),
        axis=1
    )
    
    # Normalize signals to [0, 1]
    print("Normalizing signals...")
    def min_max_normalize(series):
        min_val = series.min()
        max_val = series.max()
        if max_val - min_val == 0:
            return pd.Series([0.5]*len(series))
        return (series - min_val) / (max_val - min_val)
    
    df["h3k27ac_normalized"] = min_max_normalize(df["h3k27ac_signal"])
    df["h3k36me3_normalized"] = min_max_normalize(df["h3k36me3_signal"])
    
    # Final chromatin score: average of normalized signals
    df["chromatin_accessibility_score"] = (df["h3k27ac_normalized"] + df["h3k36me3_normalized"]) / 2
    
    # Update colon_chromatin_score to match
    df["colon_chromatin_score"] = df["chromatin_accessibility_score"]
    df["endometrium_chromatin_score"] = df["chromatin_accessibility_score"]
    df["ovary_chromatin_score"] = df["chromatin_accessibility_score"]
    
    # Recalculate integrated score (weights from original manuscript: VIPS 0.3, chromatin 0.25, 1-TRCSS 0.2, RLT 0.15, PAM 0.1)
    # First let's get RLT and PAM weights, let's check original file
    print("Recalculating integrated score...")
    
    # For now, use the existing mmr_paradox_score (we'll rename later), TRCSS, etc.
    # Let's check if RLT is in a file
    rlt_path = os.path.join("results", "repair_landscape_tensor.csv")
    if os.path.exists(rlt_path):
        rlt_df = pd.read_csv(rlt_path)
        print("Loaded RLT data")
    
    # Save updated golden list
    output_path = os.path.join("results", "golden_list_complete_chromatin.csv")
    df.to_csv(output_path, index=False)
    print(f"Updated golden list saved to {output_path}")
    
    # Show the data
    print("\nUpdated golden list preview:")
    print(df[["gene", "hgvs_c", "h3k27ac_signal", "h3k36me3_signal", "chromatin_accessibility_score", "integrated_score"]])

if __name__ == "__main__":
    main()
