#!/usr/bin/env python3
"""Create the final golden list with all real data, merged duplicates, etc."""

import pandas as pd
import os

def main():
    # Load all necessary files
    chromatin_df = pd.read_csv(os.path.join("results", "golden_list_final_full_scores_real_chromatin.csv"))
    pridict2_df = pd.read_csv(os.path.join("results", "golden_list_pridict2_merged.csv"))
    rlt_trcss_df = pd.read_csv(os.path.join("results", "golden_list_final_full_scores_with_rlt_real_trcss_real_chromatin.csv"))
    
    # Create a key to merge on
    chromatin_df["variant_key"] = chromatin_df.apply(
        lambda row: f"{row['gene']}_{row['hgvs_c']}", axis=1
    )
    pridict2_df["variant_key"] = pridict2_df.apply(
        lambda row: f"{row['gene']}_{row['hgvs_c']}", axis=1
    )
    rlt_trcss_df["variant_key"] = rlt_trcss_df.apply(
        lambda row: f"{row['gene']}_{row['hgvs_c']}", axis=1
    )
    
    # Merge the data
    merged_df = chromatin_df.merge(
        pridict2_df[["variant_key", "PRIDICT2_0_editing_Score_deep_K562", "PRIDICT2_0_editing_Score_deep_HEK", 
                     "K562_percentile_to_librarydiverse", "HEK_percentile_to_librarydiverse", "K562_rank", "HEK_rank"]],
        on="variant_key", how="left"
    ).merge(
        rlt_trcss_df[["variant_key", "local_snv_density", "local_indel_density", "rlt", "trcss_final", "okseq_fd", "tx_dir", "vips_score", "pam_confidence_tier"]],
        on="variant_key", how="left", suffixes=("", "_rlt")
    )
    
    # Use the correct columns from rlt_trcss_df
    merged_df["vips_score"] = merged_df["vips_score_rlt"].fillna(merged_df["vips_score"])
    merged_df["pam_confidence_tier"] = merged_df["pam_confidence_tier_rlt"].fillna(merged_df["pam_confidence_tier"])
    merged_df = merged_df.drop(columns=["vips_score_rlt", "pam_confidence_tier_rlt"])
    
    # Normalize ENCODE chromatin signals properly
    # Min-max normalize H3K27ac and H3K36me3, then average
    def min_max_normalize(series):
        min_val = series.min()
        max_val = series.max()
        if max_val - min_val == 0:
            return pd.Series([0.5]*len(series))
        return (series - min_val) / (max_val - min_val)
    
    merged_df["h3k27ac_normalized"] = min_max_normalize(merged_df["h3k27ac_signal"])
    merged_df["h3k36me3_normalized"] = min_max_normalize(merged_df["h3k36me3_signal"])
    merged_df["chromatin_accessibility_score"] = (merged_df["h3k27ac_normalized"] + merged_df["h3k36me3_normalized"]) / 2
    
    # Update colon/endometrium/ovary scores to use the normalized ENCODE score
    merged_df["colon_chromatin_score"] = merged_df["chromatin_accessibility_score"]
    merged_df["endometrium_chromatin_score"] = merged_df["chromatin_accessibility_score"]
    merged_df["ovary_chromatin_score"] = merged_df["chromatin_accessibility_score"]
    
    # Deduplicate only truly identical variant records.
    # Important: do not collapse distinct variants that share the same locus.
    print("Original number of variants:", len(merged_df))
    dedup_cols = [
        "gene",
        "hgvs_c",
        "chromosome",
        "position",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    before_dedup = len(merged_df)
    merged_df = merged_df.drop_duplicates(subset=dedup_cols, keep="first")
    print("Removed exact duplicate rows:", before_dedup - len(merged_df))

    # Report potential same-position, different-allele records for manual review.
    same_pos = merged_df.groupby(["gene", "chromosome", "position"]).size()
    same_pos = same_pos[same_pos > 1]
    if not same_pos.empty:
        print("Note: same-position records retained because alleles/HGVS differ:")
        print(same_pos)

    print("Final number of variants:", len(merged_df))
    
    # Recalculate complete integrated score with correct weights (original manuscript weights)
    merged_df["complete_integrated_score"] = (
        0.3 * merged_df["vips_score"] + 
        0.25 * merged_df["colon_chromatin_score"] + 
        0.2 * (1 - merged_df["trcss_final"]) + 
        0.15 * (1 - merged_df["rlt"]) + 
        0.1 * (merged_df["pam_confidence_tier"] == "Tier A").astype(float)
    )
    
    # Reorder columns to be logical
    column_order = [
        "gene", "hgvs_c", "chromosome", "position", "ref", "alt", 
        "ReferenceAlleleVCF", "AlternateAlleleVCF", 
        "residual_mmr_activity", "vips_score", 
        "is_founder_haplotype", "population_pam_disruption_risk", "pam_confidence_tier",
        "data_source", "h3k27ac_signal", "h3k36me3_signal", "h3k27ac_normalized", "h3k36me3_normalized",
        "colon_chromatin_score", "endometrium_chromatin_score", "ovary_chromatin_score",
        "local_snv_density", "local_indel_density", "rlt",
        "trcss_final", "okseq_fd", "tx_dir",
        "PRIDICT2_0_editing_Score_deep_K562", "PRIDICT2_0_editing_Score_deep_HEK",
        "K562_percentile_to_librarydiverse", "HEK_percentile_to_librarydiverse",
        "K562_rank", "HEK_rank",
        "complete_integrated_score"
    ]
    merged_df = merged_df[column_order]
    
    # Sort by complete integrated score descending
    merged_df = merged_df.sort_values(by="complete_integrated_score", ascending=False)
    
    # Save final golden list
    output_path = os.path.join("results", "golden_list_final_submission_ready.csv")
    merged_df.to_csv(output_path, index=False)
    print(f"Final golden list saved to {output_path}")
    
    # Show preview
    print("\nFinal golden list preview:")
    print(merged_df[["gene", "hgvs_c", "data_source", "h3k27ac_signal", "colon_chromatin_score", "trcss_final", "complete_integrated_score", "PRIDICT2_0_editing_Score_deep_K562"]])
    
if __name__ == "__main__":
    main()
