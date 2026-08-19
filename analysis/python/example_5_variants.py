#!/usr/bin/env python3
"""
Example script for PLANTI: Prime Editing for Lynch Syndrome Target Identification
Runs on 5 sample variants from the Golden List
"""

import os
import pandas as pd

def main():
    print("=== PLANTI Example: 5 Variants ===")
    
    # Load golden list
    golden_list_path = os.path.join("results", "golden_list_final_submission_ready.csv")
    if not os.path.exists(golden_list_path):
        print(f"Error: {golden_list_path} not found!")
        return
    
    golden_list = pd.read_csv(golden_list_path)
    
    # Take first 5 variants
    example_variants = golden_list.head(5).copy()
    print(f"\nLoaded 5 example variants from Golden List (total {len(golden_list)} variants):")
    print(example_variants[["gene", "hgvs_c", "chromosome", "position", "complete_integrated_score"]].to_string(index=False))
    
    # Save example output
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    example_output_path = os.path.join(output_dir, "example_5_variants_output.csv")
    example_variants.to_csv(example_output_path, index=False)
    
    print(f"\nExample output saved to: {example_output_path}")
    
    print("\n=== Key PLANTI Score Components for Example Variants ===")
    key_cols = [
        "gene", "hgvs_c", 
        "vips_score", "trcss_final", 
        "colon_chromatin_score", "rlt", 
        "complete_integrated_score"
    ]
    print(example_variants[key_cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    
    print("\n=== Done! ===")
    print("\nFor full analysis, see the manuscript in docs/PLANTI_Q1_MANUSCRIPT.Rmd")

if __name__ == "__main__":
    main()
