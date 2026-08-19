
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

def main():
    base_dir = Path(__file__).parent
    
    # Load all necessary data
    print("Loading expanded golden list...")
    expanded_golden = pd.read_csv(base_dir / "results/expanded_golden_list_final.csv")
    print("Loading expanded RLT...")
    rlt_expanded = pd.read_csv(base_dir / "results/repair_landscape_tensor_expanded.csv")
    print("Loading expanded TRCSS...")
    trcss_expanded = pd.read_csv(base_dir / "results/trcss_real_data_expanded.csv")
    print("Loading expanded ePRIDICT...")
    epridict_expanded = pd.read_csv(base_dir / "results/expanded_golden_list_nearest_epridict_predictions.csv")
    
    # Merge everything together
    print("Merging all data...")
    merged = expanded_golden.merge(rlt_expanded, on=['gene', 'hgvs_c', 'chromosome', 'position'], how='left')
    merged = merged.merge(trcss_expanded, on=['gene', 'hgvs_c'], how='left')
    merged = merged.merge(epridict_expanded, on=['gene', 'hgvs_c'], how='left')
    
    # Define VIPS score (from manuscript: 1.0 for frameshift/nonsense, 0 otherwise)
    # From expanded_golden: vips_score was called mmr_paradox_score originally?
    # Wait let's check what's in expanded_golden
    print("\nColumns in expanded golden list:")
    print(merged.columns.tolist())
    
    # Let's define components properly - use the _x suffix since those are from the original expanded list
    merged['vips_score'] = merged['mmr_paradox_score_x']  # Based on previous scripts
    # PAM confidence: 1 - population_pam_disruption_risk
    merged['pam_confidence'] = 1 - merged['population_pam_disruption_risk_x'].fillna(0)
    # 1 - TRCSS for integrated score
    merged['one_minus_trcss'] = 1 - merged['trcss_final'].fillna(0.5)
    # 1 - RLT for integrated score
    merged['one_minus_rlt'] = 1 - merged['rlt'].fillna(0)
    # Colon chromatin score
    merged['colon_chromatin_score'] = merged['colon_chromatin_score_x']
    
    # Define score components for correlation
    components = [
        'vips_score',
        'pam_confidence',
        'colon_chromatin_score',
        'one_minus_trcss',
        'one_minus_rlt'
    ]
    
    # Keep only components that exist in the dataframe
    available = [c for c in components if c in merged.columns]
    print(f"\nUsing components: {available}")
    
    # Compute correlation matrix
    n = len(available)
    corr_matrix = np.zeros((n, n))
    p_matrix = np.zeros((n, n))
    
    for i, col1 in enumerate(available):
        for j, col2 in enumerate(available):
            if i == j:
                corr_matrix[i, j] = 1.0
                p_matrix[i, j] = 1.0
            else:
                # Drop NA values for pairwise comparison
                valid = merged[[col1, col2]].dropna()
                if len(valid) < 2:
                    corr_matrix[i, j] = np.nan
                    p_matrix[i, j] = np.nan
                else:
                    r, p = spearmanr(valid[col1], valid[col2])
                    corr_matrix[i, j] = r
                    p_matrix[i, j] = p
    
    # Create DataFrames for output
    corr_df = pd.DataFrame(corr_matrix, index=available, columns=available)
    p_df = pd.DataFrame(p_matrix, index=available, columns=available)
    
    print("\nScore Component Correlation Matrix (Spearman's ρ):")
    print(corr_df.to_string())
    print("\n\np-values:")
    print(p_df.to_string())
    
    # Save to files
    corr_df.to_csv(base_dir / "results/score_component_orthogonality_correlation_matrix.csv")
    p_df.to_csv(base_dir / "results/score_component_orthogonality_p_values.csv")
    
    print("\n\nSaved to results/score_component_orthogonality_correlation_matrix.csv and results/score_component_orthogonality_p_values.csv")
    
    # Also save the merged expanded list with all components
    merged.to_csv(base_dir / "results/expanded_golden_list_final_full_scores.csv", index=False)
    print("\nSaved full expanded list with all scores to results/expanded_golden_list_final_full_scores.csv")

if __name__ == "__main__":
    main()

