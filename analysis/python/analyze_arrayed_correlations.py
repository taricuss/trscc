
import pandas as pd
from scipy.stats import spearmanr, norm
import numpy as np
from pathlib import Path
from statsmodels.stats.multitest import multipletests

def spearman_ci(r, n, confidence=0.95):
    if abs(r) >= 1:
        return (r, r)
    # Fisher z-transformation
    z = 0.5 * np.log((1 + r) / (1 - r))
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf((1 + confidence) / 2)
    z_low = z - z_crit * se
    z_high = z + z_crit * se
    # Back to r
    r_low = (np.exp(2 * z_low) - 1) / (np.exp(2 * z_low) + 1)
    r_high = (np.exp(2 * z_high) - 1) / (np.exp(2 * z_high) + 1)
    return (r_low, r_high)

if __name__ == "__main__":
    print("=== Analyzing Arrayed Editing Correlations ===")
    
    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss.csv")
    
    # Filter to non-NA K562 data
    df_k562 = df.dropna(subset=['K562_edited_percentage_endogenous', 'PRIDICT2_0_editing_Score_deep_K562', 'TRCSS'])
    print(f"\nNumber of variants with complete K562 data: {len(df_k562)}")
    
    print("\n=== K562 Correlations ===")
    
    # TRCSS vs K562 efficiency
    r_trcss_k562, p_trcss_k562 = spearmanr(df_k562['TRCSS'], df_k562['K562_edited_percentage_endogenous'])
    ci_trcss_k562 = spearman_ci(r_trcss_k562, len(df_k562))
    print(f"TRCSS vs K562 editing efficiency: ρ = {r_trcss_k562:.3f} [95% CI: {ci_trcss_k562[0]:.3f}, {ci_trcss_k562[1]:.3f}], p = {p_trcss_k562:.3f}")
    
    # PRIDICT2 vs K562 efficiency
    r_pridict2_k562, p_pridict2_k562 = spearmanr(df_k562['PRIDICT2_0_editing_Score_deep_K562'], df_k562['K562_edited_percentage_endogenous'])
    ci_pridict2_k562 = spearman_ci(r_pridict2_k562, len(df_k562))
    print(f"PRIDICT2 vs K562 editing efficiency: ρ = {r_pridict2_k562:.3f} [95% CI: {ci_pridict2_k562[0]:.3f}, {ci_pridict2_k562[1]:.3f}], p = {p_pridict2_k562:.3f}")
    
    # TRCSS vs PRIDICT2
    r_trcss_pridict2, p_trcss_pridict2 = spearmanr(df_k562['TRCSS'], df_k562['PRIDICT2_0_editing_Score_deep_K562'])
    ci_trcss_pridict2 = spearman_ci(r_trcss_pridict2, len(df_k562))
    print(f"TRCSS vs PRIDICT2 score: ρ = {r_trcss_pridict2:.3f} [95% CI: {ci_trcss_pridict2[0]:.3f}, {ci_trcss_pridict2[1]:.3f}], p = {p_trcss_pridict2:.3f}")
    
    # Now HEK293T
    print("\n=== HEK293T Correlations ===")
    df_hek = df.dropna(subset=['HEK293T_edited_percentage_endogenous', 'PRIDICT2_0_editing_Score_deep_HEK', 'TRCSS'])
    print(f"Number of variants with complete HEK293T data: {len(df_hek)}")
    
    r_trcss_hek, p_trcss_hek = spearmanr(df_hek['TRCSS'], df_hek['HEK293T_edited_percentage_endogenous'])
    ci_trcss_hek = spearman_ci(r_trcss_hek, len(df_hek))
    print(f"TRCSS vs HEK editing efficiency: ρ = {r_trcss_hek:.3f} [95% CI: {ci_trcss_hek[0]:.3f}, {ci_trcss_hek[1]:.3f}], p = {p_trcss_hek:.3f}")
    
    r_pridict2_hek, p_pridict2_hek = spearmanr(df_hek['PRIDICT2_0_editing_Score_deep_HEK'], df_hek['HEK293T_edited_percentage_endogenous'])
    ci_pridict2_hek = spearman_ci(r_pridict2_hek, len(df_hek))
    print(f"PRIDICT2 vs HEK editing efficiency: ρ = {r_pridict2_hek:.3f} [95% CI: {ci_pridict2_hek[0]:.3f}, {ci_pridict2_hek[1]:.3f}], p = {p_pridict2_hek:.3f}")
    
    # Collect all p-values for multiple testing correction
    p_values = [p_trcss_k562, p_pridict2_k562, p_trcss_pridict2, p_trcss_hek, p_pridict2_hek]
    
    # Apply Benjamini-Hochberg FDR correction
    _, p_adj, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Save correlation results
    results = {
        'dataset': ['K562', 'K562', 'K562', 'HEK293T', 'HEK293T'],
        'x': ['TRCSS', 'PRIDICT2', 'TRCSS', 'TRCSS', 'PRIDICT2'],
        'y': ['editing_efficiency', 'editing_efficiency', 'PRIDICT2', 'editing_efficiency', 'editing_efficiency'],
        'n': [len(df_k562), len(df_k562), len(df_k562), len(df_hek), len(df_hek)],
        'rho': [r_trcss_k562, r_pridict2_k562, r_trcss_pridict2, r_trcss_hek, r_pridict2_hek],
        'ci_low': [ci_trcss_k562[0], ci_pridict2_k562[0], ci_trcss_pridict2[0], ci_trcss_hek[0], ci_pridict2_hek[0]],
        'ci_high': [ci_trcss_k562[1], ci_pridict2_k562[1], ci_trcss_pridict2[1], ci_trcss_hek[1], ci_pridict2_hek[1]],
        'p': p_values,
        'p_adj': p_adj  # FDR-corrected p-values
    }
    results_df = pd.DataFrame(results)
    output_path = base_dir / "results/arrayed_editing_correlations.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
    print("\n=== FDR-Corrected P-Values:")
    for i, (x, y, dataset, p, padj) in enumerate(zip(results['x'], results['y'], results['dataset'], p_values, p_adj)):
        print(f"{dataset} {x} vs {y}: p = {p:.3f}, p_adj = {padj:.3f}")
