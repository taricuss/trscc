#!/usr/bin/env python3
"""Calculate correlations with confidence intervals using bootstrapping"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def bootstrap_correlation(x, y, n_bootstraps=10000, seed=42):
    """Calculate Spearman correlation and bootstrap confidence interval"""
    x = np.array(x)
    y = np.array(y)
    # Remove NaNs
    mask = ~np.isnan(x) & ~np.isnan(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    
    # Observed correlation
    obs_corr, obs_p = spearmanr(x, y)
    
    # Bootstrap with a fixed seed for reproducibility across reruns.
    rng = np.random.default_rng(seed)
    boot_corrs = []
    for _ in range(n_bootstraps):
        idx = rng.choice(n, n, replace=True)
        boot_corr, _ = spearmanr(x[idx], y[idx])
        if not np.isnan(boot_corr):
            boot_corrs.append(boot_corr)

    if not boot_corrs:
        return obs_corr, obs_p, np.nan, np.nan
    
    # 95% confidence interval
    ci_lower = np.percentile(boot_corrs, 2.5)
    ci_upper = np.percentile(boot_corrs, 97.5)
    
    return obs_corr, obs_p, ci_lower, ci_upper

def main():
    # Load final golden list
    final_df = pd.read_csv(os.path.join("results", "golden_list_final_submission_ready.csv"))
    
    # Remove rows with missing PRIDICT2 data
    df = final_df.dropna(subset=["PRIDICT2_0_editing_Score_deep_K562"]).copy()
    print(f"Analyzing {len(df)} variants with PRIDICT2 predictions")
    
    # Calculate correlations
    print("\n--- Correlation Analysis ---")
    
    # 1. TRCSS vs PRIDICT2 K562
    corr_trcss, p_trcss, ci_trcss_lower, ci_trcss_upper = bootstrap_correlation(
        df["trcss_final"], df["PRIDICT2_0_editing_Score_deep_K562"]
    )
    print(f"TRCSS vs PRIDICT2 K562: r = {corr_trcss:.3f}, 95% CI [{ci_trcss_lower:.3f}, {ci_trcss_upper:.3f}], p = {p_trcss:.3f}")
    
    # 2. PLANTI integrated score vs PRIDICT2 K562
    corr_plant, p_plant, ci_plant_lower, ci_plant_upper = bootstrap_correlation(
        df["complete_integrated_score"], df["PRIDICT2_0_editing_Score_deep_K562"]
    )
    print(f"PLANTI vs PRIDICT2 K562: r = {corr_plant:.3f}, 95% CI [{ci_plant_lower:.3f}, {ci_plant_upper:.3f}], p = {p_plant:.3f}")
    
    # 3. PLANTI without TRCSS vs PRIDICT2 K562
    df["no_trcss_score"] = (
        0.3 * df["vips_score"] + 
        0.25 * df["colon_chromatin_score"] + 
        0.15 * (1 - df["rlt"]) + 
        0.1 * (df["pam_confidence_tier"] == "Tier A").astype(float)
    )
    corr_no_trcss, p_no_trcss, ci_no_trcss_lower, ci_no_trcss_upper = bootstrap_correlation(
        df["no_trcss_score"], df["PRIDICT2_0_editing_Score_deep_K562"]
    )
    print(f"PLANTI (no TRCSS) vs PRIDICT2 K562: r = {corr_no_trcss:.3f}, 95% CI [{ci_no_trcss_lower:.3f}, {ci_no_trcss_upper:.3f}], p = {p_no_trcss:.3f}")
    
    # Save results
    results_df = pd.DataFrame({
        "Comparison": [
            "TRCSS vs PRIDICT2 K562",
            "PLANTI vs PRIDICT2 K562",
            "PLANTI (no TRCSS) vs PRIDICT2 K562"
        ],
        "Spearman_r": [corr_trcss, corr_plant, corr_no_trcss],
        "p_value": [p_trcss, p_plant, p_no_trcss],
        "CI_lower": [ci_trcss_lower, ci_plant_lower, ci_no_trcss_lower],
        "CI_upper": [ci_trcss_upper, ci_plant_upper, ci_no_trcss_upper]
    })
    results_path = os.path.join("results", "correlation_analysis_with_ci.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nResults saved to {results_path}")

if __name__ == "__main__":
    main()
