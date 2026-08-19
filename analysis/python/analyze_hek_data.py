
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

if __name__ == "__main__":
    print("=== Analyzing HEK293T Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss_covariates.csv")
    
    df_hek = df.dropna(subset=[
        "HEK293T_edited_percentage_endogenous", 
        "TRCSS"
    ]).copy()
    
    print(f"\nNumber of samples with HEK data: {len(df_hek)}")
    
    # Compute Spearman correlation
    rho, pval = spearmanr(df_hek["TRCSS"], df_hek["HEK293T_edited_percentage_endogenous"])
    print(f"\nTRCSS vs HEK293T editing efficiency: ρ = {rho:.3f}, p = {pval:.3f}")
    
    # Check if there's any interaction with other covariates
    df_hek["log_tss_distance"] = np.log1p(df_hek["tss_distance"])
    
    print("\n=== Descriptive stats for HEK data ===")
    print(df_hek[["HEK293T_edited_percentage_endogenous", "TRCSS", "log_tss_distance"]].describe())
    
    # Save results
    results = {
        "cell_line": ["HEK293T"],
        "n_samples": [len(df_hek)],
        "rho": [rho],
        "p_value": [pval]
    }
    results_df = pd.DataFrame(results)
    output_path = base_dir / "results/arrayed_editing_hek_analysis.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
