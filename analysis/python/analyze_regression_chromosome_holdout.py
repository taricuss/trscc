
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

if __name__ == "__main__":
    print("=== Chromosome-Holdout Regression Analysis for Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss_covariates.csv")
    
    df_k562 = df.dropna(subset=[
        "K562_edited_percentage_endogenous", 
        "PRIDICT2_0_editing_Score_deep_K562", 
        "TRCSS",
        "tss_distance",
        "transcription_direction"
    ]).copy()
    
    # Log-transform TSS distance
    df_k562["log_tss_distance"] = np.log1p(df_k562["tss_distance"])
    
    # Get unique chromosomes and count number of samples per chromosome
    chrom_counts = df_k562["chromosome"].value_counts()
    print(f"\nChromosome counts:\n{chrom_counts}")
    
    # Only keep chromosomes with >= 5 samples
    chromosomes = chrom_counts[chrom_counts >= 5].index.tolist()
    print(f"\nTesting on chromosomes with >= 5 samples: {chromosomes}")
    
    features_pridict2 = ["PRIDICT2_0_editing_Score_deep_K562"]
    features_covariates_trcss = [
        "PRIDICT2_0_editing_Score_deep_K562", 
        "log_tss_distance", 
        "transcription_direction",
        "TRCSS"
    ]
    y_col = "K562_edited_percentage_endogenous"
    
    r2_pridict2_list = []
    r2_covariates_trcss_list = []
    test_chroms_used = []
    
    for test_chrom in chromosomes:
        print(f"\n=== Testing on chromosome {test_chrom} ===")
        
        # Split data
        test_df = df_k562[df_k562["chromosome"] == test_chrom].copy()
        train_df = df_k562[df_k562["chromosome"] != test_chrom].copy()
        
        if len(test_df) == 0 or len(train_df) == 0:
            print(f"Skipping {test_chrom}: not enough data")
            continue
        
        print(f"Train size: {len(train_df)}, Test size: {len(test_df)}")
        
        # Model 1: PRIDICT2 only
        model_p = LinearRegression()
        model_p.fit(train_df[features_pridict2], train_df[y_col])
        y_pred_p = model_p.predict(test_df[features_pridict2])
        r2_p = r2_score(test_df[y_col], y_pred_p)
        r2_pridict2_list.append(r2_p)
        print(f"PRIDICT2 only: R² = {r2_p:.3f}")
        
        # Model 2: PRIDICT2 + covariates + TRCSS
        model_cov_trcss = LinearRegression()
        model_cov_trcss.fit(train_df[features_covariates_trcss], train_df[y_col])
        y_pred_cov_trcss = model_cov_trcss.predict(test_df[features_covariates_trcss])
        r2_cov_trcss = r2_score(test_df[y_col], y_pred_cov_trcss)
        r2_covariates_trcss_list.append(r2_cov_trcss)
        print(f"PRIDICT2 + covariates + TRCSS: R² = {r2_cov_trcss:.3f}")
        
        test_chroms_used.append(test_chrom)
    
    print("\n=== Overall Chromosome-Holdout Results (chromosomes with >=5 samples) ===")
    print(f"PRIDICT2 only: mean R² = {np.mean(r2_pridict2_list):.3f} ± {np.std(r2_pridict2_list):.3f}")
    print(f"PRIDICT2 + covariates + TRCSS: mean R² = {np.mean(r2_covariates_trcss_list):.3f} ± {np.std(r2_covariates_trcss_list):.3f}")
    print(f"Average improvement: {np.mean(r2_covariates_trcss_list) - np.mean(r2_pridict2_list):.3f}")
    
    # Save results
    results = pd.DataFrame({
        "chromosome": test_chroms_used,
        "r2_pridict2_only": r2_pridict2_list,
        "r2_covariates_trcss": r2_covariates_trcss_list
    })
    output_path = base_dir / "results/arrayed_editing_regression_chromosome_holdout.csv"
    results.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
