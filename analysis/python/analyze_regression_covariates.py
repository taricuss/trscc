
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import GroupShuffleSplit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

if __name__ == "__main__":
    print("=== Regression Analysis with Covariates for Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss_covariates.csv")
    
    df_k562 = df.dropna(subset=[
        "K562_edited_percentage_endogenous", 
        "PRIDICT2_0_editing_Score_deep_K562", 
        "TRCSS",
        "tss_distance",
        "transcription_direction"
    ]).copy()
    
    # Log-transform TSS distance (since it's skewed)
    df_k562["log_tss_distance"] = np.log1p(df_k562["tss_distance"])
    
    # Define groups by genomic locus (chromosome + position) to prevent leakage
    df_k562['locus'] = df_k562['chromosome'] + '_' + df_k562['position'].astype(str)
    groups = df_k562['locus'].values
    
    # Define features
    features_pridict2 = ["PRIDICT2_0_editing_Score_deep_K562"]
    features_covariates = [
        "PRIDICT2_0_editing_Score_deep_K562", 
        "log_tss_distance", 
        "transcription_direction"
    ]
    features_covariates_trcss = [
        "PRIDICT2_0_editing_Score_deep_K562", 
        "log_tss_distance", 
        "transcription_direction",
        "TRCSS"
    ]
    y = df_k562["K562_edited_percentage_endogenous"]
    X_pridict2 = df_k562[features_pridict2]
    X_cov = df_k562[features_covariates]
    X_cov_trcss = df_k562[features_covariates_trcss]
    
    # Run multiple grouped train-test splits for robustness
    r2_pridict2_list = []
    r2_covariates_list = []
    r2_covariates_trcss_list = []
    gss = GroupShuffleSplit(n_splits=100, test_size=0.2, random_state=42)
    
    for train_idx, test_idx in gss.split(X_pridict2, y, groups=groups):
        X_train_p, X_test_p = X_pridict2.iloc[train_idx], X_pridict2.iloc[test_idx]
        X_train_cov, X_test_cov = X_cov.iloc[train_idx], X_cov.iloc[test_idx]
        X_train_cov_trcss, X_test_cov_trcss = X_cov_trcss.iloc[train_idx], X_cov_trcss.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        model_p = LinearRegression()
        model_p.fit(X_train_p, y_train)
        y_pred_p = model_p.predict(X_test_p)
        r2_pridict2_list.append(r2_score(y_test, y_pred_p))
        
        model_cov = LinearRegression()
        model_cov.fit(X_train_cov, y_train)
        y_pred_cov = model_cov.predict(X_test_cov)
        r2_covariates_list.append(r2_score(y_test, y_pred_cov))
        
        model_cov_trcss = LinearRegression()
        model_cov_trcss.fit(X_train_cov_trcss, y_train)
        y_pred_cov_trcss = model_cov_trcss.predict(X_test_cov_trcss)
        r2_covariates_trcss_list.append(r2_score(y_test, y_pred_cov_trcss))
    
    print("\n=== K562 Regression Results (Grouped by Genomic Locus) ===")
    print(f"Model 1 (PRIDICT2 only): R² = {np.mean(r2_pridict2_list):.3f} ± {np.std(r2_pridict2_list):.3f}")
    print(f"Model 2 (PRIDICT2 + covariates): R² = {np.mean(r2_covariates_list):.3f} ± {np.std(r2_covariates_list):.3f}")
    print(f"Model 3 (PRIDICT2 + covariates + TRCSS): R² = {np.mean(r2_covariates_trcss_list):.3f} ± {np.std(r2_covariates_trcss_list):.3f}")
    print(f"Average improvement (covariates vs PRIDICT2): {np.mean(r2_covariates_list) - np.mean(r2_pridict2_list):.3f}")
    print(f"Average improvement (TRCSS added): {np.mean(r2_covariates_trcss_list) - np.mean(r2_covariates_list):.3f}")
    
    # Also do a single grouped train-test split for example
    gss_single = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    train_idx, test_idx = next(gss_single.split(X_pridict2, y, groups=groups))
    
    X_train_p, X_test_p = X_pridict2.iloc[train_idx], X_pridict2.iloc[test_idx]
    X_train_cov, X_test_cov = X_cov.iloc[train_idx], X_cov.iloc[test_idx]
    X_train_cov_trcss, X_test_cov_trcss = X_cov_trcss.iloc[train_idx], X_cov_trcss.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    
    model_p = LinearRegression()
    model_p.fit(X_train_p, y_train)
    y_pred_p = model_p.predict(X_test_p)
    print(f"\nSingle grouped split - PRIDICT2 only: R² = {r2_score(y_test, y_pred_p):.3f}")
    print(f"Coefficients (PRIDICT2 only): intercept {model_p.intercept_:.3f}, PRIDICT2 {model_p.coef_[0]:.3f}")
    
    model_cov = LinearRegression()
    model_cov.fit(X_train_cov, y_train)
    y_pred_cov = model_cov.predict(X_test_cov)
    print(f"\nSingle grouped split - PRIDICT2 + covariates: R² = {r2_score(y_test, y_pred_cov):.3f}")
    print("Coefficients:")
    for name, coef in zip(features_covariates, model_cov.coef_):
        print(f"  {name}: {coef:.3f}")
    
    model_cov_trcss = LinearRegression()
    model_cov_trcss.fit(X_train_cov_trcss, y_train)
    y_pred_cov_trcss = model_cov_trcss.predict(X_test_cov_trcss)
    print(f"\nSingle grouped split - PRIDICT2 + covariates + TRCSS: R² = {r2_score(y_test, y_pred_cov_trcss):.3f}")
    print("Coefficients:")
    for name, coef in zip(features_covariates_trcss, model_cov_trcss.coef_):
        print(f"  {name}: {coef:.3f}")
    
    # Save results
    results = {
        "model": [
            "PRIDICT2 only", 
            "PRIDICT2 + covariates", 
            "PRIDICT2 + covariates + TRCSS"
        ],
        "mean_r2": [
            np.mean(r2_pridict2_list), 
            np.mean(r2_covariates_list), 
            np.mean(r2_covariates_trcss_list)
        ],
        "std_r2": [
            np.std(r2_pridict2_list), 
            np.std(r2_covariates_list), 
            np.std(r2_covariates_trcss_list)
        ]
    }
    results_df = pd.DataFrame(results)
    output_path = base_dir / "results/arrayed_editing_regression_covariates.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
