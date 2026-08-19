
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import GroupShuffleSplit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats

if __name__ == "__main__":
    print("=== Regression Analysis for Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss.csv")
    
    df_k562 = df.dropna(subset=['K562_edited_percentage_endogenous', 'PRIDICT2_0_editing_Score_deep_K562', 'TRCSS'])
    
    # Define groups by genomic locus (chromosome + position) to prevent leakage
    df_k562['locus'] = df_k562['chromosome'] + '_' + df_k562['position'].astype(str)
    groups = df_k562['locus'].values
    
    # Features and target
    X_pridict2 = df_k562[['PRIDICT2_0_editing_Score_deep_K562']]
    X_combined = df_k562[['PRIDICT2_0_editing_Score_deep_K562', 'TRCSS']]
    y = df_k562['K562_edited_percentage_endogenous']
    
    # Run multiple grouped train-test splits for robustness
    r2_pridict2_list = []
    r2_combined_list = []
    gss = GroupShuffleSplit(n_splits=100, test_size=0.2, random_state=42)
    
    for train_idx, test_idx in gss.split(X_pridict2, y, groups=groups):
        X_train_p, X_test_p = X_pridict2.iloc[train_idx], X_pridict2.iloc[test_idx]
        X_train_c, X_test_c = X_combined.iloc[train_idx], X_combined.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        model_p = LinearRegression()
        model_p.fit(X_train_p, y_train)
        y_pred_p = model_p.predict(X_test_p)
        r2_p = r2_score(y_test, y_pred_p)
        r2_pridict2_list.append(r2_p)
        
        model_c = LinearRegression()
        model_c.fit(X_train_c, y_train)
        y_pred_c = model_c.predict(X_test_c)
        r2_c = r2_score(y_test, y_pred_c)
        r2_combined_list.append(r2_c)
    
    # Statistical tests
    delta_r2 = np.array(r2_combined_list) - np.array(r2_pridict2_list)
    t_stat, t_pval = stats.ttest_rel(r2_combined_list, r2_pridict2_list)
    w_stat, w_pval = stats.wilcoxon(r2_combined_list, r2_pridict2_list)
    
    print("\n=== K562 Regression Results (Grouped by Genomic Locus) ===")
    print(f"Model 1 (PRIDICT2 only): Mean R² = {np.mean(r2_pridict2_list):.3f} ± {np.std(r2_pridict2_list):.3f}, Median R² = {np.median(r2_pridict2_list):.3f}, IQR = {np.percentile(r2_pridict2_list, 75) - np.percentile(r2_pridict2_list, 25):.3f}")
    print(f"Model 2 (PRIDICT2 + TRCSS): Mean R² = {np.mean(r2_combined_list):.3f} ± {np.std(r2_combined_list):.3f}, Median R² = {np.median(r2_combined_list):.3f}, IQR = {np.percentile(r2_combined_list, 75) - np.percentile(r2_combined_list, 25):.3f}")
    print(f"Mean improvement in R²: {np.mean(delta_r2):.3f}, Median improvement: {np.median(delta_r2):.3f}")
    print(f"Paired t-test: t = {t_stat:.3f}, p = {t_pval:.3e}")
    print(f"Wilcoxon signed-rank test: W = {w_stat:.3f}, p = {w_pval:.3e}")
    
    # Also do a single grouped train-test split for example
    gss_single = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    train_idx, test_idx = next(gss_single.split(X_pridict2, y, groups=groups))
    
    X_train_p, X_test_p = X_pridict2.iloc[train_idx], X_pridict2.iloc[test_idx]
    X_train_c, X_test_c = X_combined.iloc[train_idx], X_combined.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    
    model_p = LinearRegression()
    model_p.fit(X_train_p, y_train)
    y_pred_p = model_p.predict(X_test_p)
    print(f"\nSingle grouped split - PRIDICT2 only: R² = {r2_score(y_test, y_pred_p):.3f}")
    print(f"Coefficients (PRIDICT2): intercept {model_p.intercept_:.3f}, PRIDICT2 {model_p.coef_[0]:.3f}")
    
    model_c = LinearRegression()
    model_c.fit(X_train_c, y_train)
    y_pred_c = model_c.predict(X_test_c)
    print(f"Single grouped split - PRIDICT2 + TRCSS: R² = {r2_score(y_test, y_pred_c):.3f}")
    print(f"Coefficients: intercept {model_c.intercept_:.3f}, PRIDICT2 {model_c.coef_[0]:.3f}, TRCSS {model_c.coef_[1]:.3f}")
    
    # Save results
    results = {
        'model': ['PRIDICT2 only', 'PRIDICT2 + TRCSS'],
        'mean_r2': [np.mean(r2_pridict2_list), np.mean(r2_combined_list)],
        'std_r2': [np.std(r2_pridict2_list), np.std(r2_combined_list)],
        'median_r2': [np.median(r2_pridict2_list), np.median(r2_combined_list)],
        'iqr_r2': [np.percentile(r2_pridict2_list, 75) - np.percentile(r2_pridict2_list, 25), np.percentile(r2_combined_list, 75) - np.percentile(r2_combined_list, 25)]
    }
    results_df = pd.DataFrame(results)
    output_path = base_dir / "results/arrayed_editing_regression.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
