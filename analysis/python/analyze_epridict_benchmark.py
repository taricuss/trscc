
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import r2_score

# Load data
base_dir = Path(__file__).parent
data_path = base_dir / "results/arrayed_editing_full_trcss.csv"
df = pd.read_csv(data_path)

# Keep only rows with non-null K562 editing efficiency, PRIDICT2, ePRIDICT, and TRCSS
cols_needed = ["K562_edited_percentage_endogenous", "PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full", "TRCSS"]
df = df.dropna(subset=cols_needed)

# Define groups by genomic locus (chromosome + position) to prevent leakage
df['locus'] = df['chromosome'] + '_' + df['position'].astype(str)
groups = df['locus'].values

print(f"Number of samples after filtering: {len(df)}")

# Prepare feature sets
X_pridict2 = df[["PRIDICT2_0_editing_Score_deep_K562"]]
X_pridict2_epridict = df[["PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full"]]
X_pridict2_epridict_trcss = df[["PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full", "TRCSS"]]
y = df["K562_edited_percentage_endogenous"]

# Run 100 grouped train-test splits to get stable R2 estimates
np.random.seed(42)
n_splits = 100
results = []
gss = GroupShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=42)

split_idx = 0
for train_idx, test_idx in gss.split(X_pridict2_epridict_trcss, y, groups=groups):
    X_train, X_test = X_pridict2_epridict_trcss.iloc[train_idx], X_pridict2_epridict_trcss.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    
    # Model 1: PRIDICT2 only
    model1 = LinearRegression()
    model1.fit(X_train[["PRIDICT2_0_editing_Score_deep_K562"]], y_train)
    y_pred1 = model1.predict(X_test[["PRIDICT2_0_editing_Score_deep_K562"]])
    r2_1 = r2_score(y_test, y_pred1)
    
    # Model 2: PRIDICT2 + ePRIDICT
    model2 = LinearRegression()
    model2.fit(X_train[["PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full"]], y_train)
    y_pred2 = model2.predict(X_test[["PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full"]])
    r2_2 = r2_score(y_test, y_pred2)
    
    # Model 3: PRIDICT2 + ePRIDICT + TRCSS
    model3 = LinearRegression()
    model3.fit(X_train, y_train)
    y_pred3 = model3.predict(X_test)
    r2_3 = r2_score(y_test, y_pred3)
    
    results.append({
        "split": split_idx,
        "r2_pridict2": r2_1,
        "r2_pridict2_epridict": r2_2,
        "r2_pridict2_epridict_trcss": r2_3,
        "delta_trcss_over_baseline": r2_3 - r2_2
    })
    split_idx += 1

# Convert to dataframe and compute summary stats
results_df = pd.DataFrame(results)
summary_stats = results_df.describe().T

print("\n--- Summary Statistics (Grouped by Genomic Locus) ---")
print(summary_stats[["mean", "std", "min", "25%", "50%", "75%", "max"]])

# Test if the improvement from adding TRCSS is statistically significant
# Using paired t-test and Wilcoxon signed-rank test
from scipy import stats
t_stat, t_pval = stats.ttest_rel(results_df["r2_pridict2_epridict_trcss"], results_df["r2_pridict2_epridict"])
w_stat, w_pval = stats.wilcoxon(results_df["r2_pridict2_epridict_trcss"], results_df["r2_pridict2_epridict"])
print(f"\nPaired t-test (TRCSS improvement): t={t_stat:.3f}, p={t_pval:.4g}")
print(f"Wilcoxon signed-rank test (TRCSS improvement): W={w_stat:.3f}, p={w_pval:.4g}")

# Save results
results_df.to_csv(base_dir / "results/arrayed_editing_epridict_benchmark.csv", index=False)
summary_stats.to_csv(base_dir / "results/arrayed_editing_epridict_benchmark_summary.csv")

print("\n--- Results saved to results/ directory ---")

