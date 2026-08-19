
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy import stats
from scipy.stats import spearmanr

base_dir = Path(__file__).parent
data_path = base_dir / "results/arrayed_editing_full_trcss.csv"
df = pd.read_csv(data_path)

cols_needed = ["K562_edited_percentage_endogenous",
               "PRIDICT2_0_editing_Score_deep_K562",
               "ePRIDICT_prediction_full",
               "TRCSS"]
df = df.dropna(subset=cols_needed)

df['locus'] = df['chromosome'].astype(str) + '_' + df['position'].astype(str)
groups = df['locus'].values
print(f"Number of samples after filtering: {len(df)}")
print(f"Number of unique loci: {df['locus'].nunique()}")

X_cols_baseline = ["PRIDICT2_0_editing_Score_deep_K562", "ePRIDICT_prediction_full"]
X_cols_full = X_cols_baseline + ["TRCSS"]
X = df[X_cols_full].copy()
y = df["K562_edited_percentage_endogenous"].values

def ndcg_at_k(y_true, y_pred, k=20):
    n = len(y_true)
    if n == 0:
        return 0.0
    k_use = min(k, n)
    order_pred = np.argsort(y_pred)[::-1]
    dcg = np.sum((2.0 ** y_true[order_pred[:k_use]] - 1.0) /
                 np.log2(np.arange(2, k_use + 2)))
    order_true = np.argsort(y_true)[::-1]
    idcg = np.sum((2.0 ** y_true[order_true[:k_use]] - 1.0) /
                  np.log2(np.arange(2, k_use + 2)))
    if idcg == 0:
        return 0.0
    return dcg / idcg

def precision_at_k(y_true, y_pred, k=10):
    n = len(y_true)
    if n == 0:
        return 0.0
    k_use = min(k, n)
    threshold = np.percentile(y_true, 80)
    order_pred = np.argsort(y_pred)[::-1]
    top_k_true = y_true[order_pred[:k_use]]
    return np.mean(top_k_true >= threshold)

def run_model(X_train_in, X_test_in, y_train_in, model_type="ridge"):
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train_in)
    X_test_sc = scaler.transform(X_test_in)
    if model_type == "ridge":
        alphas = np.logspace(-3, 3, 20)
        model = RidgeCV(alphas=alphas, cv=None)
        model.fit(X_train_sc, y_train_in)
    else:
        model = LinearRegression()
        model.fit(X_train_sc, y_train_in)
    y_pred = model.predict(X_test_sc)
    return y_pred, model

np.random.seed(42)
n_splits = 100
results = []
gss = GroupShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=42)

for split_idx, (train_idx, test_idx) in enumerate(gss.split(X, y, groups=groups)):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    for model_type in ["ridge", "ols"]:
        y_pred_baseline, _ = run_model(
            X_train[X_cols_baseline], X_test[X_cols_baseline], y_train, model_type=model_type
        )
        y_pred_full, _ = run_model(
            X_train[X_cols_full], X_test[X_cols_full], y_train, model_type=model_type
        )

        r2_b = r2_score(y_test, y_pred_baseline)
        r2_f = r2_score(y_test, y_pred_full)
        mae_b = mean_absolute_error(y_test, y_pred_baseline)
        mae_f = mean_absolute_error(y_test, y_pred_full)
        rmse_b = np.sqrt(mean_squared_error(y_test, y_pred_baseline))
        rmse_f = np.sqrt(mean_squared_error(y_test, y_pred_full))
        spear_b, _ = spearmanr(y_test, y_pred_baseline)
        spear_f, _ = spearmanr(y_test, y_pred_full)
        ndcg_b = ndcg_at_k(y_test, y_pred_baseline, k=20)
        ndcg_f = ndcg_at_k(y_test, y_pred_full, k=20)
        prec_b = precision_at_k(y_test, y_pred_baseline, k=10)
        prec_f = precision_at_k(y_test, y_pred_full, k=10)

        results.append({
            "split": split_idx,
            "model_type": model_type,
            "r2_baseline": r2_b,
            "r2_full": r2_f,
            "r2_delta": r2_f - r2_b,
            "mae_baseline": mae_b,
            "mae_full": mae_f,
            "mae_delta": mae_f - mae_b,
            "rmse_baseline": rmse_b,
            "rmse_full": rmse_f,
            "rmse_delta": rmse_f - rmse_b,
            "spearman_baseline": spear_b,
            "spearman_full": spear_f,
            "spearman_delta": spear_f - spear_b,
            "ndcg20_baseline": ndcg_b,
            "ndcg20_full": ndcg_f,
            "ndcg20_delta": ndcg_f - ndcg_b,
            "prec10_baseline": prec_b,
            "prec10_full": prec_f,
            "prec10_delta": prec_f - prec_b,
        })

results_df = pd.DataFrame(results)
results_df.to_csv(base_dir / "results/arrayed_editing_epridict_benchmark_ridge_per_split.csv",
                  index=False)

for mt in ["ridge", "ols"]:
    sub = results_df[results_df["model_type"] == mt].copy()
    print(f"\n{'='*70}")
    print(f"  MODEL TYPE: {mt.upper()} — PRIDICT2 + ePRIDICT (baseline) vs +TRCSS (full)")
    print(f"  n_splits = {len(sub)} GroupShuffleSplit (test_size=0.2, grouped by locus)")
    print(f"{'='*70}")

    for metric in ["r2", "mae", "rmse", "spearman", "ndcg20", "prec10"]:
        b_col = f"{metric}_baseline"
        f_col = f"{metric}_full"
        d_col = f"{metric}_delta"
        b_mean = sub[b_col].mean()
        b_med = sub[b_col].median()
        f_mean = sub[f_col].mean()
        f_med = sub[f_col].median()
        d_mean = sub[d_col].mean()
        d_med = sub[d_col].median()
        _, t_p = stats.ttest_rel(sub[f_col], sub[b_col])
        _, w_p = stats.wilcoxon(sub[f_col], sub[b_col])
        n_improved = int((sub[d_col] > 0).sum())

        higher_better = metric in ["r2", "spearman", "ndcg20", "prec10"]
        if higher_better:
            sig = ""
            if d_mean > 0 and t_p < 0.05:
                sig = " ***" if t_p < 1e-8 else " **" if t_p < 1e-4 else " *"
        else:
            sig = ""
            if d_mean < 0 and t_p < 0.05:
                sig = " ***" if t_p < 1e-8 else " **" if t_p < 1e-4 else " *"

        print(f"\n  [{metric.upper()}]")
        print(f"    Baseline: mean={b_mean:+.4f}  median={b_med:+.4f}")
        print(f"    Full:     mean={f_mean:+.4f}  median={f_med:+.4f}")
        print(f"    Δ (F-B):  mean={d_mean:+.4f}  median={d_med:+.4f}  "
              f"t_p={t_p:.4g}  Wilcox_p={w_p:.4g}  splits improved: {n_improved}/100{sig}")

summary_rows = []
for mt in ["ridge", "ols"]:
    sub = results_df[results_df["model_type"] == mt]
    for metric in ["r2", "spearman", "ndcg20", "prec10"]:
        b_col = f"{metric}_baseline"
        f_col = f"{metric}_full"
        d_col = f"{metric}_delta"
        _, t_p = stats.ttest_rel(sub[f_col], sub[b_col])
        _, w_p = stats.wilcoxon(sub[f_col], sub[b_col])
        summary_rows.append({
            "estimator": mt,
            "metric": metric,
            "baseline_mean": sub[b_col].mean(),
            "baseline_median": sub[b_col].median(),
            "full_mean": sub[f_col].mean(),
            "full_median": sub[f_col].median(),
            "delta_mean": sub[d_col].mean(),
            "delta_median": sub[d_col].median(),
            "t_pvalue": t_p,
            "wilcoxon_pvalue": w_p,
            "splits_improved": int((sub[d_col] > 0).sum()),
        })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(base_dir / "results/arrayed_editing_epridict_benchmark_ridge_summary.csv",
                  index=False)

print("\n\n--- Summary table saved to results/arrayed_editing_epridict_benchmark_ridge_summary.csv ---")
print(f"--- Per-split results saved to results/arrayed_editing_epridict_benchmark_ridge_per_split.csv ---")
