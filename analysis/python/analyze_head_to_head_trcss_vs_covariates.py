
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler
from scipy import stats

ALPHAS = np.logspace(-3, 3, 20)

def ndcg_at_k(y_true, y_pred, k):
    n = len(y_true)
    k = min(k, n)
    if k == 0:
        return 0.0
    order_true = np.argsort(y_true)[::-1]
    order_pred = np.argsort(y_pred)[::-1]
    true_ranks = np.zeros(n, dtype=int)
    true_ranks[order_true] = np.arange(n)
    dcg = 0.0
    idcg = 0.0
    for i in range(k):
        idx = order_pred[i]
        rank_of_pred = true_ranks[idx] + 1
        dcg += (2 ** (n - rank_of_pred) - 1) / np.log2(i + 2)
    for i in range(k):
        idx = order_true[i]
        rank_of_true = true_ranks[idx] + 1
        idcg += (2 ** (n - rank_of_true) - 1) / np.log2(i + 2)
    if idcg == 0:
        return 0.0
    return dcg / idcg

def precision_at_k(y_true, y_pred, k):
    n = len(y_true)
    k = min(k, n)
    if k == 0:
        return 0.0
    top_k_true_idx = set(np.argsort(y_true)[::-1][:k])
    top_k_pred_idx = set(np.argsort(y_pred)[::-1][:k])
    return len(top_k_true_idx & top_k_pred_idx) / k

def run_model(X_train, X_test, y_train, model_type="ridge"):
    if model_type == "ridge":
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        m = RidgeCV(alphas=ALPHAS)
        m.fit(X_train_s, y_train)
        y_pred = m.predict(X_test_s)
        alpha_used = m.alpha_
    else:
        m = LinearRegression()
        m.fit(X_train, y_train)
        y_pred = m.predict(X_test)
        alpha_used = None
    return y_pred, alpha_used

if __name__ == "__main__":
    print("=== HEAD-TO-HEAD CV: PRIDICT2+ePRIDICT+TSS+TxDir (NO TRCSS) vs +TRCSS ===\n")

    base_dir = Path(__file__).parent
    df = pd.read_csv(base_dir / "results/arrayed_editing_full_trcss_covariates.csv")

    cols_req = [
        "K562_edited_percentage_endogenous",
        "PRIDICT2_0_editing_Score_deep_K562",
        "ePRIDICT_prediction_full",
        "TRCSS",
        "tss_distance",
        "transcription_direction",
        "chromosome",
        "position"
    ]
    df = df.dropna(subset=cols_req).copy()
    df["log_tss_distance"] = np.log1p(df["tss_distance"])
    df['locus'] = df['chromosome'] + '_' + df['position'].astype(str)
    groups = df['locus'].values

    FEATURES_BASELINE = [
        "PRIDICT2_0_editing_Score_deep_K562",
        "ePRIDICT_prediction_full",
        "log_tss_distance",
        "transcription_direction"
    ]
    FEATURES_FULL = FEATURES_BASELINE + ["TRCSS"]

    y = df["K562_edited_percentage_endogenous"].values
    X_base = df[FEATURES_BASELINE].values
    X_full = df[FEATURES_FULL].values

    print(f"N samples: {len(df)}  |  N unique loci: {df['locus'].nunique()}")
    print(f"Baseline features: {FEATURES_BASELINE}")
    print(f"Full features (+TRCSS): {FEATURES_FULL}\n")

    np.random.seed(42)
    n_splits = 100
    gss = GroupShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=42)

    records = []
    for split_idx, (train_idx, test_idx) in enumerate(gss.split(X_full, y, groups=groups)):
        X_tr_b, X_te_b = X_base[train_idx], X_base[test_idx]
        X_tr_f, X_te_f = X_full[train_idx], X_full[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]

        for mtype in ["ridge", "ols"]:
            yp_b, a_b = run_model(X_tr_b, X_te_b, y_tr, mtype)
            yp_f, a_f = run_model(X_tr_f, X_te_f, y_tr, mtype)
            r2_b = r2_score(y_te, yp_b)
            r2_f = r2_score(y_te, yp_f)
            mae_b = mean_absolute_error(y_te, yp_b)
            mae_f = mean_absolute_error(y_te, yp_f)
            rmse_b = np.sqrt(mean_squared_error(y_te, yp_b))
            rmse_f = np.sqrt(mean_squared_error(y_te, yp_f))
            rho_b, _ = stats.spearmanr(y_te, yp_b) if len(y_te) >= 3 else (np.nan, np.nan)
            rho_f, _ = stats.spearmanr(y_te, yp_f) if len(y_te) >= 3 else (np.nan, np.nan)
            ndcg20_b = ndcg_at_k(y_te, yp_b, 20)
            ndcg20_f = ndcg_at_k(y_te, yp_f, 20)
            p10_b = precision_at_k(y_te, yp_b, 10)
            p10_f = precision_at_k(y_te, yp_f, 10)
            n_test = len(y_te)
            n_train = len(y_tr)
            records.append({
                "split": split_idx,
                "model": mtype,
                "n_train": n_train,
                "n_test": n_test,
                "r2_baseline": r2_b,
                "r2_full": r2_f,
                "delta_r2": r2_f - r2_b,
                "mae_baseline": mae_b,
                "mae_full": mae_f,
                "delta_mae": mae_f - mae_b,
                "rmse_baseline": rmse_b,
                "rmse_full": rmse_f,
                "delta_rmse": rmse_f - rmse_b,
                "rho_baseline": rho_b,
                "rho_full": rho_f,
                "delta_rho": rho_f - rho_b,
                "ndcg20_baseline": ndcg20_b,
                "ndcg20_full": ndcg20_f,
                "delta_ndcg20": ndcg20_f - ndcg20_b,
                "p10_baseline": p10_b,
                "p10_full": p10_f,
                "delta_p10": p10_f - p10_b,
                "alpha_baseline": a_b,
                "alpha_full": a_f,
            })

    res = pd.DataFrame(records)
    res.to_csv(base_dir / "results/head_to_head_trcss_vs_covariates_per_split.csv", index=False)

    # Summaries
    for mtype in ["ridge", "ols"]:
        rr = res[res["model"] == mtype].copy()
        print("=" * 78)
        print(f"  MODEL: {mtype.upper()}  (n_splits={len(rr)})")
        print("=" * 78)
        metrics = [
            ("R²", "r2_baseline", "r2_full", "delta_r2"),
            ("MAE", "mae_baseline", "mae_full", "delta_mae"),
            ("RMSE", "rmse_baseline", "rmse_full", "delta_rmse"),
            ("Spearman ρ", "rho_baseline", "rho_full", "delta_rho"),
            ("NDCG@20", "ndcg20_baseline", "ndcg20_full", "delta_ndcg20"),
            ("Precision@10", "p10_baseline", "p10_full", "delta_p10"),
        ]
        print(f"{'Metric':<16}{'Baseline (no TRCSS)':<26}{'Full (+TRCSS)':<20}{'Δ (Full - Baseline)':<18}")
        print("-" * 78)
        rows_summ = []
        for mname, bcol, fcol, dcol in metrics:
            b_arr = rr[bcol].values
            f_arr = rr[fcol].values
            d_arr = rr[dcol].values
            b_m, b_med = np.nanmean(b_arr), np.nanmedian(b_arr)
            f_m, f_med = np.nanmean(f_arr), np.nanmedian(f_arr)
            d_m, d_med = np.nanmean(d_arr), np.nanmedian(d_arr)
            d_q1, d_q3 = np.nanpercentile(d_arr, [25, 75])
            print(f"{mname:<16}"
                  f"mean={b_m:>7.3f} med={b_med:>7.3f}  "
                  f"mean={f_m:>7.3f} med={f_med:>7.3f}  "
                  f"mean={d_m:>+7.3f} med={d_med:>+7.3f}")
            rows_summ.append({
                "metric": mname,
                "baseline_mean": b_m,
                "baseline_median": b_med,
                "full_mean": f_m,
                "full_median": f_med,
                "delta_mean": d_m,
                "delta_median": d_med,
                "delta_q1": d_q1,
                "delta_q3": d_q3,
                "delta_iqr": d_q3 - d_q1,
            })
        pd.DataFrame(rows_summ).to_csv(
            base_dir / f"results/head_to_head_summary_{mtype}.csv", index=False
        )

        # Significance tests on ΔR² and ΔNDCG and ΔPrecision@10
        def sigtest(label, col):
            arr = rr[col].dropna().values
            t, tp = stats.ttest_1samp(arr, 0.0)
            w, wp = stats.wilcoxon(arr, zero_method='wilcox')
            print(f"\n  {label}  Δ mean={np.mean(arr):+.4f}  median={np.median(arr):+.4f}")
            print(f"     one-sample t: t={t:+.3f}  p={tp:.4g}")
            print(f"     Wilcoxon:  W={w:.1f}  p={wp:.4g}")
            n_pos = np.sum(arr > 0)
            n_neg = np.sum(arr < 0)
            print(f"     Splits improved: {n_pos}/{len(arr)}  (worsened: {n_neg})")
        print()
        sigtest("R² improvement", "delta_r2")
        sigtest("NDCG@20 improvement", "delta_ndcg20")
        sigtest("Precision@10 improvement", "delta_p10")
        sigtest("Spearman ρ improvement", "delta_rho")

        if mtype == "ridge":
            a_b = rr["alpha_baseline"].dropna()
            a_f = rr["alpha_full"].dropna()
            print(f"\n  Ridge alpha baseline: median={np.median(a_b):.4g}  IQR=[{np.percentile(a_b,25):.4g},{np.percentile(a_b,75):.4g}]")
            print(f"  Ridge alpha full:     median={np.median(a_f):.4g}  IQR=[{np.percentile(a_f,25):.4g},{np.percentile(a_f,75):.4g}]")

    print("\n[DONE] Per-split and summary CSVs written to results/")
