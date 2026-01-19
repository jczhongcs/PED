import pandas as pd
import numpy as np
from itertools import combinations
import os
import warnings

warnings.filterwarnings('ignore')


file_before = './data/data/expfiles/quartet/confounded/protein/maxlfq/expdata_log.csv'
file_after = './data/data/expfiles/quartet/confounded/protein/maxlfq/expdata_ratio.csv'
file_meta = './data/data/expfiles/quartet/confounded/meta.csv'
output_csv = ('./PED_V4.2_Robust_Report_ratio_confounded.csv')


def load_data(file_path):
    print(f"📂 Loading: {os.path.basename(file_path)} ...")
    try:
        df = pd.read_csv(file_path, index_col=0)
        df = df.apply(pd.to_numeric, errors='coerce')
        return df
    except Exception as e:
        print(f"❌ Loading Error: {e}")
        return None


def calculate_ped_v4_robust(df, meta, b1, b2):
    """
    PED V4.2: Added robustness processing to prevent outliers from causing score explosion.
    """
    # 1. Get sample columns for the corresponding batches
    sub_meta = meta[meta['batch'].isin([b1, b2])]
    cols1 = [c for c in meta[meta['batch'] == b1]['run_id'] if c in df.columns]
    cols2 = [c for c in meta[meta['batch'] == b2]['run_id'] if c in df.columns]

    if not cols1 or not cols2: return None

    # --- A. Distribution Score (Dist) ---
    def get_stats(row):
        v = row.dropna()
        if len(v) < 3: return pd.Series([np.nan, np.nan])
        return pd.Series([v.mean(), v.std()])

    stats1 = df[cols1].apply(get_stats, axis=1)
    stats2 = df[cols2].apply(get_stats, axis=1)

    valid = stats1[0].notna() & stats2[0].notna()
    if valid.sum() < 10: return None

    # Filter out extreme outlier rows (e.g., abnormal proteins with SD > 10)
    valid = valid & (stats1[1] < 10) & (stats2[1] < 10)

    mu1, sd1 = stats1.loc[valid, 0], stats1.loc[valid, 1]
    mu2, sd2 = stats2.loc[valid, 0], stats2.loc[valid, 1]
    epsilon = 1e-6

    s_loc = ((mu1 - mu2).abs() / ((mu1 - mu2).abs() + sd1 + sd2 + epsilon)).mean()
    s_scale = ((sd1 - sd2).abs() / (sd1 + sd2 + epsilon)).mean()

    # Distribution Score (0.3/0.7 weights)
    score_dist = 0.3 * s_loc + 0.7 * s_scale

    # --- B. Stability (CV Penalty) ---
    cv_values = []
    for s in sub_meta['sample'].unique():
        s_reps = [r for r in sub_meta[sub_meta['sample'] == s]['run_id'] if r in df.columns]
        if len(s_reps) > 1:
            g_data = df.loc[valid, s_reps]
            g_cv = g_data.std(axis=1)  # Use SD to approximate CV for log data
            cv_values.append(g_cv.median())

    if cv_values:
        raw_cv = np.mean(cv_values)
        # ⚠️ Clipping mechanism: If CV > 1.0 (extremely abnormal), force it to 1.0
        score_cv = min(raw_cv, 1.0)
    else:
        score_cv = 0.2

    # --- C. Accuracy (SNR Penalty) ---
    # Calculate Signal-to-Noise Ratio (SNR)
    noise = score_cv
    sample_means = pd.DataFrame()
    for s in sub_meta['sample'].unique():
        s_reps = [r for r in sub_meta[sub_meta['sample'] == s]['run_id'] if r in df.columns]
        if s_reps:
            sample_means[s] = df.loc[valid, s_reps].mean(axis=1)

    if sample_means.shape[1] > 1:
        signal = sample_means.std(axis=1).median()
    else:
        signal = 0.0

    snr = signal / (noise + epsilon)

    # ⚠️ Optimized penalty function: Make scores for good data lower, approaching 0
    # New formula: 1/(1 + 2*SNR)
    score_snr = 1.0 / (1.0 + 2.0 * snr)

    # --- D. Total Score Synthesis ---
    # Weights: 20% Distribution + 50% Stability + 30% Accuracy
    total_score = 0.2* score_dist + 0.5 * score_cv + 0.3 * score_snr

    return {
        'Total': total_score,
        'Dist': score_dist,
        'CV': score_cv,
        'SNR_Penalty': score_snr,
        'Real_SNR': snr
    }


# ==========================================
# 3. Main Program (Modified: Save all components for Before and After)
# ==========================================
if __name__ == "__main__":
    df_before = load_data(file_before)
    df_after = load_data(file_after)
    meta = pd.read_csv(file_meta)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    if df_before is not None and df_after is not None:
        batches = sorted(meta['batch'].unique())
        pairs = list(combinations(batches, 2))

        print(f"\n{'Batch Pair':<20} | {'Before':<8} | {'After':<8} | {'Imp%':<6} | {'Real_SNR(Before -> After)':<25}")
        print("-" * 100)

        results = []
        for b1, b2 in pairs:
            # Calculate
            res_b = calculate_ped_v4_robust(df_before, meta, b1, b2)
            res_a = calculate_ped_v4_robust(df_after, meta, b1, b2)

            if res_b and res_a:
                imp = (res_b['Total'] - res_a['Total']) / res_b['Total'] * 100

                print(
                    f"{b1} vs {b2:<12} | {res_b['Total']:.4f}   | {res_a['Total']:.4f}   | {imp:5.1f}% | {res_b['Real_SNR']:.2f} -> {res_a['Real_SNR']:.2f}")

                results.append({
                    'Batch_A': b1, 'Batch_B': b2,
                    # --- Core Metrics (Before) ---
                    'Total_Before': res_b['Total'],
                    'Dist_Before': res_b['Dist'],
                    'CV_Before': res_b['CV'],
                    'SNR_Penalty_Before': res_b['SNR_Penalty'],
                    'Real_SNR_Before': res_b['Real_SNR'],

                    # --- Core Metrics (After) ---
                    'Total_After': res_a['Total'],
                    'Dist_After': res_a['Dist'],
                    'CV_After': res_a['CV'],
                    'SNR_Penalty_After': res_a['SNR_Penalty'],
                    'Real_SNR_After': res_a['Real_SNR'],

                    # --- Improvement Percentage ---
                    'Improvement_Pct': imp
                })

        if results:
            df_res = pd.DataFrame(results)
            # Adjust column order for better intuition
            cols = ['Batch_A', 'Batch_B',
                    'Total_Before', 'Total_After', 'Improvement_Pct',
                    'Dist_Before', 'Dist_After',
                    'CV_Before', 'CV_After',
                    'Real_SNR_Before', 'Real_SNR_After']
            # Keep only existing columns (in case keys are modified in the future)
            cols = [c for c in cols if c in df_res.columns] + [c for c in df_res.columns if c not in cols]

            df_res[cols].to_csv(output_csv, index=False)
            print(f"\n✅ Detailed comparison report saved: {output_csv}")
            print("   (The file now contains all original metrics like Dist_Before, CV_Before, etc., for plotting comparison graphs)")