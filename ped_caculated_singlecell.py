import pandas as pd
import numpy as np
import itertools

# ==========================================
# 1. Configuration & Core Functions (PED V4.2)
# ==========================================

# Weight Configuration (V4.2 Optimized)
PED_WEIGHTS = {
    'w_dist': 0.2,  # Distribution consistency weight
    'w_cv': 0.5,    # Stability weight (Core)
    'w_snr': 0.3    # Accuracy weight
}


def calculate_ped_v4_2(df1, df2, meta_df, b1_cols, b2_cols):
    """
    Calculate PED V4.2 score and component metrics between two batches.
    """
    # --- A. Valid Protein Filtering (Pre-check) ---
    # Single-cell data is sparse; require values in at least 3 cells
    n1 = df1.notna().sum(axis=1)
    n2 = df2.notna().sum(axis=1)
    valid_mask = (n1 >= 3) & (n2 >= 3)

    # If too few valid proteins, cannot calculate
    if valid_mask.sum() < 10:
        return None

    d1 = df1.loc[valid_mask]
    d2 = df2.loc[valid_mask]

    # --- B. Dist Score (Distribution Consistency) ---
    mu1, mu2 = d1.mean(axis=1), d2.mean(axis=1)
    std1, std2 = d1.std(axis=1), d2.std(axis=1)
    epsilon = 1e-6

    s_loc = (np.abs(mu1 - mu2)) / (np.abs(mu1 - mu2) + std1 + std2 + epsilon)
    s_scale = (np.abs(std1 - std2)) / (std1 + std2 + epsilon)
    score_dist = 0.3 * s_loc.mean() + 0.7 * s_scale.mean()

    # --- C. CV Penalty (Stability) ---
    # Calculate variation within the same cell type (sc_m0)
    target_cell = 'sc_m0'
    cv_scores = []

    for cols in [b1_cols, b2_cols]:
        # Filter cells belonging to sc_m0 in this batch
        target_cols = [c for c in cols if meta_df.loc[c, 'celltype'] == target_cell]
        if len(target_cols) < 3: continue

        # In Log space, CV is approximated by SD
        sub_data = df1[target_cols] if cols[0] in df1.columns else df2[target_cols]
        sub_data = sub_data.loc[valid_mask]

        # Calculate Median SD for this batch
        sd_vals = sub_data.std(axis=1)
        cv_scores.append(sd_vals.median())

    # Take average of two batches, set circuit breaker cap at 1.0
    if not cv_scores:
        penalty_cv = 1.0
    else:
        penalty_cv = min(np.mean(cv_scores), 1.0)

    # --- D. SNR Penalty (Accuracy) ---
    # Measure discrimination between sc_m0 vs sc_u
    combined = pd.concat([d1, d2], axis=1)
    all_cols = b1_cols + b2_cols

    grp_m0 = [c for c in all_cols if meta_df.loc[c, 'celltype'] == 'sc_m0']
    grp_u = [c for c in all_cols if meta_df.loc[c, 'celltype'] == 'sc_u']

    if len(grp_m0) < 3 or len(grp_u) < 3:
        penalty_snr = 1.0
        real_snr = 0.0
    else:
        # Signal: Difference in means between groups
        signal = np.abs(combined[grp_m0].mean(axis=1) - combined[grp_u].mean(axis=1))
        # Noise: Mean of within-group SDs
        noise = (combined[grp_m0].std(axis=1) + combined[grp_u].std(axis=1)) / 2

        # Real Signal-to-Noise Ratio
        real_snr = (signal / (noise + epsilon)).median()
        # Penalty conversion
        penalty_snr = 1 / (1 + 2 * real_snr)

    # --- E. Total Score Calculation ---
    total = (PED_WEIGHTS['w_dist'] * score_dist +
             PED_WEIGHTS['w_cv'] * penalty_cv +
             PED_WEIGHTS['w_snr'] * penalty_snr)

    return {
        'Total': total,
        'Dist': score_dist,
        'CV': penalty_cv,
        'SNR_P': penalty_snr,
        'Real_SNR': real_snr,
        'Valid_N': valid_mask.sum()
    }


def get_detailed_metrics(protein_df, meta_df):
    """
    Iterate through all batch pairs, calculate weighted average metrics.
    """
    batches = meta_df['batch_chromatography'].unique()
    pairs = list(itertools.combinations(batches, 2))
    results = []

    for b1, b2 in pairs:
        # Find cell IDs belonging to these two batches
        c1 = meta_df[meta_df['batch_chromatography'] == b1].index
        c2 = meta_df[meta_df['batch_chromatography'] == b2].index

        # Take intersection (ensure these columns exist in data)
        cols1 = [c for c in c1 if c in protein_df.columns]
        cols2 = [c for c in c2 if c in protein_df.columns]

        # Sample size check
        if len(cols1) < 5 or len(cols2) < 5: continue

        # Calculate
        res = calculate_ped_v4_2(protein_df[cols1], protein_df[cols2], meta_df, cols1, cols2)
        if res:
            results.append(res)

    if not results: return None

    # Aggregate and calculate weighted averages
    df_res = pd.DataFrame(results)
    w = df_res['Valid_N']  # Weight by number of valid proteins

    return {
        'Total': np.average(df_res['Total'], weights=w),
        'Dist': np.average(df_res['Dist'], weights=w),
        'CV': np.average(df_res['CV'], weights=w),
        'SNR_P': np.average(df_res['SNR_P'], weights=w),
        'Real_SNR': np.average(df_res['Real_SNR'], weights=w)
    }


# ==========================================
# 2. Data Loading & Processing
# ==========================================

print("Loading Data...")
# 1. Load metadata
cells = pd.read_csv('data/Cells.csv', index_col=0)
if 'batch_chromatography' not in cells.columns:
    cells = cells.T
cells.index.name = 'cell_id'

# 2. Prepare Raw Data (Aggregation)
print("Processing Raw Data...")
peptides = pd.read_csv('data/Peptides-raw.csv')
# Remove peptide string column, group by protein and take median
proteins_raw = peptides.drop(columns=['peptide'], errors='ignore').groupby('protein').median()

# 3. Prepare Processed Data
print("Processing Processed Data...")
proteins_proc = pd.read_csv('data/Proteins-processed.csv', index_col=0)


# ==========================================
# 3. Run Evaluation & Print Report (Modified: Added Save Function)
# ==========================================

print("\nRunning PED Evaluation...")
metrics_raw = get_detailed_metrics(proteins_raw, cells)
metrics_proc = get_detailed_metrics(proteins_proc, cells)

# --- Prepare a list to collect data for saving ---
table_data = []

# Print headers
print("\n" + "=" * 80)
print(" PED V4.2 Detailed Component Breakdown (Raw vs Processed)")
print("=" * 80)
headers = ["Metric", "Raw Data", "Processed Data", "Change", "Status"]
print(f"{headers[0]:<20} | {headers[1]:<15} | {headers[2]:<15} | {headers[3]:<10} | {headers[4]:<10}")
print("-" * 80)


def print_and_save_row(label, val_raw, val_proc, higher_better=False):
    """
    Print to screen and save to list.
    """
    diff = val_proc - val_raw

    # Determine status (Better/Worse)
    if higher_better:
        status = "BETTER" if diff > 0 else "WORSE"
    else:  # Lower is better (for Penalties and Total Score)
        status = "BETTER" if diff < 0 else "WORSE"

    # 1. Print to console (maintain original format)
    print(f"{label:<20} | {val_raw:.4f}          | {val_proc:.4f}          | {diff:+.4f}     | {status}")

    # 2. Add to list (for saving)
    table_data.append({
        "Metric": label,
        "Raw Data": val_raw,
        "Processed Data": val_proc,
        "Change": diff,
        "Status": status
    })


# --- Execute output and collect data ---

# Total Score
print_and_save_row("PED Total Score", metrics_raw['Total'], metrics_proc['Total'])
print("-" * 80)

# Component Metrics
print_and_save_row("1. Dist Score", metrics_raw['Dist'], metrics_proc['Dist'])
print_and_save_row("2. CV Penalty", metrics_raw['CV'], metrics_proc['CV'])
print_and_save_row("3. SNR Penalty", metrics_raw['SNR_P'], metrics_proc['SNR_P'])
print("-" * 80)

# Real SNR (Higher is better)
print_and_save_row("Real SNR (Signal)", metrics_raw['Real_SNR'], metrics_proc['Real_SNR'], higher_better=True)
print("=" * 80)

# --- Save to file ---
output_filename = "PED_Validation_Results.csv"
df_results = pd.DataFrame(table_data)
df_results.to_csv(output_filename, index=False)

print(f"\n[Success] Result table successfully saved as: {output_filename}")
print(f"You can open this file directly in Excel to view detailed data.")