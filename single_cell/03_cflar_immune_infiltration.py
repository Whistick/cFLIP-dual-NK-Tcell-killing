"""Tumour CFLAR (or MHC-I) level vs NK/T-cell infiltration, per cancer type.

For each sample, tumour cells are split into CFLAR-high / CFLAR-low groups
(threshold computed within each cancer type), and the infiltration fraction of a
target microenvironment cell type (default: NK) among non-tumour cells is compared
between groups (Mann-Whitney U).

`adata` is the pre-processed pan-cancer AnnData, already in the session, with
log-normalised expression in `adata.X`.

Output (written to `OUTDIR`):
  - plot_df_for_R.csv          per-sample table for the R plotting script
  - sample_level_summary.csv   full per-sample summary
  - mhc_thresholds_per_cancer.csv
  - per_cancer_stats.csv       per-cancer test statistics
  - box_<cancer>.png           quick-look box plots
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import matplotlib.pyplot as plt

# -------- configuration ------------------------------------------------------
SCORE_GENES = ["CFLAR"]              # gene set scored on tumour cells
SAMPLE_COL = "Sample"
CELLTYPE_COL = "major_celltypes_tier2"
CANCER_COL = "CancerType_final"
TUMOR_LABEL = "Tumor"
TARGET_CELL_LABELS = ["NK"]         # microenvironment cell type to quantify
MIN_TUMOR_CELLS = 100               # drop samples with fewer tumour cells
GROUP_METHOD = "median"             # 'median' | 'quantile' | 'custom'
CUSTOM_THRESHOLD = 1.0
MIN_SAMPLES_PER_GROUP = 3
OUTDIR = "../data/cflar_infiltration"
# -----------------------------------------------------------------------------


def get_matrix(adata, var_names):
    """Return a dense expression DataFrame for `var_names` (handles sparse / missing)."""
    present = [g for g in var_names if g in adata.var_names]
    missing = [g for g in var_names if g not in adata.var_names]
    if not present:
        raise ValueError(f"None of the requested genes are in adata.var_names. Missing: {missing}")
    sub = adata[:, present].X
    if sp.issparse(sub):
        sub = sub.toarray()
    return pd.DataFrame(sub, index=adata.obs_names, columns=present), present, missing


def compute_score(adata, genes=None, agg="mean"):
    """Per-cell score = mean (or sum) expression of `genes`. Assumes normalised X."""
    if genes is None:
        genes = SCORE_GENES
    expr_df, _present, missing = get_matrix(adata, genes)
    score = expr_df.mean(axis=1) if agg == "mean" else expr_df.sum(axis=1)
    if missing:
        print(f"Warning: genes missing from adata.var_names: {missing}")
    return score


def fraction_target_in_microenv(adata, target_labels,
                                sample_col=SAMPLE_COL, celltype_col=CELLTYPE_COL,
                                tumor_label=TUMOR_LABEL):
    """Fraction of target cells among non-tumour cells, per sample."""
    rows = []
    for s, grp in adata.obs.groupby(sample_col):
        non_tumor = grp[grp[celltype_col] != tumor_label]
        n_total = len(non_tumor)
        if n_total == 0:
            n_target, frac = 0, np.nan
        else:
            n_target = non_tumor[celltype_col].isin(target_labels).sum()
            frac = n_target / n_total
        rows.append({"Sample": s, "n_total_non_tumor": n_total,
                     "n_target": n_target, "fraction_target": frac})
    return pd.DataFrame(rows).set_index("Sample")


def analyze(adata, score_genes=SCORE_GENES, group_method=GROUP_METHOD,
            custom_threshold=CUSTOM_THRESHOLD, min_tumor_cells=MIN_TUMOR_CELLS,
            target_labels=TARGET_CELL_LABELS, sample_col=SAMPLE_COL,
            celltype_col=CELLTYPE_COL, cancer_col=CANCER_COL, tumor_label=TUMOR_LABEL,
            min_samples_per_group=MIN_SAMPLES_PER_GROUP, outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)

    print("Computing per-cell tumour score ...")
    adata.obs["tumor_score"] = compute_score(adata, genes=score_genes).reindex(adata.obs_names)

    obs = adata.obs
    tumor_mask = obs[celltype_col] == tumor_label
    tumor_df = obs[tumor_mask].groupby(sample_col).agg(
        tumor_cells=("tumor_score", "size"),
        tumor_mhc_mean=("tumor_score", "mean"),
    )
    valid_samples = tumor_df[tumor_df["tumor_cells"] >= min_tumor_cells].index.tolist()
    print(f"{len(valid_samples)} samples pass tumour-cell cutoff (>= {min_tumor_cells})")

    frac_df = fraction_target_in_microenv(adata, target_labels, sample_col,
                                          celltype_col, tumor_label)
    sample_info = obs[[sample_col, cancer_col]].drop_duplicates().set_index(sample_col)
    merged = (tumor_df.join(sample_info, how="inner")
                      .join(frac_df, how="left")
                      .loc[valid_samples]
                      .dropna(subset=["fraction_target"]))

    # Split high/low within each cancer type
    if group_method == "median":
        compute_thr = lambda s: s.median()
    elif group_method == "quantile":
        compute_thr = lambda s: s.quantile(0.25)
    elif group_method == "custom":
        compute_thr = lambda s: custom_threshold
    else:
        raise ValueError("Unsupported group_method")

    merged["MHC_group"] = np.nan
    thresholds = {}
    for ct, sub in merged.groupby(cancer_col):
        thr = compute_thr(sub["tumor_mhc_mean"])
        thresholds[ct] = thr
        merged.loc[sub.index, "MHC_group"] = np.where(
            sub["tumor_mhc_mean"] <= thr, "MHC_I_low", "MHC_I_high")

    pd.Series(thresholds, name="mhc_threshold").to_csv(
        os.path.join(outdir, "mhc_thresholds_per_cancer.csv"))
    merged.to_csv(os.path.join(outdir, "sample_level_summary.csv"))

    plot_df = (merged.reset_index()
               .rename(columns={sample_col: "Sample", cancer_col: "CancerType"})
               [["Sample", "CancerType", "tumor_cells", "tumor_mhc_mean",
                 "n_total_non_tumor", "n_target", "fraction_target", "MHC_group"]])
    plot_df["percentage"] = plot_df["fraction_target"] * 100
    plot_df.to_csv(os.path.join(outdir, "plot_df_for_R.csv"), index=False)
    print("Wrote plot_df_for_R.csv (input for the R plotting script)")

    # Per-cancer test + quick-look box plots
    results = []
    for ct in merged[cancer_col].unique():
        sub = merged[merged[cancer_col] == ct]
        grp_low = sub.loc[sub["MHC_group"] == "MHC_I_low", "fraction_target"].dropna()
        grp_high = sub.loc[sub["MHC_group"] == "MHC_I_high", "fraction_target"].dropna()
        n_low, n_high, pval, test_name = len(grp_low), len(grp_high), np.nan, ""
        if n_low >= min_samples_per_group and n_high >= min_samples_per_group:
            _stat, pval = stats.mannwhitneyu(grp_low, grp_high, alternative="two-sided")
            test_name = "Mann-Whitney U"
        else:
            print(f"Cancer {ct}: too few samples (low={n_low}, high={n_high}); skipped")
        results.append({"CancerType": ct, "n_low": n_low, "n_high": n_high,
                        "pval": pval, "test": test_name})

        if len(sub) > 0:
            fig, ax = plt.subplots(figsize=(4, 5))
            ax.boxplot([grp_low.values, grp_high.values],
                       labels=["MHC_I_low", "MHC_I_high"])
            ax.set_title(str(ct))
            ax.set_ylabel("fraction of target cells (non-tumour denominator)")
            if not np.isnan(pval):
                ymax = max(grp_low.max() if n_low else 0, grp_high.max() if n_high else 0)
                ax.text(1.5, ymax * 0.95,
                        f"p={pval:.3e}\n(n_low={n_low}, n_high={n_high})", ha="center")
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, f"box_{ct}.png"), dpi=150)
            plt.close(fig)

    pd.DataFrame(results).to_csv(os.path.join(outdir, "per_cancer_stats.csv"), index=False)
    print("Done. Results in:", outdir)
    return merged, pd.DataFrame(results)


if __name__ == "__main__":
    # NK infiltration vs tumour CFLAR level
    merged, stats_df = analyze(adata, score_genes=["CFLAR"],
                               group_method="median", min_tumor_cells=100,
                               target_labels=["NK"])

    # Example variants used in the paper:
    #   target_labels=["NK", "T"]                      # combined NK + T
    #   score_genes=["HLA-A","HLA-B","HLA-C","B2M",...]  # MHC-I gene set instead of CFLAR
