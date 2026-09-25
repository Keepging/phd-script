import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(r"D:\博士\Protein contour\Phospho")
OUT = ROOT / "anomaly_first_codex"
OUT.mkdir(exist_ok=True)


def read_csv(path, **kwargs):
    return pd.read_csv(ROOT / path, **kwargs)


def cliff_delta(a, b):
    a = pd.Series(a).dropna().astype(float).values
    b = pd.Series(b).dropna().astype(float).values
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    d = 2 * u / (len(a) * len(b)) - 1
    return float(d), float(p)


lines = []
def log(s=""):
    print(s)
    lines.append(str(s))


log("ANOMALY-FIRST SCAN (Codex)")
log("=" * 80)

# ---------- Layer B 7-group presence-pattern stats ----------
presence_stats = read_csv("martinez_network_check/outputs/presence_pattern_stats.csv")
presence_sites = read_csv("martinez_network_check/outputs/presence_pattern_sites.csv")
log("\n[Layer B presence-pattern stats]")
log(presence_stats.to_string(index=False))
log("\nLayer B group counts:")
log(presence_sites["LocalizationGroup"].value_counts().reindex(
    ["C-only", "M-only", "N-only", "C&M", "C&N", "M&N", "C&M&N"]
).to_string())

# identify feature extrema per group from medians
median_cols = [c for c in presence_stats.columns if c.startswith("med_")]
extrema_rows = []
for _, row in presence_stats.iterrows():
    meds = row[median_cols].rename(lambda x: x.replace("med_", ""))
    extrema_rows.append({
        "feature": row["feature"],
        "eps2": row["eps2"],
        "max_group": meds.astype(float).idxmax(),
        "max_median": meds.astype(float).max(),
        "min_group": meds.astype(float).idxmin(),
        "min_median": meds.astype(float).min(),
        "max_minus_min": meds.astype(float).max() - meds.astype(float).min(),
    })
extrema = pd.DataFrame(extrema_rows)
extrema.to_csv(OUT / "layerB_feature_extrema.csv", index=False)
log("\nLayer B median extrema by feature:")
log(extrema.to_string(index=False))

# ---------- Stage6 mover/non-mover biophysics ----------
stage6 = read_csv("martinez_network_check/stage6_feature_comparison.csv")
reloc_pool = read_csv("martinez_network_check/stage6_reloc_resid_pool.csv")
stage6_feats = read_csv("martinez_network_check/stage6_biophysics_features.csv")
log("\n[Stage6 mover/non-mover biophysics]")
log(stage6.to_string(index=False))
log("\nStage6 reloc direction counts:")
log(reloc_pool["dir_resid"].value_counts(dropna=False).to_string())
log("\nStage6 mover/non-mover group counts:")
log(reloc_pool["grp"].value_counts(dropna=False).to_string())

confound_rows = []
for col in ["length", "abundance", "n_sites"]:
    mov = stage6_feats.loc[stage6_feats["grp"] == "mover", col]
    non = stage6_feats.loc[stage6_feats["grp"] == "non-mover", col]
    d, p = cliff_delta(mov, non)
    confound_rows.append({
        "variable": col,
        "mover_median": mov.median(),
        "non_mover_median": non.median(),
        "cliff_delta_mover_minus_non": d,
        "p": p,
        "n_mover": mov.notna().sum(),
        "n_non_mover": non.notna().sum(),
    })
stage6_confound = pd.DataFrame(confound_rows)
stage6_confound.to_csv(OUT / "stage6_confounder_scan.csv", index=False)
log("\nStage6 confounder scan:")
log(stage6_confound.to_string(index=False))

# ---------- Layer A raw vs controlled deltas ----------
layerA = read_csv("layerA_full/layerA_full_stats_summary.csv")
layerA["eps_drop_abs"] = layerA["raw_epsilon_squared"] - layerA["controlled_epsilon_squared"]
layerA["eps_retention"] = layerA["controlled_epsilon_squared"] / layerA["raw_epsilon_squared"]
layerA.to_csv(OUT / "layerA_raw_controlled_delta.csv", index=False)
log("\n[Layer A raw vs controlled]")
log(layerA[[
    "feature", "raw_epsilon_squared", "controlled_epsilon_squared",
    "eps_drop_abs", "eps_retention", "raw_p_BH", "controlled_p_BH"
]].to_string(index=False))

layerA_counts = read_csv("layerA_full/group_counts.csv")
log("\nLayer A group counts:")
log(layerA_counts.to_string(index=False))

# ---------- Compare Layer A controlled with Layer B ----------
name_map = {
    "backbone_dynamics": "backbone_dynamics",
    "sidechain_dynamics": "sidechain_dynamics",
    "disorder_propensity": "disorder_propensity",
    "helix_propensity": "helix_propensity",
    "sheet_propensity": "sheet_propensity",
    "coil_propensity": "coil_propensity",
    "earlyFolding": "earlyFolding",
}
comp = presence_stats[["feature", "eps2"]].rename(columns={"eps2": "layerB_eps2"}).merge(
    layerA[["feature", "controlled_epsilon_squared", "raw_epsilon_squared"]],
    on="feature",
    how="outer",
)
comp["A_control_minus_B"] = comp["controlled_epsilon_squared"] - comp["layerB_eps2"]
comp["A_raw_minus_B"] = comp["raw_epsilon_squared"] - comp["layerB_eps2"]
comp.to_csv(OUT / "layerA_layerB_eps_compare.csv", index=False)
log("\n[Layer A vs Layer B eps2 comparison]")
log(comp.to_string(index=False))

# ---------- Kinase correlation matrices ----------
kc = read_csv("kinase_correlation/kinase_substrate_correlation.csv")
kc_ns = kc[(kc["Kinase"] != kc["Substrate"]) & kc["Overall_Pearson_r"].notna()].copy()
kc_ns["is_cdk"] = kc_ns["Kinase"].str.upper().str.startswith("CDK")
log("\n[Kinase correlation]")
log(f"non-self valid pairs: {len(kc_ns)}")
for thr in [0.7, 0.8, 0.9, 0.95]:
    sub = kc_ns[kc_ns["Overall_Pearson_r"] > thr]
    log(f"r>{thr}: {len(sub)} pairs; CDK {int(sub['is_cdk'].sum())}/{len(sub)}")
log("Overall_r quantiles:")
log(kc_ns["Overall_Pearson_r"].quantile([0, .1, .25, .5, .75, .9, .95, .99, 1]).to_string())

# multi-kinase substrates
mk = pd.read_csv(ROOT / "multi_kinase_analysis/task3_multi_kinase_substrates.csv")
log("\nMulti-kinase cluster summary:")
log(mk["N_competing_kinases"].describe().to_string())
log("Top multi-kinase substrates:")
log(mk.head(15)[["Substrate", "N_competing_kinases", "Top_r", "Spread", "Kinases"]].to_string(index=False))

# relationship between multi-kinase count and abundance variability
hub = pd.read_csv(ROOT / "multi_kinase_analysis/hub_spatial_variability.csv")
merged_hub = mk.merge(hub, on="Substrate", how="left", suffixes=("_task3", "_hub"))
rho, p = spearmanr(
    merged_hub["N_competing_kinases_task3"],
    merged_hub["cross_compartment_SD"],
    nan_policy="omit",
)
log(f"\nSpearman N_competing_kinases vs cross_compartment_SD: rho={rho:.3f}, p={p:.3g}, n={merged_hub['cross_compartment_SD'].notna().sum()}")
log("Hub max_compartment counts:")
log(merged_hub["max_compartment"].value_counts(dropna=False).to_string())
merged_hub.to_csv(OUT / "multi_kinase_hub_join.csv", index=False)

# ---------- Original matrices shape and missingness ----------
for path in [
    "kinase_correlation/protein_abundance_matrix.csv",
    "kinase_correlation/substrate_abundance_matrix.csv",
    "kinase_correlation/kinase_substrate_correlation.csv",
]:
    df = read_csv(path)
    log(f"\nMatrix {path}: shape={df.shape}")
    log("columns: " + ", ".join(df.columns[:12].astype(str)) + (" ..." if len(df.columns) > 12 else ""))
    miss = df.isna().mean().sort_values(ascending=False).head(8)
    log("top missingness:")
    log(miss.to_string())

with open(OUT / "anomaly_scan_summary.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
