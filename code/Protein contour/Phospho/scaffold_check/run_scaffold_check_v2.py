"""Step-2 scaffold finding: protein-level confounder control.

Read-only for project inputs. Outputs are written only under scaffold_check/.

This reconstructs the whole-protein table used by D:\博士\Phospho\0331.py because
the original notebook/script produced figures but did not save prot_df as a CSV.
The comparison unit is one protein, not residues or phosphosites.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, t


ROOT = Path.cwd()
OUT = ROOT / "scaffold_check"
BIOPHYS = ROOT / "biophys_json"
RNG = np.random.default_rng(20260609)

HIGH_SITES = ROOT / "high_mobility_sites.csv"
LOW_SITES = ROOT / "low_mobility_sites.csv"
GENE_TO_UNIPROT = ROOT / "gene_to_uniprot.json"

FIGURE_PATHS = [
    r"D:\博士\Phospho\Fig_ProteinLevel_AllFeatures.png",
    r"D:\博士\Phospho\Fig_ProteinLevel_AllFeatures.tiff",
    r"D:\博士\Phospho\Fig_ProteinLevel_Confounding_v3.png",
    r"D:\博士\Phospho\Fig_ProteinLevel_Confounding_v3.tiff",
    r"D:\博士\Phospho\Fig_WithinProtein_Paired.png",
    r"D:\博士\Phospho\Fig_WithinProtein_Paired.tiff",
]

FEATURE_KEYS = {
    "backbone": "Backbone",
    "sidechain": "Sidechain",
    "disoMine": "Disorder",
    "helix": "Helix",
    "sheet": "Sheet",
    "coil": "Coil",
    "earlyFolding": "EarlyFolding",
}
FEATURES = list(FEATURE_KEYS.values())
VALUE_COLUMNS = [f"avg_{x}" for x in FEATURES]


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Cliff's delta from Mann-Whitney U; positive means mobile > static."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    stat, pvalue = mannwhitneyu(a, b, alternative="two-sided")
    delta = 2 * stat / (len(a) * len(b)) - 1
    return float(delta), float(pvalue)


def ols_fit(y: np.ndarray, x: np.ndarray) -> dict[str, float | np.ndarray]:
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    n, p = x.shape
    rss = float(np.sum(resid**2))
    df_resid = n - p
    sigma2 = rss / df_resid
    xtx_inv = np.linalg.inv(x.T @ x)
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    tvals = beta / se
    pvals = 2 * t.sf(np.abs(tvals), df_resid)
    return {
        "beta": beta,
        "resid": resid,
        "pvals": pvals,
        "rss": rss,
        "df_resid": df_resid,
    }


def residualize(df: pd.DataFrame, col: str, controls: list[str]) -> pd.Series:
    sub = df[["group", col] + controls].dropna()
    x = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy(float) for c in controls])
    y = sub[col].to_numpy(float)
    fit = ols_fit(y, x)
    out = pd.Series(np.nan, index=df.index, dtype=float)
    out.loc[sub.index] = fit["resid"]
    return out


def adjusted_ols_group_effect(df: pd.DataFrame, col: str, controls: list[str]) -> tuple[float, float]:
    """Return adjusted group beta and partial eta squared for group."""
    sub = df[["group", col] + controls].dropna()
    group = (sub["group"] == "mobile").astype(float).to_numpy()
    x_full = np.column_stack(
        [np.ones(len(sub)), group] + [sub[c].to_numpy(float) for c in controls]
    )
    x_red = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy(float) for c in controls])
    y = sub[col].to_numpy(float)
    full = ols_fit(y, x_full)
    red = ols_fit(y, x_red)
    ss_group = red["rss"] - full["rss"]
    partial_eta2 = ss_group / (ss_group + full["rss"])
    group_p = float(full["pvals"][1])
    return float(partial_eta2), group_p


def permutation_p(values: np.ndarray, labels: np.ndarray, observed_delta: float, n_perm: int = 5000) -> float:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    ok = np.isfinite(values)
    values = values[ok]
    labels = labels[ok].copy()
    ranks = rankdata(values, method="average")
    n = len(values)
    n_mobile = int(np.sum(labels == "mobile"))
    n_static = n - n_mobile
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        RNG.shuffle(labels)
        rank_sum_mobile = float(np.sum(ranks[labels == "mobile"]))
        u_mobile = rank_sum_mobile - n_mobile * (n_mobile + 1) / 2
        null[i] = 2 * u_mobile / (n_mobile * n_static) - 1
    return float((np.sum(np.abs(null) >= abs(observed_delta)) + 1) / (n_perm + 1))


def matched_indices(df: pd.DataFrame, caliper: float = 0.25) -> tuple[list[int], list[int]]:
    """Greedy 1:1 nearest-neighbour match on standardized log length and pct_dis."""
    match_cols = ["log10_length", "pct_dis"]
    d2 = df.dropna(subset=match_cols)
    z = (d2[match_cols] - d2[match_cols].mean()) / d2[match_cols].std(ddof=0)
    mobile_idx = d2[d2["group"] == "mobile"].index.tolist()
    static_pool = d2[d2["group"] == "static"].index.tolist()
    z_static = z.loc[static_pool].to_numpy(float)
    used: set[int] = set()
    pairs: list[tuple[int, int]] = []
    for mi in mobile_idx:
        distances = np.sqrt(np.sum((z_static - z.loc[mi].to_numpy(float)) ** 2, axis=1))
        for order_i in np.argsort(distances):
            si = static_pool[order_i]
            if si in used:
                continue
            if distances[order_i] > caliper:
                break
            used.add(si)
            pairs.append((mi, si))
            break
    return [m for m, _ in pairs], [s for _, s in pairs]


def reconstruct_table() -> pd.DataFrame:
    gene_to_uniprot = json.loads(GENE_TO_UNIPROT.read_text(encoding="utf-8"))
    high = pd.read_csv(HIGH_SITES)
    low = pd.read_csv(LOW_SITES)
    high_genes = set(high["Gene"].dropna().astype(str))
    low_genes = set(low["Gene"].dropna().astype(str))

    mobile_genes = high_genes & low_genes
    static_genes = low_genes - high_genes

    rows: list[dict[str, float | str | int]] = []
    for gene in sorted(mobile_genes | static_genes):
        uid = gene_to_uniprot.get(gene)
        if not uid:
            continue
        json_path = BIOPHYS / f"{uid}.json"
        if not json_path.exists():
            continue
        data = json.loads(json_path.read_text(encoding="utf-8"))
        residues = data.get("residues", [])
        if not residues:
            continue
        row: dict[str, float | str | int] = {
            "gene": gene,
            "uniprot": uid,
            "group": "mobile" if gene in mobile_genes else "static",
            "length": len(residues),
        }
        for json_key, label in FEATURE_KEYS.items():
            vals = [r.get(json_key) for r in residues if r.get(json_key) is not None]
            row[f"avg_{label}"] = float(np.mean(vals)) if vals else math.nan
            if json_key == "disoMine":
                row["pct_dis"] = float(100 * sum(v > 0.5 for v in vals) / len(residues))
        rows.append(row)

    df = pd.DataFrame(rows)
    df["log10_length"] = np.log10(df["length"])
    return df


def main() -> None:
    OUT.mkdir(exist_ok=True)
    df = reconstruct_table()
    df.to_csv(OUT / "per_protein_table_v2.csv", index=False)

    source_lines = [
        "Source-data locator",
        "===================",
        "",
        "Figures located:",
    ]
    source_lines.extend([f"- {p} | exists={Path(p).exists()}" for p in FIGURE_PATHS])
    source_lines.extend(
        [
            "",
            "Original figure script:",
            "- D:\\博士\\Phospho\\0331.py | lines around Cell 1 reconstruct prot_df and save Fig_ProteinLevel_AllFeatures; lines around Cell 2 save Fig_ProteinLevel_Confounding_v3 and Fig_WithinProtein_Paired.",
            "",
            "Underlying per-protein table:",
            "- Original prot_df was not saved by 0331.py.",
            f"- Reconstructed table: {OUT / 'per_protein_table_v2.csv'}",
            f"- Inputs: {HIGH_SITES.name}, {LOW_SITES.name}, {GENE_TO_UNIPROT.name}, biophys_json/*.json.",
            f"- Protein rows: {len(df)}; group counts: {df['group'].value_counts().to_dict()}",
            "- Unit of analysis: one row per protein.",
        ]
    )
    (OUT / "source_locator_v2.md").write_text("\n".join(source_lines), encoding="utf-8")

    group_summary = []
    mobile = df[df["group"] == "mobile"]
    static = df[df["group"] == "static"]
    for col in ["length", "pct_dis"]:
        d, p = cliffs_delta(mobile[col].to_numpy(), static[col].to_numpy())
        group_summary.append(
            {
                "variable": col,
                "mobile_n": int(mobile[col].notna().sum()),
                "static_n": int(static[col].notna().sum()),
                "mobile_median": float(mobile[col].median()),
                "static_median": float(static[col].median()),
                "cliff_delta_mobile_minus_static": d,
                "mannwhitney_p": p,
            }
        )
    pd.DataFrame(group_summary).to_csv(OUT / "confounder_group_differences_v2.csv", index=False)

    matched_mobile, matched_static = matched_indices(df)
    results = []
    for feature in FEATURES:
        col = f"avg_{feature}"
        mob_values = mobile[col].dropna().to_numpy(float)
        sta_values = static[col].dropna().to_numpy(float)
        raw_delta, raw_p = cliffs_delta(mob_values, sta_values)

        labels = df.loc[df[col].notna(), "group"].to_numpy()
        raw_perm = permutation_p(df.loc[df[col].notna(), col].to_numpy(float), labels, raw_delta)

        primary_controls = ["log10_length", "pct_dis"]
        note = "length+pct_dis"
        if feature == "Disorder":
            # avg_Disorder and pct_dis are two summaries of the same predicted disorder track.
            # Keep a length-only control as the interpretable primary view, and report
            # length+pct_dis as an explicit overcontrolled sensitivity analysis.
            primary_controls = ["log10_length"]
            note = "length only; length+pct_dis sensitivity also reported"

        resid = residualize(df, col, primary_controls)
        res_mobile = resid[df["group"] == "mobile"].dropna().to_numpy(float)
        res_static = resid[df["group"] == "static"].dropna().to_numpy(float)
        resid_delta, resid_p = cliffs_delta(res_mobile, res_static)
        resid_labels = df.loc[resid.notna(), "group"].to_numpy()
        resid_perm = permutation_p(resid.loc[resid.notna()].to_numpy(float), resid_labels, resid_delta)
        eta2, ols_p = adjusted_ols_group_effect(df, col, primary_controls)

        matched_delta, matched_p = cliffs_delta(
            df.loc[matched_mobile, col].dropna().to_numpy(float),
            df.loc[matched_static, col].dropna().to_numpy(float),
        )

        sens_delta = sens_p = sens_eta2 = sens_ols_p = math.nan
        if feature == "Disorder":
            sens = residualize(df, col, ["log10_length", "pct_dis"])
            sens_delta, sens_p = cliffs_delta(
                sens[df["group"] == "mobile"].dropna().to_numpy(float),
                sens[df["group"] == "static"].dropna().to_numpy(float),
            )
            sens_eta2, sens_ols_p = adjusted_ols_group_effect(df, col, ["log10_length", "pct_dis"])

        results.append(
            {
                "feature": feature,
                "mobile_n": int(mobile[col].notna().sum()),
                "static_n": int(static[col].notna().sum()),
                "mobile_median": float(mobile[col].median()),
                "static_median": float(static[col].median()),
                "raw_cliff_delta": raw_delta,
                "raw_mannwhitney_p": raw_p,
                "raw_permutation_p": raw_perm,
                "controlled_for": "+".join(primary_controls),
                "residualized_cliff_delta": resid_delta,
                "residualized_mannwhitney_p": resid_p,
                "residualized_permutation_p": resid_perm,
                "adjusted_ols_partial_eta2_group": eta2,
                "adjusted_ols_group_p": ols_p,
                "matched_pairs_n": len(matched_mobile),
                "matched_cliff_delta": matched_delta,
                "matched_mannwhitney_p": matched_p,
                "disorder_overcontrolled_len_pctdis_cliff_delta": sens_delta,
                "disorder_overcontrolled_len_pctdis_p": sens_p,
                "disorder_overcontrolled_len_pctdis_partial_eta2": sens_eta2,
                "disorder_overcontrolled_len_pctdis_ols_p": sens_ols_p,
                "note": note,
            }
        )

    out = pd.DataFrame(results)
    out.to_csv(OUT / "scaffold_confounder_results_v2.csv", index=False)

    print((OUT / "source_locator_v2.md").read_text(encoding="utf-8"))
    print("\nConfounder group differences:")
    print(pd.DataFrame(group_summary).to_string(index=False))
    print("\nBefore/after scaffold results:")
    display_cols = [
        "feature",
        "raw_cliff_delta",
        "raw_mannwhitney_p",
        "residualized_cliff_delta",
        "residualized_mannwhitney_p",
        "adjusted_ols_partial_eta2_group",
        "matched_pairs_n",
        "matched_cliff_delta",
        "matched_mannwhitney_p",
    ]
    print(out[display_cols].to_string(index=False))
    print("\nSaved v2 outputs under scaffold_check/.")


if __name__ == "__main__":
    main()
