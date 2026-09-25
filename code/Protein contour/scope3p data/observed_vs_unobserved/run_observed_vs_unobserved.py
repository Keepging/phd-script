import sys, os, re, glob, json, csv, math, gzip, shutil
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, mannwhitneyu


JSON_DIR, EVIDENCE_PATH, OUT = sys.argv[1:4]
TABLES = os.path.join(OUT, "tables")
FIGURES = os.path.join(OUT, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)
shutil.rmtree(os.path.join(OUT, "_pooled_binary_tmp"), ignore_errors=True)

FEATURES = ["backbone", "sidechain", "disoMine", "helix", "sheet", "coil", "earlyFolding"]
RESIDUES = {"S", "T", "Y"}
MAIN = "main_n_projects>=1"
SENS = "sensitivity_n_projects>=3"
STRATA = [f"residue_{aa}_n_projects>=1" for aa in "STY"]
ANALYSES = [MAIN, SENS] + STRATA


def canonical_accession(value):
    return re.sub(r"-\d+$", "", str(value).strip().upper())


def finite_float(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else np.nan
    except Exception:
        return np.nan


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    result = np.full(len(pvalues), np.nan)
    valid = np.where(np.isfinite(pvalues))[0]
    if not len(valid):
        return result
    values = pvalues[valid]
    order = np.argsort(values)
    ranked = values[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty(m)
    restored[order] = adjusted
    result[valid] = restored
    return result


if not os.path.isfile(EVIDENCE_PATH):
    raise RuntimeError("Core assumption failed: observed-site evidence table is missing")

evidence = pd.read_csv(EVIDENCE_PATH, encoding="utf-8-sig", low_memory=False)
required_evidence = {"accession", "position", "residue", "n_projects"}
if not required_evidence.issubset(evidence.columns):
    raise RuntimeError(f"Core assumption failed: evidence columns missing: {sorted(required_evidence-set(evidence.columns))}")
evidence["accession"] = evidence["accession"].astype("string").map(canonical_accession)
evidence["position"] = pd.to_numeric(evidence["position"], errors="coerce")
evidence["residue"] = evidence["residue"].astype("string").str.strip().str.upper()
evidence["n_projects"] = pd.to_numeric(evidence["n_projects"], errors="coerce").fillna(0).astype(int)
valid_evidence = (
    evidence["accession"].notna() & evidence["accession"].ne("") &
    evidence["position"].notna() & np.isfinite(evidence["position"]) &
    (np.floor(evidence["position"]) == evidence["position"]) & (evidence["position"] > 0) &
    evidence["residue"].isin(RESIDUES) & evidence["n_projects"].ge(1)
)
if not valid_evidence.all():
    raise RuntimeError(f"Core assumption failed: {(~valid_evidence).sum()} invalid evidence keys")
evidence["position"] = evidence["position"].astype(int)
if evidence.duplicated(["accession", "position", "residue"]).any():
    raise RuntimeError("Core assumption failed: evidence keys are not unique")

evidence_by_accession = defaultdict(lambda: defaultdict(dict))
for row in evidence[["accession", "position", "residue", "n_projects"]].itertuples(index=False):
    evidence_by_accession[row.accession][row.position][row.residue] = row.n_projects

json_files = sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
if not json_files:
    raise RuntimeError("Core assumption failed: no JSON files found")

candidate_temp = os.path.join(TABLES, "_candidate_sites_labeled.tmp.csv")
candidate_csv = os.path.join(TABLES, "candidate_sites_labeled.csv")
candidate_gz = candidate_csv + ".gz"
for path in [candidate_temp, candidate_csv, candidate_gz]:
    if os.path.exists(path):
        os.remove(path)

pairs = {analysis: {feature: [] for feature in FEATURES} for analysis in ANALYSES}
main_exclusions = Counter()
sensitivity_exclusions = Counter()
stratum_qualifying = Counter()
qualifying_site_counts = []
warnings_log = []

counts = Counter()
feature_missing = Counter()


def add_pairs(analysis, candidates, observed_rule, unobserved_rule):
    observed = [row for row in candidates if observed_rule(row)]
    unobserved = [row for row in candidates if unobserved_rule(row)]
    if not observed or not unobserved:
        return
    for index, feature in enumerate(FEATURES):
        observed_values = np.asarray([row[2][index] for row in observed], dtype=float)
        unobserved_values = np.asarray([row[2][index] for row in unobserved], dtype=float)
        observed_values = observed_values[np.isfinite(observed_values)]
        unobserved_values = unobserved_values[np.isfinite(unobserved_values)]
        if len(observed_values) and len(unobserved_values):
            pairs[analysis][feature].append((np.median(observed_values), np.median(unobserved_values)))


with open(candidate_temp, "w", newline="", encoding="utf-8-sig") as output_handle:
    writer = csv.writer(output_handle)
    writer.writerow(["accession", "position", "residue"] + FEATURES + ["observed", "n_projects"])

    for file_number, json_path in enumerate(json_files, start=1):
        accession = canonical_accession(os.path.splitext(os.path.basename(json_path))[0])
        try:
            with open(json_path, "r", encoding="utf-8-sig") as json_handle:
                payload = json.load(json_handle)
            residues = payload.get("residues")
            if not isinstance(residues, list):
                raise ValueError("residues is not a list")
        except Exception as error:
            counts["malformed_files"] += 1
            if len(warnings_log) < 10:
                warnings_log.append(f"Skipped {os.path.basename(json_path)}: {type(error).__name__}: {error}")
            continue

        counts["parsed_files"] += 1
        if payload.get("protein_id") is not None and canonical_accession(payload["protein_id"]) != accession:
            counts["protein_id_mismatch_files"] += 1

        position_to_residue = {}
        protein_candidates = []
        candidate_rows = []
        counts["total_residue_records"] += len(residues)

        for record in residues:
            if not isinstance(record, dict):
                counts["invalid_residue_records"] += 1
                continue
            try:
                position_number = float(record.get("seqpos"))
                if not math.isfinite(position_number) or position_number <= 0 or math.floor(position_number) != position_number:
                    raise ValueError
                position = int(position_number)
            except Exception:
                counts["invalid_residue_records"] += 1
                continue
            residue = str(record.get("aa", "")).strip().upper()
            if len(residue) != 1:
                counts["invalid_residue_records"] += 1
                continue
            if position in position_to_residue:
                counts["duplicate_seqpos_records"] += 1
                continue
            position_to_residue[position] = residue
            if residue not in RESIDUES:
                continue

            values = tuple(finite_float(record.get(feature)) for feature in FEATURES)
            for index, feature in enumerate(FEATURES):
                if not np.isfinite(values[index]):
                    feature_missing[feature] += 1
            n_projects = int(evidence_by_accession.get(accession, {}).get(position, {}).get(residue, 0))
            observed = n_projects >= 1
            counts["candidate_sites"] += 1
            counts["observed_candidate_sites"] += int(observed)
            counts["unobserved_candidate_sites"] += int(not observed)
            protein_candidates.append((position, residue, values, n_projects))
            candidate_rows.append(
                [accession, position, residue] +
                ["" if not np.isfinite(value) else f"{value:.10g}" for value in values] +
                [bool(observed), n_projects]
            )
        writer.writerows(candidate_rows)

        for position, residue_map in evidence_by_accession.get(accession, {}).items():
            for evidence_residue in residue_map:
                counts["evidence_sites_on_parsed_accessions"] += 1
                if position not in position_to_residue:
                    counts["evidence_positions_missing"] += 1
                elif position_to_residue[position] == evidence_residue:
                    counts["evidence_exact_matches"] += 1
                else:
                    counts["evidence_residue_mismatches"] += 1

        observed_count = sum(row[3] >= 1 for row in protein_candidates)
        unobserved_count = sum(row[3] == 0 for row in protein_candidates)
        if not protein_candidates:
            main_exclusions["no_candidate_STY"] += 1
        elif observed_count == 0:
            main_exclusions["no_observed_site"] += 1
        elif unobserved_count == 0:
            main_exclusions["no_unobserved_candidate_site"] += 1
        else:
            main_exclusions["qualifying"] += 1
            counts["qualifying_observed_sites"] += observed_count
            counts["qualifying_unobserved_sites"] += unobserved_count
            qualifying_site_counts.append((accession, observed_count, unobserved_count))
            add_pairs(MAIN, protein_candidates, lambda row: row[3] >= 1, lambda row: row[3] == 0)

        strict_observed_count = sum(row[3] >= 3 for row in protein_candidates)
        middle_count = sum(1 <= row[3] <= 2 for row in protein_candidates)
        if not protein_candidates:
            sensitivity_exclusions["no_candidate_STY"] += 1
        elif strict_observed_count == 0:
            sensitivity_exclusions["no_observed_site_ge3"] += 1
        elif unobserved_count == 0:
            sensitivity_exclusions["no_unobserved_candidate_site"] += 1
        else:
            sensitivity_exclusions["qualifying"] += 1
            counts["sensitivity_observed_sites"] += strict_observed_count
            counts["sensitivity_unobserved_sites"] += unobserved_count
            counts["sensitivity_excluded_1_2_sites"] += middle_count
            add_pairs(SENS, protein_candidates, lambda row: row[3] >= 3, lambda row: row[3] == 0)

        for residue_type in "STY":
            subset = [row for row in protein_candidates if row[1] == residue_type]
            if any(row[3] >= 1 for row in subset) and any(row[3] == 0 for row in subset):
                stratum_qualifying[residue_type] += 1
                add_pairs(
                    f"residue_{residue_type}_n_projects>=1", subset,
                    lambda row: row[3] >= 1, lambda row: row[3] == 0,
                )

        if file_number % 500 == 0 or file_number == len(json_files):
            output_handle.flush()
            print(
                f"[parse] {file_number:,}/{len(json_files):,} parsed={counts['parsed_files']:,} "
                f"malformed={counts['malformed_files']:,} candidates={counts['candidate_sites']:,}",
                flush=True,
            )

addressable = counts["evidence_exact_matches"] + counts["evidence_residue_mismatches"]
mismatch_rate = counts["evidence_residue_mismatches"] / addressable if addressable else 1.0
overlap_rate = counts["evidence_exact_matches"] / len(evidence)
if addressable == 0 or counts["evidence_exact_matches"] < max(100, int(0.01 * len(evidence))):
    os.remove(candidate_temp)
    raise RuntimeError(f"Core assumption failed: near-zero key overlap ({counts['evidence_exact_matches']:,}/{len(evidence):,}, {overlap_rate:.2%})")
if mismatch_rate > 0.05:
    os.remove(candidate_temp)
    raise RuntimeError(f"Core assumption failed: residue mismatch rate {mismatch_rate:.2%} exceeds 5%")

os.replace(candidate_temp, candidate_csv)
candidate_path = candidate_csv
candidate_size = os.path.getsize(candidate_csv)
candidate_compressed = False
if candidate_size > 500 * 1024 * 1024:
    with open(candidate_csv, "rb") as source, gzip.open(candidate_gz, "wb", compresslevel=6) as destination:
        shutil.copyfileobj(source, destination, length=1024 * 1024)
    os.remove(candidate_csv)
    candidate_path = candidate_gz
    candidate_compressed = True


def calculate_paired_stats(analysis):
    rows = []
    for feature in FEATURES:
        values = np.asarray(pairs[analysis][feature], dtype=float)
        if values.size == 0:
            rows.append({
                "analysis": analysis, "feature": feature, "n_proteins": 0,
                "median_observed": np.nan, "median_unobserved": np.nan,
                "median_paired_diff": np.nan, "frac_proteins_observed_higher": np.nan,
                "wilcoxon_stat": np.nan, "p_raw": np.nan,
            })
            continue
        observed_values, unobserved_values = values[:, 0], values[:, 1]
        differences = observed_values - unobserved_values
        if np.allclose(differences, 0):
            statistic, pvalue = 0.0, 1.0
        else:
            test = wilcoxon(observed_values, unobserved_values, alternative="two-sided", zero_method="wilcox", method="auto")
            statistic, pvalue = float(test.statistic), float(test.pvalue)
        rows.append({
            "analysis": analysis, "feature": feature, "n_proteins": len(values),
            "median_observed": float(np.median(observed_values)),
            "median_unobserved": float(np.median(unobserved_values)),
            "median_paired_diff": float(np.median(differences)),
            "frac_proteins_observed_higher": float(np.mean(differences > 0)),
            "wilcoxon_stat": statistic, "p_raw": pvalue,
        })
    for row, corrected in zip(rows, bh_adjust([row["p_raw"] for row in rows])):
        row["p_BH"] = corrected
    return rows


main_rows = calculate_paired_stats(MAIN)
sensitivity_rows = calculate_paired_stats(SENS)
stratified_rows = []
for analysis in STRATA:
    stratified_rows.extend(calculate_paired_stats(analysis))

pd.DataFrame(main_rows).drop(columns=["analysis"]).to_csv(
    os.path.join(TABLES, "within_protein_paired_stats.csv"), index=False, encoding="utf-8-sig"
)
pd.DataFrame(sensitivity_rows + stratified_rows).to_csv(
    os.path.join(TABLES, "sensitivity_and_stratified_stats.csv"), index=False, encoding="utf-8-sig"
)

pooled_rows = []
qualifying_accessions = {row[0] for row in qualifying_site_counts}
for feature in FEATURES:
    observed_parts, unobserved_parts = [], []
    for chunk in pd.read_csv(
        candidate_path, usecols=["accession", feature, "observed"], chunksize=500000,
        encoding="utf-8-sig", low_memory=False,
    ):
        values = pd.to_numeric(chunk[feature], errors="coerce").to_numpy(float)
        observed_mask = chunk["observed"].astype(str).str.lower().eq("true").to_numpy()
        qualifying_mask = chunk["accession"].isin(qualifying_accessions).to_numpy()
        observed_parts.append(values[qualifying_mask & observed_mask & np.isfinite(values)])
        unobserved_parts.append(values[qualifying_mask & (~observed_mask) & np.isfinite(values)])
    observed_values = np.concatenate(observed_parts)
    unobserved_values = np.concatenate(unobserved_parts)
    test = mannwhitneyu(observed_values, unobserved_values, alternative="two-sided", method="asymptotic")
    pooled_rows.append({
        "analysis": "pooled_site_level_confounded_reference", "feature": feature,
        "n_observed_sites": len(observed_values),
        "n_unobserved_candidate_sites": len(unobserved_values),
        "median_observed": float(np.median(observed_values)),
        "median_unobserved": float(np.median(unobserved_values)),
        "median_difference_observed_minus_unobserved": float(np.median(observed_values) - np.median(unobserved_values)),
        "mannwhitney_U": float(test.statistic), "p_raw": float(test.pvalue),
    })
    print(f"[pooled] {feature} observed={len(observed_values):,} unobserved={len(unobserved_values):,}", flush=True)
for row, corrected in zip(pooled_rows, bh_adjust([row["p_raw"] for row in pooled_rows])):
    row["p_BH_reference"] = corrected
pd.DataFrame(pooled_rows).to_csv(
    os.path.join(TABLES, "pooled_sitelevel_reference.csv"), index=False, encoding="utf-8-sig"
)

# Static figure contracts: distributions, paired relationship, and within-protein candidate counts.
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 13,
    "axes.labelsize": 10, "axes.edgecolor": "#4A5568", "axes.labelcolor": "#1F2937",
    "xtick.color": "#374151", "ytick.color": "#374151",
})
blue, ink, grid = "#2F6B9A", "#1F2937", "#D9DEE5"

difference_data = []
for feature in FEATURES:
    values = np.asarray(pairs[MAIN][feature], dtype=float)
    difference_data.append(values[:, 0] - values[:, 1])
fig, ax = plt.subplots(figsize=(8.2, 5.1), facecolor="white")
violins = ax.violinplot(
    difference_data, positions=np.arange(1, len(FEATURES) + 1), vert=False,
    showmeans=False, showmedians=True, showextrema=False, widths=0.8,
)
for body in violins["bodies"]:
    body.set_facecolor(blue)
    body.set_edgecolor(ink)
    body.set_alpha(0.58)
violins["cmedians"].set_color(ink)
violins["cmedians"].set_linewidth(1.8)
ax.axvline(0, color=ink, linestyle="--", linewidth=1.2)
ax.set_yticks(np.arange(1, len(FEATURES) + 1), FEATURES)
ax.set_xlabel("Per-protein median difference (observed - unobserved candidate)")
ax.set_ylabel("Predicted feature")
ax.set_title("Within-protein paired differences by predicted feature", loc="left", color=ink, pad=10)
ax.grid(axis="x", color=grid, linewidth=0.7, alpha=0.8)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "fig_paired_diff_by_feature.png"), dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)

disorder_pairs = np.asarray(pairs[MAIN]["disoMine"], dtype=float)
disorder_row = next(row for row in main_rows if row["feature"] == "disoMine")
fig, ax = plt.subplots(figsize=(6.4, 5.5), facecolor="white")
ax.scatter(disorder_pairs[:, 1], disorder_pairs[:, 0], s=10, alpha=0.24, color=blue, edgecolors="none")
low, high = float(np.nanmin(disorder_pairs)), float(np.nanmax(disorder_pairs))
padding = 0.03 * (high - low if high > low else 1)
ax.plot([low-padding, high+padding], [low-padding, high+padding], color=ink, linestyle="--", linewidth=1.2, label="Equal medians")
ax.set_xlim(low-padding, high+padding)
ax.set_ylim(low-padding, high+padding)
ax.set_xlabel("Median disoMine: unobserved candidate sites")
ax.set_ylabel("Median disoMine: observed sites")
ax.set_title("Paired within-protein disorder predictions", loc="left", color=ink, pad=10)
display_p = "<1e-300" if disorder_row["p_BH"] == 0 else (f"{disorder_row['p_BH']:.2e}" if disorder_row["p_BH"] < 0.001 else f"{disorder_row['p_BH']:.3f}")
p_label = f"BH p{display_p}" if display_p.startswith("<") else f"BH p={display_p}"
ax.text(0.03, 0.97, f"n={len(disorder_pairs):,} proteins; {p_label}", transform=ax.transAxes, ha="left", va="top", color=ink)
ax.grid(color=grid, linewidth=0.7, alpha=0.7)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "fig_disorder_paired.png"), dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)

site_counts = np.asarray([(row[1], row[2]) for row in qualifying_site_counts], dtype=float)
fig, ax = plt.subplots(figsize=(6.7, 5.5), facecolor="white")
ax.scatter(site_counts[:, 1], site_counts[:, 0], s=10, alpha=0.25, color=blue, edgecolors="none")
ax.set_xscale("log")
ax.set_yscale("log")
minimum = max(1, float(site_counts.min()))
maximum = float(site_counts.max())
ax.plot([minimum, maximum], [minimum, maximum], color=ink, linestyle="--", linewidth=1.2, label="Equal counts")
ax.set_xlabel("Unobserved candidate S/T/Y sites per protein")
ax.set_ylabel("Observed S/T/Y sites per protein")
ax.set_title("Candidate-site counts within qualifying proteins", loc="left", color=ink, pad=10)
ax.text(0.03, 0.97, f"n={len(site_counts):,} proteins", transform=ax.transAxes, ha="left", va="top", color=ink)
ax.grid(which="both", color=grid, linewidth=0.6, alpha=0.6)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "fig_candidate_counts.png"), dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)

main_lookup = {row["feature"]: row for row in main_rows}
sensitivity_lookup = {row["feature"]: row for row in sensitivity_rows}
main_significant = [feature for feature in FEATURES if main_lookup[feature]["p_BH"] < 0.05]
sensitivity_significant = [feature for feature in FEATURES if sensitivity_lookup[feature]["p_BH"] < 0.05]
sensitivity_direction_agreement = sum(
    np.sign(main_lookup[feature]["median_paired_diff"]) == np.sign(sensitivity_lookup[feature]["median_paired_diff"])
    for feature in FEATURES
)
stratum_direction_agreement, stratum_significant = {}, {}
for residue_type in "STY":
    lookup = {
        row["feature"]: row for row in stratified_rows
        if row["analysis"] == f"residue_{residue_type}_n_projects>=1"
    }
    stratum_direction_agreement[residue_type] = sum(
        np.sign(main_lookup[feature]["median_paired_diff"]) == np.sign(lookup[feature]["median_paired_diff"])
        for feature in FEATURES
    )
    stratum_significant[residue_type] = [feature for feature in FEATURES if lookup[feature]["p_BH"] < 0.05]


def format_pvalue(value):
    if value == 0:
        return "<1e-300"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3f}"


def feature_sentence(row):
    difference = row["median_paired_diff"]
    return (
        f"- For **{row['feature']}**, observed minus unobserved candidate was **{difference:.4f}**; "
        f"**{100*row['frac_proteins_observed_higher']:.1f}%** of proteins had a higher observed-site median; "
        f"BH p was **{format_pvalue(row['p_BH'])}**."
    )


candidate_filename = os.path.basename(candidate_path)
compression_sentence = (
    f"The candidate table exceeded 500 MB and was written as {candidate_filename}."
    if candidate_compressed else
    f"The candidate table was {candidate_size/1024/1024:.1f} MB and was written as {candidate_filename}."
)
main_significant_text = ", ".join(main_significant) if main_significant else "none"
sensitivity_significant_text = ", ".join(sensitivity_significant) if sensitivity_significant else "none"
stratum_text = "; ".join(
    f"{residue_type}: {stratum_direction_agreement[residue_type]}/7 directions agreed, BH-significant features "
    f"{', '.join(stratum_significant[residue_type]) if stratum_significant[residue_type] else 'none'}"
    for residue_type in "STY"
)

findings = f"""# Observed versus unobserved candidate phosphosites within proteins

## Technical summary

This analysis compared **observed sites** with **unobserved candidate sites** only within the same protein. An observed site was an S/T/Y residue present in the 116-project evidence table with n_projects at least 1. An unobserved candidate site was another S/T/Y residue on that same protein with n_projects equal to 0.

The JSON collection yielded **{counts['candidate_sites']:,} candidate S/T/Y sites from {counts['parsed_files']:,} parsed proteins**. The main paired analysis retained **{main_exclusions['qualifying']:,} proteins**, containing **{counts['qualifying_observed_sites']:,} observed sites** and **{counts['qualifying_unobserved_sites']:,} unobserved candidate sites**. It excluded **{main_exclusions['no_observed_site']:,} proteins with no observed site**, **{main_exclusions['no_unobserved_candidate_site']:,} with no unobserved candidate site**, and **{main_exclusions['no_candidate_STY']:,} with no candidate S/T/Y residue**.

## The seven within-protein results

Each protein contributed the median feature value among its observed sites and the median among its unobserved candidate sites. The reported difference is observed minus unobserved candidate.

{chr(10).join(feature_sentence(row) for row in main_rows)}

The BH-significant main features were **{main_significant_text}**. These tests describe consistent within-protein associations; they do not establish a causal effect.

## Sensitivity and residue-type checks

With the stricter threshold n_projects at least 3, sites seen in only 1-2 projects were excluded rather than relabelled. This sensitivity analysis retained **{sensitivity_exclusions['qualifying']:,} proteins**, **{counts['sensitivity_observed_sites']:,} observed sites**, and **{counts['sensitivity_unobserved_sites']:,} unobserved candidate sites**; **{counts['sensitivity_excluded_1_2_sites']:,} sites with n_projects 1-2 were excluded within those proteins**. The direction agreed with the main analysis for **{sensitivity_direction_agreement}/7 features**. Its BH-significant features were **{sensitivity_significant_text}**.

The residue-stratified checks retained **{stratum_qualifying['S']:,} S proteins**, **{stratum_qualifying['T']:,} T proteins**, and **{stratum_qualifying['Y']:,} Y proteins**. Relative to the main direction, **{stratum_text}**. Full effect sizes and corrected p-values are in the stratified table.

## Data checks and interpretation

The canonical-key overlap matched **{counts['evidence_exact_matches']:,} evidence sites exactly**. At positions assessable in parsed JSONs, **{counts['evidence_residue_mismatches']:,} residue letters disagreed out of {addressable:,} ({100*mismatch_rate:.3f}%)**, below the 5% stop threshold; mismatches were not corrected. Another **{counts['evidence_positions_missing']:,} evidence positions on parsed accessions were absent from their JSON residue lists**, and **{len(evidence)-counts['evidence_sites_on_parsed_accessions']:,} evidence sites belonged to accessions without a JSON file**.

The pooled site-level comparison is supplied only as a **confounded reference** because it ignores protein identity. The within-protein paired analysis is the core result. {compression_sentence}

## Caveats

- **Unobserved candidate does not mean non-phosphorylated.** Sampling, abundance, digestion, detectability, and identification confidence all affect observation.
- The results are **descriptive associations, not causal effects**.
- The seven biophysical values are **predictions, not measurements**.
- {counts['malformed_files']:,} unreadable or malformed JSON files were skipped; feature-specific missing predictions were omitted only from the affected feature pair.
- Extremely small p-values that underflowed double precision are displayed as **less than 1e-300** rather than as zero.
"""

findings_path = os.path.join(OUT, "FINDINGS_observed_vs_unobserved.md")
with open(findings_path, "w", encoding="utf-8") as handle:
    handle.write(findings)

print("\n[run log]", flush=True)
for key in [
    "parsed_files", "malformed_files", "protein_id_mismatch_files", "invalid_residue_records",
    "duplicate_seqpos_records", "candidate_sites", "observed_candidate_sites",
    "unobserved_candidate_sites", "qualifying_observed_sites", "qualifying_unobserved_sites",
    "evidence_exact_matches", "evidence_residue_mismatches", "evidence_positions_missing",
]:
    print(f"{key}={counts[key]}", flush=True)
print(f"json_files_found={len(json_files)}", flush=True)
print(f"proteins_qualifying_main={main_exclusions['qualifying']}", flush=True)
print(f"proteins_excluded_no_observed={main_exclusions['no_observed_site']}", flush=True)
print(f"proteins_excluded_no_unobserved_candidate={main_exclusions['no_unobserved_candidate_site']}", flush=True)
print(f"proteins_excluded_no_candidate_STY={main_exclusions['no_candidate_STY']}", flush=True)
print(f"residue_mismatch_rate={mismatch_rate:.6%}", flush=True)
print("warnings:", flush=True)
if warnings_log:
    for warning in warnings_log:
        print(f"- {warning}", flush=True)
else:
    print("- none", flush=True)
print("\n[FINDINGS_observed_vs_unobserved.md]\n", flush=True)
print(findings, flush=True)
