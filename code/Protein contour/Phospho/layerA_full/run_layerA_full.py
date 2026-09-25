import csv
import io
import json
import math
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests


ROOT = r"D:\博士\Protein contour\Phospho"
OUT = os.path.join(ROOT, "layerA_full")
BIO = os.path.join(ROOT, "biophys_json")

HPA_URL = "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
SCOP3P_URL = "https://iomics.ugent.be/scop3p/api/get-all-modifications"
UNIPROT_IDMAPPING = "https://rest.uniprot.org/idmapping"

HPA_ZIP = os.path.join(OUT, "subcellular_location_25.1.tsv.zip")
HPA_TSV = os.path.join(OUT, "subcellular_location_25.1.tsv")
SCOP3P_JSON = os.path.join(OUT, "scop3p_all_modifications.json")
ENSEMBL_MAP_JSON = os.path.join(OUT, "uniprot_idmapping_ensembl_v2.json")
GENE_MAP_JSON = os.path.join(OUT, "uniprot_idmapping_gene_name_fallback_v2.json")
LOG_PATH = os.path.join(OUT, "run_log.txt")

FEATURES = [
    ("backbone", "backbone_dynamics"),
    ("sidechain", "sidechain_dynamics"),
    ("disoMine", "disorder_propensity"),
    ("helix", "helix_propensity"),
    ("sheet", "sheet_propensity"),
    ("coil", "coil_propensity"),
    ("earlyFolding", "earlyFolding"),
]
FEATURE_COLS = [label for _, label in FEATURES]
GROUP_ORDER = ["C-only", "M-only", "N-only", "C&M", "C&N", "M&N", "C&M&N"]

# HPA Main-location tokens are mapped onto the Layer B biochemical bins.
# Proteins with any Main-location token outside these bins are excluded from the
# clean 7-pattern analysis and counted in QC.
LOCATION_TO_BIN = {
    # Nucleus
    "Nucleoplasm": "N",
    "Nucleoli": "N",
    "Nucleoli fibrillar center": "N",
    "Nucleoli rim": "N",
    "Nuclear bodies": "N",
    "Nuclear speckles": "N",
    "Nuclear membrane": "N",
    "Kinetochore": "N",
    "Mitotic chromosome": "N",
    # Cytosolic/cytoskeletal
    "Cytosol": "C",
    "Actin filaments": "C",
    "Aggresome": "C",
    "Annulus": "C",
    "Basal body": "C",
    "Calyx": "C",
    "Centriolar satellite": "C",
    "Centrosome": "C",
    "Cleavage furrow": "C",
    "Connecting piece": "C",
    "Cytokinetic bridge": "C",
    "Cytoplasmic bodies": "C",
    "End piece": "C",
    "Equatorial segment": "C",
    "Flagellar centriole": "C",
    "Focal adhesion sites": "C",
    "Intermediate filaments": "C",
    "Microtubules": "C",
    "Microtubule ends": "C",
    "Mid piece": "C",
    "Midbody": "C",
    "Midbody ring": "C",
    "Mitotic spindle": "C",
    "Perinuclear theca": "C",
    "Primary cilium": "C",
    "Primary cilium tip": "C",
    "Primary cilium transition zone": "C",
    "Principal piece": "C",
    "Rods & Rings": "C",
    # Membrane/organelle-associated
    "Acrosome": "M",
    "Cell Junctions": "M",
    "Endoplasmic reticulum": "M",
    "Endosomes": "M",
    "Golgi apparatus": "M",
    "Lipid droplets": "M",
    "Lysosomes": "M",
    "Mitochondria": "M",
    "Peroxisomes": "M",
    "Plasma membrane": "M",
    "Vesicles": "M",
}


def log(message=""):
    print(message, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as handle:
        handle.write(str(message) + "\n")


def ensure_outdir():
    os.makedirs(OUT, exist_ok=True)
    with open(LOG_PATH, "w", encoding="utf-8") as handle:
        handle.write("Layer A full run log\n")


def request_json(url, data=None, method=None, timeout=120):
    headers = {"Accept": "application/json", "User-Agent": "ProteinContourLayerA/1.0"}
    body = None
    if data is not None:
        body = urllib.parse.urlencode(data).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8")), response.headers


def download_hpa():
    if not os.path.exists(HPA_ZIP):
        log(f"Downloading HPA 25.1 TSV zip: {HPA_URL}")
        req = urllib.request.Request(HPA_URL, headers={"User-Agent": "ProteinContourLayerA/1.0"})
        with urllib.request.urlopen(req, timeout=180) as response:
            data = response.read()
        with open(HPA_ZIP, "wb") as handle:
            handle.write(data)
    else:
        log(f"Using cached HPA zip: {HPA_ZIP}")

    if not os.path.exists(HPA_TSV):
        with zipfile.ZipFile(HPA_ZIP) as zf:
            members = [name for name in zf.namelist() if name.endswith(".tsv")]
            if len(members) != 1:
                raise RuntimeError(f"Expected one TSV in HPA zip, found {members}")
            with zf.open(members[0]) as src, open(HPA_TSV, "wb") as dst:
                dst.write(src.read())
    log(f"HPA TSV ready: {HPA_TSV}")
    return HPA_TSV


def load_hpa_table():
    hpa_path = download_hpa()
    df = pd.read_csv(hpa_path, sep="\t", dtype=str).fillna("")
    log("\n=== HPA COLUMN CONFIRMATION ===")
    log(", ".join(df.columns.tolist()))
    required = {"Gene", "Gene name", "Reliability", "Main location"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise RuntimeError(f"HPA required columns missing: {missing}")
    return df


def pattern_from_locations(main_location):
    tokens = [token.strip() for token in str(main_location).split(";") if token.strip()]
    if not tokens:
        return None, tokens, ["<empty>"]
    bins = []
    unmapped = []
    for token in tokens:
        mapped = LOCATION_TO_BIN.get(token)
        if mapped is None:
            unmapped.append(token)
        else:
            bins.append(mapped)
    if unmapped or not bins:
        return None, tokens, unmapped
    present = set(bins)
    if present == {"C"}:
        return "C-only", tokens, []
    if present == {"M"}:
        return "M-only", tokens, []
    if present == {"N"}:
        return "N-only", tokens, []
    return "&".join([letter for letter in ["C", "M", "N"] if letter in present]), tokens, []


def build_clean_hpa(df):
    rows = []
    qc = Counter()
    unmapped_tokens = Counter()
    for _, row in df.iterrows():
        if row["Reliability"] == "Uncertain":
            qc["hpa_reliability_uncertain"] += 1
            continue
        pattern, tokens, unmapped = pattern_from_locations(row["Main location"])
        if pattern is None:
            qc["hpa_main_location_not_clean_CMN"] += 1
            for token in unmapped:
                unmapped_tokens[token] += 1
            continue
        qc["hpa_clean_CMN"] += 1
        rows.append(
            {
                "ensembl_gene": row["Gene"],
                "gene_name": row["Gene name"],
                "reliability": row["Reliability"],
                "main_location": row["Main location"],
                "presence_pattern": pattern,
                "main_location_tokens": ";".join(tokens),
            }
        )
    clean = pd.DataFrame(rows)
    return clean, qc, unmapped_tokens


def read_biophys_accessions():
    accessions = set()
    for name in os.listdir(BIO):
        if name.endswith(".json"):
            accessions.add(name[:-5])
    return accessions


def run_uniprot_job(ids, from_db, cache_path, batch_size=400):
    ids = sorted({item for item in ids if item})
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as handle:
            cached = json.load(handle)
        log(f"Using cached UniProt {from_db} mapping: {cache_path}")
        return cached

    all_results = []
    all_failed = []
    for start in range(0, len(ids), batch_size):
        chunk = ids[start : start + batch_size]
        data = {"from": from_db, "to": "UniProtKB", "ids": ",".join(chunk)}
        if from_db == "Gene_Name":
            data["taxId"] = "9606"
        payload, _ = request_json(f"{UNIPROT_IDMAPPING}/run", data=data, method="POST", timeout=120)
        job_id = payload.get("jobId")
        if not job_id:
            raise RuntimeError(f"UniProt mapping did not return jobId for {from_db}: {payload}")

        status_url = f"{UNIPROT_IDMAPPING}/status/{job_id}"
        while True:
            status, _ = request_json(status_url, timeout=120)
            if "jobStatus" not in status:
                break
            if status["jobStatus"] == "FINISHED":
                break
            if status["jobStatus"] in {"FAILED", "ERROR"}:
                raise RuntimeError(f"UniProt mapping failed for job {job_id}: {status}")
            time.sleep(2)

        results_url = f"{UNIPROT_IDMAPPING}/results/{job_id}?format=json&size=500"
        while results_url:
            page, headers = request_json(results_url, timeout=120)
            all_results.extend(page.get("results", []))
            all_failed.extend(page.get("failedIds", []))
            next_link = None
            for part in headers.get("Link", "").split(","):
                if 'rel="next"' in part:
                    next_link = part[part.find("<") + 1 : part.find(">")]
                    break
            results_url = next_link
        log(f"UniProt {from_db}: mapped chunk {start + 1}-{start + len(chunk)} / {len(ids)}")

    out = {"from_db": from_db, "results": all_results, "failedIds": all_failed}
    with open(cache_path, "w", encoding="utf-8") as handle:
        json.dump(out, handle, ensure_ascii=False)
    return out


def entry_accession(entry):
    if isinstance(entry, str):
        return entry
    return entry.get("primaryAccession") or entry.get("uniProtkbId") or ""


def entry_is_reviewed(entry):
    if isinstance(entry, str):
        return False
    return "reviewed" in str(entry.get("entryType", "")).lower() or str(entry.get("reviewed", "")).lower() == "true"


def entry_gene_names(entry):
    if isinstance(entry, str):
        return set()
    names = set()
    for gene in entry.get("genes", []) or []:
        for key in ["geneName", "orderedLocusNames", "orfNames", "synonyms"]:
            value = gene.get(key)
            values = value if isinstance(value, list) else [value]
            for item in values:
                if isinstance(item, dict) and item.get("value"):
                    names.add(str(item["value"]).upper())
    return names


def choose_candidate(candidates, local_accessions, expected_gene=None):
    def score(entry):
        acc = entry_accession(entry)
        reviewed = entry_is_reviewed(entry)
        local = acc in local_accessions
        gene_match = expected_gene and expected_gene.upper() in entry_gene_names(entry)
        return (
            1 if local and reviewed else 0,
            1 if reviewed else 0,
            1 if local else 0,
            1 if gene_match else 0,
            acc,
        )

    return sorted(candidates, key=score, reverse=True)[0]


def build_mapping(mapping_payload, local_accessions, expected_gene_by_from=None):
    grouped = defaultdict(list)
    for item in mapping_payload.get("results", []):
        source = item.get("from")
        if not source:
            continue
        grouped[source].append(item.get("to"))

    rows = []
    chosen = {}
    qc = Counter()
    for source, candidates in grouped.items():
        expected_gene = None
        if expected_gene_by_from:
            expected_gene = expected_gene_by_from.get(source)
        picked = choose_candidate(candidates, local_accessions, expected_gene=expected_gene)
        acc = entry_accession(picked)
        reviewed = entry_is_reviewed(picked)
        in_local = acc in local_accessions
        chosen[source] = acc
        qc["mapped_sources"] += 1
        if len(candidates) > 1:
            qc["multi_hit_sources"] += 1
        if reviewed:
            qc["chosen_reviewed"] += 1
        if in_local:
            qc["chosen_has_biophys_json"] += 1
        rows.append(
            {
                "source_id": source,
                "chosen_uniprot": acc,
                "candidate_count": len(candidates),
                "chosen_reviewed": reviewed,
                "chosen_has_biophys_json": in_local,
                "candidate_accessions": ";".join(sorted({entry_accession(c) for c in candidates if entry_accession(c)})),
            }
        )
    return chosen, pd.DataFrame(rows), qc


def map_hpa_to_uniprot(clean_hpa, local_accessions):
    ensembl_ids = clean_hpa["ensembl_gene"].dropna().astype(str).unique().tolist()
    ensembl_payload = run_uniprot_job(ensembl_ids, "Ensembl", ENSEMBL_MAP_JSON)
    expected_gene_by_ens = clean_hpa.set_index("ensembl_gene")["gene_name"].to_dict()
    ens_map, ens_map_rows, ens_qc = build_mapping(ensembl_payload, local_accessions, expected_gene_by_from=expected_gene_by_ens)
    ens_map_rows.to_csv(os.path.join(OUT, "uniprot_mapping_ensembl_resolved.csv"), index=False)

    mapped_accs = clean_hpa["ensembl_gene"].map(ens_map)
    needs_fallback = clean_hpa.loc[mapped_accs.isna(), "gene_name"].dropna().astype(str).unique().tolist()
    gene_map = {}
    gene_map_rows = pd.DataFrame()
    gene_qc = Counter()
    if needs_fallback:
        gene_payload = run_uniprot_job(needs_fallback, "Gene_Name", GENE_MAP_JSON)
        expected_gene_by_symbol = {name: name for name in needs_fallback}
        gene_map, gene_map_rows, gene_qc = build_mapping(
            gene_payload, local_accessions, expected_gene_by_from=expected_gene_by_symbol
        )
        gene_map_rows.to_csv(os.path.join(OUT, "uniprot_mapping_gene_name_fallback_resolved.csv"), index=False)

    resolved = clean_hpa.copy()
    resolved["uniprot"] = resolved["ensembl_gene"].map(ens_map)
    fallback_mask = resolved["uniprot"].isna()
    resolved.loc[fallback_mask, "uniprot"] = resolved.loc[fallback_mask, "gene_name"].map(gene_map)
    resolved["mapping_source"] = np.where(fallback_mask & resolved["uniprot"].notna(), "Gene_Name_taxId9606", "Ensembl")
    resolved["has_biophys_json"] = resolved["uniprot"].isin(local_accessions)
    resolved.to_csv(os.path.join(OUT, "hpa_clean_CMN_uniprot_resolved.csv"), index=False)
    return resolved, ens_qc, gene_qc, len(needs_fallback)


def fetch_scop3p():
    if os.path.exists(SCOP3P_JSON):
        log(f"Using cached Scop3P bulk JSON: {SCOP3P_JSON}")
        with open(SCOP3P_JSON, encoding="utf-8") as handle:
            return json.load(handle)
    log(f"Fetching Scop3P bulk endpoint: {SCOP3P_URL}")
    payload, _ = request_json(SCOP3P_URL, timeout=240)
    with open(SCOP3P_JSON, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False)
    return payload


def scop3p_modifications(payload):
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        for key in ["modifications", "data", "results"]:
            if isinstance(payload.get(key), list):
                return payload[key]
    raise RuntimeError("Could not find modifications list in Scop3P payload")


def build_phosphosite_index(payload):
    by_acc = defaultdict(dict)
    total = 0
    phospho = 0
    for mod in scop3p_modifications(payload):
        total += 1
        name = str(mod.get("name") or "").strip().lower()
        if name != "phosphorylation":
            continue
        phospho += 1
        acc = mod.get("uniprotId") or mod.get("accession")
        pos = mod.get("position")
        if not acc or pos in (None, ""):
            continue
        try:
            pos = int(pos)
        except ValueError:
            continue
        by_acc[str(acc)][pos] = mod.get("residue")
    return by_acc, {"scop3p_total_modifications": total, "scop3p_phosphorylation_records": phospho}


def read_biophys(acc):
    with open(os.path.join(BIO, acc + ".json"), encoding="utf-8") as handle:
        data = json.load(handle)
    residues = data.get("residues", [])
    by_pos = {int(row["seqpos"]): row for row in residues if "seqpos" in row}
    dis = [float(row["disoMine"]) for row in residues if row.get("disoMine") is not None]
    return by_pos, len(residues), float(np.mean(dis)) if dis else np.nan


def build_site_table(hpa_resolved, phosphosites):
    eligible = hpa_resolved[hpa_resolved["has_biophys_json"] & hpa_resolved["uniprot"].notna()].copy()
    eligible = eligible.drop_duplicates(subset=["uniprot", "presence_pattern"])

    rows = []
    qc = Counter()
    protein_with_sites = set()
    cache = {}
    for _, prot in eligible.iterrows():
        acc = prot["uniprot"]
        sites = phosphosites.get(acc)
        if not sites:
            qc["protein_no_scop3p_phosphosite"] += 1
            continue
        if acc not in cache:
            cache[acc] = read_biophys(acc)
        residues, length, overall_disorder = cache[acc]
        for pos, scop_residue in sorted(sites.items()):
            rr = residues.get(pos)
            if rr is None:
                qc["site_no_biophys_position"] += 1
                continue
            if rr.get("aa") not in {"S", "T", "Y"}:
                qc["site_not_STY_in_biophys"] += 1
                continue
            if any(rr.get(src) is None for src, _ in FEATURES):
                qc["site_missing_feature"] += 1
                continue
            protein_with_sites.add(acc)
            row = {
                "ensembl_gene": prot["ensembl_gene"],
                "gene_name": prot["gene_name"],
                "uniprot": acc,
                "presence_pattern": prot["presence_pattern"],
                "main_location": prot["main_location"],
                "mapping_source": prot["mapping_source"],
                "position": pos,
                "aa": rr.get("aa"),
                "scop3p_residue": scop_residue,
                "protein_length": length,
                "protein_overall_disorder": overall_disorder,
                "log10_protein_length": math.log10(length) if length > 0 else np.nan,
            }
            for src, label in FEATURES:
                row[label] = float(rr[src])
            rows.append(row)
    site_df = pd.DataFrame(rows)
    return eligible, site_df, qc, protein_with_sites


def eps2_kw(h_stat, n, k):
    if n <= k:
        return np.nan
    return float((h_stat - k + 1) / (n - k))


def residualize(y, controls):
    sub = pd.concat([y, controls], axis=1).dropna()
    if len(sub) < controls.shape[1] + 2:
        out = pd.Series(np.nan, index=y.index)
        return out
    yv = sub.iloc[:, 0].astype(float).values
    x = sub.iloc[:, 1:].astype(float).values
    x = np.column_stack([np.ones(len(x)), x])
    beta, _, _, _ = np.linalg.lstsq(x, yv, rcond=None)
    residuals = yv - x @ beta
    out = pd.Series(np.nan, index=y.index)
    out.loc[sub.index] = residuals
    return out


def run_kw_table(site_df, value_cols, prefix):
    rows = []
    for feature in value_cols:
        sub = site_df[["presence_pattern", feature]].dropna()
        samples = [sub.loc[sub["presence_pattern"] == group, feature].astype(float).values for group in GROUP_ORDER]
        valid = [(group, sample) for group, sample in zip(GROUP_ORDER, samples) if len(sample) > 0]
        if len(valid) < 2:
            h_stat, p_value, eps = np.nan, np.nan, np.nan
        else:
            h_stat, p_value = kruskal(*[sample for _, sample in valid])
            eps = eps2_kw(h_stat, sum(len(sample) for _, sample in valid), len(valid))
        medians = {
            f"{prefix}_median_{group}": float(np.median(sample)) if len(sample) else np.nan
            for group, sample in zip(GROUP_ORDER, samples)
        }
        rows.append(
            {
                "feature": feature.replace("_residualized", ""),
                f"{prefix}_H": float(h_stat) if pd.notna(h_stat) else np.nan,
                f"{prefix}_p": float(p_value) if pd.notna(p_value) else np.nan,
                f"{prefix}_epsilon_squared": eps,
                f"{prefix}_n_tested": int(sum(len(sample) for _, sample in valid)),
                f"{prefix}_k_groups": int(len(valid)),
                **medians,
            }
        )
    out = pd.DataFrame(rows)
    valid_p = out[f"{prefix}_p"].notna()
    out[f"{prefix}_p_BH"] = np.nan
    out[f"{prefix}_BH_significant"] = False
    if valid_p.any():
        _, p_bh, _, _ = multipletests(out.loc[valid_p, f"{prefix}_p"], method="fdr_bh")
        out.loc[valid_p, f"{prefix}_p_BH"] = p_bh
        out.loc[valid_p, f"{prefix}_BH_significant"] = p_bh < 0.05
    return out


def run_stats(site_df):
    controls = site_df[["log10_protein_length", "protein_overall_disorder"]]
    controlled = site_df.copy()
    residual_cols = []
    for feature in FEATURE_COLS:
        col = feature + "_residualized"
        controlled[col] = residualize(site_df[feature], controls)
        residual_cols.append(col)

    raw = run_kw_table(site_df, FEATURE_COLS, "raw")
    ctrl = run_kw_table(controlled, residual_cols, "controlled")
    summary = raw.merge(ctrl, on="feature", how="outer")
    ordered_cols = [
        "feature",
        "raw_n_tested",
        "raw_k_groups",
        "raw_H",
        "raw_p",
        "raw_p_BH",
        "raw_BH_significant",
        "raw_epsilon_squared",
        "controlled_n_tested",
        "controlled_k_groups",
        "controlled_H",
        "controlled_p",
        "controlled_p_BH",
        "controlled_BH_significant",
        "controlled_epsilon_squared",
    ]
    return raw, ctrl, summary[ordered_cols], controlled


def write_group_counts(eligible, site_df, protein_with_sites):
    protein_counts = eligible.groupby("presence_pattern")["uniprot"].nunique().reindex(GROUP_ORDER, fill_value=0)
    site_counts = site_df.groupby("presence_pattern").size().reindex(GROUP_ORDER, fill_value=0)
    protein_site_counts = (
        site_df.groupby("presence_pattern")["uniprot"].nunique().reindex(GROUP_ORDER, fill_value=0)
        if len(site_df)
        else pd.Series(0, index=GROUP_ORDER)
    )
    out = pd.DataFrame(
        {
            "presence_pattern": GROUP_ORDER,
            "eligible_hpa_mapped_biophys_proteins": protein_counts.values,
            "proteins_with_scop3p_phosphosites": protein_site_counts.values,
            "site_count": site_counts.values,
        }
    )
    out.to_csv(os.path.join(OUT, "group_counts.csv"), index=False)
    return out


def write_qc_tables(hpa_qc, unmapped_tokens, ens_qc, gene_qc, fallback_n, site_qc, scop_qc):
    qc_rows = []
    for scope, counter in [
        ("hpa", hpa_qc),
        ("uniprot_ensembl", ens_qc),
        ("uniprot_gene_name_fallback", gene_qc),
        ("site_build", site_qc),
        ("scop3p", Counter(scop_qc)),
    ]:
        for key, value in sorted(counter.items()):
            qc_rows.append({"scope": scope, "metric": key, "value": value})
    qc_rows.append({"scope": "uniprot_gene_name_fallback", "metric": "fallback_query_count", "value": fallback_n})
    pd.DataFrame(qc_rows).to_csv(os.path.join(OUT, "qc_counts.csv"), index=False)
    pd.DataFrame(
        [{"token": token, "count": count} for token, count in unmapped_tokens.most_common()]
    ).to_csv(os.path.join(OUT, "hpa_unmapped_main_location_tokens.csv"), index=False)


def write_summary_md(group_counts, stats_summary, site_df):
    raw_max = stats_summary["raw_epsilon_squared"].max()
    ctrl_max = stats_summary["controlled_epsilon_squared"].max()
    top_ctrl = stats_summary.sort_values("controlled_epsilon_squared", ascending=False).iloc[0]
    disorder = stats_summary.loc[stats_summary["feature"] == "disorder_propensity"].iloc[0]
    lines = [
        "# Layer A Full Static Landscape",
        "",
        "Pinned interfaces used:",
        f"- Scop3P bulk: `{SCOP3P_URL}`",
        f"- HPA 25.1 TSV: `{HPA_URL}`",
        "- UniProt idmapping: `Ensembl -> UniProtKB`; fallback `Gene_Name -> UniProtKB` with `taxId=9606`.",
        "",
        "Mapping rule: prefer reviewed UniProtKB entries with local b2bTools JSON; otherwise prefer reviewed, then local JSON, then gene-name match, then sorted accession.",
        "",
        "Clean HPA rule: non-Uncertain rows using `Main location`; every Main-location token must map cleanly to C/M/N. Multiple mapped bins form the 7 presence-pattern groups.",
        "",
        "## Group counts",
        "",
        group_counts.to_markdown(index=False),
        "",
        f"Total analyzed phosphosites: {len(site_df)}",
        "",
        "## Raw and controlled epsilon-squared",
        "",
        stats_summary.to_markdown(index=False),
        "",
        "## Headline",
        "",
        f"Raw max epsilon-squared = {raw_max:.6g}; controlled max epsilon-squared = {ctrl_max:.6g}.",
        f"Top controlled feature: {top_ctrl['feature']} (epsilon-squared={top_ctrl['controlled_epsilon_squared']:.6g}, p_BH={top_ctrl['controlled_p_BH']:.6g}, n={int(top_ctrl['controlled_n_tested'])}).",
        f"Disorder after length + whole-protein disorder control: epsilon-squared={disorder['controlled_epsilon_squared']:.6g}, p_BH={disorder['controlled_p_BH']:.6g}.",
        "",
        "Interpretation is correlational only.",
    ]
    with open(os.path.join(OUT, "README_results.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main():
    ensure_outdir()
    local_accessions = read_biophys_accessions()
    log(f"Local b2bTools JSON proteins: {len(local_accessions)}")

    hpa = load_hpa_table()
    clean_hpa, hpa_qc, unmapped_tokens = build_clean_hpa(hpa)
    clean_hpa.to_csv(os.path.join(OUT, "hpa_clean_CMN_before_uniprot.csv"), index=False)
    log("\n=== HPA CLEANING COUNTS ===")
    for key, value in sorted(hpa_qc.items()):
        log(f"{key}: {value}")
    log("Top unmapped Main-location tokens: " + json.dumps(dict(unmapped_tokens.most_common(12)), ensure_ascii=False))

    resolved, ens_qc, gene_qc, fallback_n = map_hpa_to_uniprot(clean_hpa, local_accessions)
    log("\n=== UNIPROT MAPPING COUNTS ===")
    for label, counter in [("Ensembl", ens_qc), ("Gene_Name fallback", gene_qc)]:
        log(label)
        for key, value in sorted(counter.items()):
            log(f"  {key}: {value}")
    log(f"Gene_Name fallback query count: {fallback_n}")
    log(f"Resolved HPA proteins with local b2bTools JSON: {int(resolved['has_biophys_json'].sum())}")

    scop_payload = fetch_scop3p()
    phosphosites, scop_qc = build_phosphosite_index(scop_payload)
    log("\n=== SCOP3P COUNTS ===")
    for key, value in sorted(scop_qc.items()):
        log(f"{key}: {value}")
    log(f"Proteins with Scop3P phosphorylation records: {len(phosphosites)}")

    eligible, site_df, site_qc, protein_with_sites = build_site_table(resolved, phosphosites)
    site_df.to_csv(os.path.join(OUT, "layerA_full_sites.csv"), index=False)
    log("\n=== SITE BUILD COUNTS ===")
    for key, value in sorted(site_qc.items()):
        log(f"{key}: {value}")
    log(f"Analyzed phosphosites total: {len(site_df)}")
    log(f"Analyzed proteins with phosphosites: {len(protein_with_sites)}")

    group_counts = write_group_counts(eligible, site_df, protein_with_sites)
    log("\n=== GROUP COUNTS (7 presence-pattern groups) ===")
    log(group_counts.to_string(index=False))

    raw, ctrl, summary, controlled_site_df = run_stats(site_df)
    raw.to_csv(os.path.join(OUT, "layerA_full_stats_raw.csv"), index=False)
    ctrl.to_csv(os.path.join(OUT, "layerA_full_stats_controlled.csv"), index=False)
    summary.to_csv(os.path.join(OUT, "layerA_full_stats_summary.csv"), index=False)
    controlled_cols = [
        "gene_name",
        "uniprot",
        "presence_pattern",
        "position",
        "aa",
        "protein_length",
        "protein_overall_disorder",
    ] + FEATURE_COLS + [feature + "_residualized" for feature in FEATURE_COLS]
    controlled_site_df[controlled_cols].to_csv(os.path.join(OUT, "layerA_full_sites_with_residuals.csv"), index=False)

    log("\n=== RAW AND LENGTH+OVERALL-DISORDER CONTROLLED KW / EPSILON-SQUARED ===")
    log(summary.to_string(index=False))

    write_qc_tables(hpa_qc, unmapped_tokens, ens_qc, gene_qc, fallback_n, site_qc, scop_qc)
    write_summary_md(group_counts, summary, site_df)

    log("\nOutputs written to:")
    for name in [
        "group_counts.csv",
        "layerA_full_sites.csv",
        "layerA_full_sites_with_residuals.csv",
        "layerA_full_stats_summary.csv",
        "layerA_full_stats_raw.csv",
        "layerA_full_stats_controlled.csv",
        "qc_counts.csv",
        "README_results.md",
        "run_log.txt",
    ]:
        log(f"  {os.path.join(OUT, name)}")


if __name__ == "__main__":
    main()
