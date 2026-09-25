import csv
import json
import os
import urllib.request
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests


ROOT = r"D:\博士\Protein contour\Phospho"
OUT = os.path.join(ROOT, "pilot_layerA_codex")
HPA = os.path.join(ROOT, "pilot_layerA", "subcellular_location.tsv")
BIO = os.path.join(ROOT, "biophys_json")
GENE_MAP = os.path.join(ROOT, "gene_uniprot_map.csv")
SCOP3P_CACHE = os.path.join(OUT, "scop3p_all_modifications.json")
SITES_OUT = os.path.join(OUT, "layerA_pilot_sites.csv")
STATS_OUT = os.path.join(OUT, "layerA_pilot_stats.csv")

FEATURES = [
    "backbone",
    "sidechain",
    "disoMine",
    "helix",
    "sheet",
    "coil",
    "earlyFolding",
]

# Clean single-compartment HPA main locations, aligned to Layer B's 3 bins.
COMPARTMENT_TOKENS = {
    "Nucleus": {
        "Nucleoplasm",
        "Nucleoli",
        "Nucleoli fibrillar center",
        "Nucleoli rim",
        "Nuclear bodies",
        "Nuclear speckles",
        "Nuclear membrane",
        "Kinetochore",
        "Mitotic chromosome",
    },
    "Cytosol": {"Cytosol"},
    "Membrane": {"Plasma membrane", "Endoplasmic reticulum"},
}


def read_gene_map():
    gene_to_uniprot = {}
    with open(GENE_MAP, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            gene = row.get("gene")
            acc = row.get("uniprot")
            if gene and acc and gene not in gene_to_uniprot:
                gene_to_uniprot[gene] = acc
    return gene_to_uniprot


def load_hpa_assignments(gene_to_uniprot):
    assignments = {}
    stage = Counter()
    with open(HPA, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Reliability"] == "Uncertain":
                stage["hpa_uncertain"] += 1
                continue
            main = [token.strip() for token in row["Main location"].split(";") if token.strip()]
            if not main:
                stage["hpa_no_main"] += 1
                continue

            comp = None
            for candidate, tokens in COMPARTMENT_TOKENS.items():
                if set(main) <= tokens:
                    comp = candidate
                    break
            if comp is None:
                stage["hpa_not_clean_3comp"] += 1
                continue

            stage["hpa_clean_3comp"] += 1
            acc = gene_to_uniprot.get(row["Gene name"])
            if not acc:
                stage["no_uniprot_map"] += 1
                continue
            if not os.path.exists(os.path.join(BIO, acc + ".json")):
                stage["no_biophys_json"] += 1
                continue
            if acc in assignments and assignments[acc] != comp:
                stage["conflicting_assignment"] += 1
                continue
            assignments[acc] = comp
    return assignments, stage


def fetch_scop3p_all():
    if os.path.exists(SCOP3P_CACHE):
        with open(SCOP3P_CACHE, encoding="utf-8") as handle:
            return json.load(handle), "cache"

    url = "https://iomics.ugent.be/scop3p/api/get-all-modifications"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as response:
        payload = json.loads(response.read().decode("utf-8"))
    with open(SCOP3P_CACHE, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    return payload, "live"


def build_phosphosite_index(payload):
    by_acc = defaultdict(dict)
    for mod in payload.get("modifications", []):
        name = str(mod.get("name") or "").lower()
        if "phospho" not in name:
            continue
        acc = mod.get("uniprotId")
        pos = mod.get("position")
        residue = mod.get("residue")
        if not acc or pos is None:
            continue
        by_acc[acc][int(pos)] = residue
    return by_acc


def read_residue_features(acc):
    with open(os.path.join(BIO, acc + ".json"), encoding="utf-8") as handle:
        data = json.load(handle)
    return {int(row["seqpos"]): row for row in data["residues"]}


def build_site_table(assignments, phosphosites):
    rows = []
    dropped = Counter()
    for acc, comp in sorted(assignments.items()):
        sites = phosphosites.get(acc)
        if not sites:
            dropped["protein_no_scop3p_phosphosite"] += 1
            continue
        residues = read_residue_features(acc)
        for pos, scop_residue in sorted(sites.items()):
            rr = residues.get(pos)
            if rr is None:
                dropped["site_no_biophys_position"] += 1
                continue
            if rr.get("aa") not in {"S", "T", "Y"}:
                dropped["site_not_STY_in_biophys"] += 1
                continue
            if any(rr.get(feature) is None for feature in FEATURES):
                dropped["site_missing_feature"] += 1
                continue
            row = {
                "uniprot": acc,
                "compartment": comp,
                "position": pos,
                "aa": rr["aa"],
                "scop3p_residue": scop_residue,
            }
            for feature in FEATURES:
                row[feature] = float(rr[feature])
            rows.append(row)
    return rows, dropped


def eps2_kw(samples, h_stat):
    k = len(samples)
    n = sum(len(sample) for sample in samples)
    return (h_stat - k + 1) / (n - k) if n > k else float("nan")


def run_stats(rows):
    groups = ["Nucleus", "Cytosol", "Membrane"]
    stats = []
    for feature in FEATURES:
        samples = [[row[feature] for row in rows if row["compartment"] == group] for group in groups]
        h_stat, p_value = kruskal(*samples)
        medians = {
            "median_" + group: float(np.median(sample))
            for group, sample in zip(groups, samples)
        }
        stats.append(
            {
                "feature": feature,
                "H": float(h_stat),
                "p": float(p_value),
                "epsilon_squared": float(eps2_kw(samples, h_stat)),
                "n_total": sum(len(sample) for sample in samples),
                **medians,
            }
        )

    _, p_bh, _, _ = multipletests([row["p"] for row in stats], method="fdr_bh")
    for row, adjusted in zip(stats, p_bh):
        row["p_BH"] = float(adjusted)
        row["BH_significant"] = bool(adjusted < 0.05)
    return stats


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    os.makedirs(OUT, exist_ok=True)

    gene_to_uniprot = read_gene_map()
    assignments, hpa_stage = load_hpa_assignments(gene_to_uniprot)
    by_comp_proteins = Counter(assignments.values())

    payload, scop3p_source = fetch_scop3p_all()
    phosphosites = build_phosphosite_index(payload)
    site_rows, dropped = build_site_table(assignments, phosphosites)
    stats = run_stats(site_rows)

    write_csv(
        SITES_OUT,
        site_rows,
        ["uniprot", "compartment", "position", "aa", "scop3p_residue"] + FEATURES,
    )
    write_csv(
        STATS_OUT,
        stats,
        [
            "feature",
            "H",
            "p",
            "p_BH",
            "BH_significant",
            "epsilon_squared",
            "n_total",
            "median_Nucleus",
            "median_Cytosol",
            "median_Membrane",
        ],
    )

    print("=== VERIFIED INPUTS ===")
    print("Scop3P endpoint: https://iomics.ugent.be/scop3p/api/get-all-modifications")
    print(f"Scop3P payload source: {scop3p_source}; total modifications: {len(payload.get('modifications', []))}")
    print(f"HPA file: {HPA}")
    print(f"b2bTools JSON directory: {BIO}")

    print("\n=== HPA / MAPPING INTERMEDIATES ===")
    for key in sorted(hpa_stage):
        print(f"{key}: {hpa_stage[key]}")
    for comp in ["Nucleus", "Cytosol", "Membrane"]:
        print(f"mapped proteins {comp}: {by_comp_proteins[comp]}")

    print("\n=== SITE COUNTS PER COMPARTMENT ===")
    site_counts = Counter(row["compartment"] for row in site_rows)
    protein_with_sites = defaultdict(set)
    for row in site_rows:
        protein_with_sites[row["compartment"]].add(row["uniprot"])
    for comp in ["Nucleus", "Cytosol", "Membrane"]:
        print(f"{comp}: {site_counts[comp]} sites from {len(protein_with_sites[comp])} proteins")
    print(f"total mapped sites: {len(site_rows)}")
    for key in sorted(dropped):
        print(f"dropped {key}: {dropped[key]}")

    print("\n=== KRUSKAL-WALLIS + BH + EPSILON-SQUARED ===")
    print("feature\tH\tp\tp_BH\tBH_significant\tepsilon_squared\tmedian_Nucleus\tmedian_Cytosol\tmedian_Membrane")
    for row in sorted(stats, key=lambda item: item["epsilon_squared"], reverse=True):
        print(
            f"{row['feature']}\t{row['H']:.6g}\t{row['p']:.6g}\t{row['p_BH']:.6g}\t"
            f"{row['BH_significant']}\t{row['epsilon_squared']:.6g}\t"
            f"{row['median_Nucleus']:.6g}\t{row['median_Cytosol']:.6g}\t{row['median_Membrane']:.6g}"
        )
    print(f"\nHEADLINE max epsilon_squared: {max(row['epsilon_squared'] for row in stats):.6g}")
    print(f"Outputs: {SITES_OUT}; {STATS_OUT}; {SCOP3P_CACHE}")


if __name__ == "__main__":
    main()
