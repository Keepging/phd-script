#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷酸化位点 ELM 模体 / NLS-NES / 已知结合注释。

用法:
    python annotate_sites.py --site-meta site_meta.tsv --data-dir data --out site_annotation.tsv

输入 site_meta.tsv 至少要有列 site_id、protein、gene；site_id 形如 P04637_S15。
同一 site_id 的多个 multiplicity 行只算一次。

所有数据文件及版本见 README.md；只依赖 Python 标准库 (3.8+)。
"""
import argparse
import collections
import csv
import gzip
import os
import re
import sys

PHOS = set("STY")
WINDOW = 15          # NLS / NES 搜索窗口 (±残基)
FLANK = 7            # 输出的侧翼长度


# ----------------------------------------------------------------------------- 读取
def read_site_meta(path):
    """按 site_id 去重，保留首次出现的行。"""
    sites = collections.OrderedDict()
    with open(path, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            sid = r["site_id"].strip()
            if sid in sites:
                continue
            m = re.match(r"^([A-Z0-9]+(?:-\d+)?)_([A-Z])(\d+)$", sid)
            if not m:
                sys.stderr.write("跳过无法解析的 site_id: %s\n" % sid)
                continue
            sites[sid] = {
                "site_id": sid,
                "protein": r.get("protein", m.group(1)).strip() or m.group(1),
                "gene": (r.get("gene") or "").strip(),
                "residue": m.group(2),
                "position": int(m.group(3)),
            }
    return sites


def read_fasta_gz(path, wanted=None):
    """UniProt FASTA (>sp|ACC|NAME ...) -> {ACC: seq}。wanted 为需要的 accession 集合。"""
    seqs = {}
    opener = gzip.open if path.endswith(".gz") else open
    acc, buf = None, []
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc and (wanted is None or acc in wanted):
                    seqs[acc] = "".join(buf)
                hdr = line[1:].split()[0]
                acc = hdr.split("|")[1] if "|" in hdr else hdr
                buf = []
            else:
                buf.append(line.strip())
        if acc and (wanted is None or acc in wanted):
            seqs[acc] = "".join(buf)
    return seqs


def read_sec_ac(path):
    """UniProt sec_ac.txt: 次级 accession -> 主 accession。"""
    sec = {}
    if not path or not os.path.exists(path):
        return sec
    started = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("Secondary AC"):
                started = True
                continue
            if started:
                p = line.split()
                if len(p) == 2 and not line.startswith("_"):
                    sec[p[0]] = p[1]
    return sec


def read_uniprot_entries(path):
    """protmapper 的 uniprot_entries.tsv.gz: accession -> gene_primary (用于把伙伴 accession 翻译成基因名)。"""
    genes = {}
    if not path or not os.path.exists(path):
        return genes
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r.get("gene_primary"):
                genes[r["accession"]] = r["gene_primary"]
    return genes


def read_elm_classes(motifs_path, classes_tsv_path=None):
    """
    ELM 类正则。优先 SLiMSuite 的 elm2019.motifs (NAME  REGEX  # 描述)，
    再用 elm_classes TSV (ELM 官方格式) 补充缺失的类。
    返回 {ELMIdentifier: {"regex":..., "desc":..., "source":...}}
    """
    classes = {}
    if motifs_path and os.path.exists(motifs_path):
        with open(motifs_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                name, rest = line.split(None, 1)
                rx, _, desc = rest.partition("#")
                classes[name] = {"regex": rx.strip(), "desc": desc.strip(),
                                 "source": os.path.basename(motifs_path)}
    if classes_tsv_path and os.path.exists(classes_tsv_path):
        with open(classes_tsv_path, encoding="utf-8") as fh:
            rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
            for r in rows:
                name = r["ELMIdentifier"]
                if name not in classes:
                    classes[name] = {"regex": r["Regex"], "desc": r.get("Description", ""),
                                     "source": os.path.basename(classes_tsv_path)}
                elif not classes[name]["desc"]:
                    classes[name]["desc"] = r.get("Description", "")
    for name, c in classes.items():
        c["pattern"] = re.compile("(?=(" + c["regex"] + "))")
    return classes


def read_elm_instances(path, organism="Homo sapiens"):
    """ELM instances TSV (官方格式, 带 # 头) -> {acc: [inst,...]}，只保留 true positive 且指定物种。"""
    by_acc = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return by_acc
    with open(path, encoding="utf-8") as fh:
        rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in rows:
            if r.get("InstanceLogic", "true positive") != "true positive":
                continue
            if organism and r.get("Organism") != organism:
                continue
            inst = {
                "id": r["Accession"], "elm": r["ELMIdentifier"], "type": r["ELMType"],
                "primary": r["Primary_Acc"], "start": int(r["Start"]), "end": int(r["End"]),
                "pmid": r.get("References", ""), "methods": r.get("Methods", ""),
            }
            accs = {r["Primary_Acc"].split("-")[0]}
            accs.update(a.split("-")[0] for a in r.get("Accessions", "").split())
            for a in accs:
                by_acc[a].append(inst)
    return by_acc


def read_elm_interactions(path):
    """ELM interactions TSV -> {(elm, acc, start, end): [(partner_acc, pfam, pmid), ...]} (仅人-人)。"""
    out = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "9606" not in r.get("taxonomyElm", "") or "9606" not in r.get("taxonomyDomain", ""):
                continue
            key = (r["Elm"], r["interactorElm"].split("-")[0], int(r["StartElm"]), int(r["StopElm"]))
            out[key].append((r["interactorDomain"], r["Domain"], r.get("PMID", "")))
    return out


def read_elm_interaction_domains(path):
    """ELM class -> [Pfam id:name] (类级别的结合结构域)。"""
    out = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["ELM identifier"]].append("%s:%s" % (r["Interaction Domain Id"], r["Interaction Domain Description"]))
    return out


def read_gmt(path):
    """PTM-POSE 风格 GMT: term \\t source \\t "ACC_S12-Phosphorylation\\t..." -> {ACC_S12: [term,...]}。"""
    out = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 2)
            if len(parts) < 3:
                continue
            term, members = parts[0], parts[2].strip().strip('"')
            for mem in members.split("\t"):
                mem = mem.strip()
                if not mem:
                    continue
                site, _, mod = mem.partition("-")
                if mod and mod != "Phosphorylation":
                    continue
                out[site].append(term)
    return out


def read_psp_sites(path, organism="human"):
    """PhosphoSitePlus Phosphorylation_site_dataset -> {ACC_S12: (LT_LIT, MS_LIT, MS_CST, SITE_GRP_ID)}。"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        # 前几行是版权说明，直到以 GENE 开头的表头
        for line in fh:
            if line.startswith("GENE\t"):
                header = line.rstrip("\n").split("\t")
                break
        else:
            return out
        rd = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        for r in rd:
            if r.get("ORGANISM") != organism:
                continue
            rsd = r["MOD_RSD"].split("-")[0]
            out[r["ACC_ID"] + "_" + rsd] = (r.get("LT_LIT", ""), r.get("MS_LIT", ""), r.get("MS_CST", ""), r.get("SITE_GRP_ID", ""))
    return out


# ----------------------------------------------------------------------------- 匹配
def scan_protein(seq, classes):
    """
    对一条序列扫描所有 ELM 类。返回 list of (elm, start1, end1, phos_positions1)
    start/end 为 1-based 闭区间；phos_positions 为正则中括号捕获的单个 S/T/Y 残基位置 (ELM 用括号标记修饰残基)。
    用 (?=(...)) 前瞻获得重叠匹配。
    """
    hits = []
    for name, c in classes.items():
        for m in c["pattern"].finditer(seq):
            s, e = m.start(1), m.end(1)
            if e <= s:
                continue
            phos = []
            for g in range(2, (m.re.groups or 1) + 1):
                txt = m.group(g)
                if txt is not None and len(txt) == 1 and txt in PHOS:
                    phos.append(m.start(g) + 1)
            hits.append((name, s + 1, e, tuple(phos)))
    return hits


def category(elm):
    if elm.startswith("LIG_SH2_"):
        return "sh2"
    if elm.startswith("LIG_PTB_"):
        return "ptb"
    if elm.startswith("LIG_14-3-3_"):
        return "1433"
    if elm.startswith("LIG_"):
        return "other_lig"
    return "other_type"


def fmt_hit(elm, s, e):
    return "%s[%d-%d]" % (elm, s, e)


def join(items):
    return ";".join(items)


# ----------------------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-meta", default="site_meta.tsv")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--fasta", default=None, help="覆盖默认的 Swiss-Prot FASTA (可用最新 UniProt 导出)")
    ap.add_argument("--elm-motifs", default=None, help="SLiMSuite .motifs 格式的 ELM 类 (默认 data/elm/elm2019.motifs)")
    ap.add_argument("--elm-classes", default=None, help="ELM 官方 elms_index.tsv / elm_classes.tsv (补充或替换)")
    ap.add_argument("--elm-instances", default=None)
    ap.add_argument("--elm-interactions", default=None)
    ap.add_argument("--elm-interaction-domains", default=None)
    ap.add_argument("--out", default="site_annotation.tsv")
    args = ap.parse_args()

    D = args.data_dir
    p = lambda *x: os.path.join(D, *x)
    fasta = args.fasta or p("uniprot", "uniprot_sprot.fasta.gz")
    # 传空字符串 '' 表示不使用该文件 (例如只用官方 elms_index.tsv 时 --elm-motifs '')
    dflt = lambda v, *x: p(*x) if v is None else v
    elm_motifs = dflt(args.elm_motifs, "elm", "elm2019.motifs")
    elm_classes_tsv = dflt(args.elm_classes, "elm", "elm_classes.2018-07-17.tsv")
    elm_instances = dflt(args.elm_instances, "elm", "elm_instances.2018-07-17.tsv")
    elm_interactions = dflt(args.elm_interactions, "elm", "elm_interactions.2018.tsv")
    elm_idom = dflt(args.elm_interaction_domains, "elm", "elm_interaction_domains.2018.tsv")

    sites = read_site_meta(args.site_meta)
    accs = {s["protein"] for s in sites.values()}
    sys.stderr.write("位点 %d (去重后), 蛋白 %d\n" % (len(sites), len(accs)))

    sec = read_sec_ac(p("uniprot", "uniprot_sec_ac.txt.gz"))
    wanted = set(accs) | {sec[a] for a in accs if a in sec}
    seqs = read_fasta_gz(fasta, wanted)
    genes = read_uniprot_entries(p("uniprot", "uniprot_entries.tsv.gz"))
    sys.stderr.write("找到序列 %d / %d\n" % (sum(1 for a in accs if a in seqs or sec.get(a) in seqs), len(accs)))

    classes = read_elm_classes(elm_motifs, elm_classes_tsv)
    instances = read_elm_instances(elm_instances)
    interactions = read_elm_interactions(elm_interactions)
    class_domains = read_elm_interaction_domains(elm_idom)
    sys.stderr.write("ELM 类 %d, 人源实例蛋白 %d\n" % (len(classes), len(instances)))

    A = lambda n: p("ptm_pose_annotations", n)
    psp_inter = read_gmt(A("PhosphoSitePlus_Interactions.gmt"))
    psp_func = read_gmt(A("PhosphoSitePlus_Function.gmt"))
    psp_proc = read_gmt(A("PhosphoSitePlus_Process.gmt"))
    psp_enz = read_gmt(A("PhosphoSitePlus_Enzyme.gmt"))
    ptmint = read_gmt(A("PTMInt_Interactions.gmt"))
    omni_enz = read_gmt(A("OmniPath_Writer_Enzyme.gmt"))
    psp_sites = read_psp_sites(p("psp", "Phosphorylation_site_dataset.tsv.gz"))

    def gene_of(acc):
        g = genes.get(acc.split("-")[0])
        return "%s(%s)" % (g, acc) if g else acc

    # 按蛋白扫描一次正则
    scan_cache = {}

    cols = [
        "site_id", "protein", "gene", "residue", "position",
        "seq_source", "seq_length", "residue_check", "flank_pm7",
        "elm_sh2_instance", "elm_sh2_regex",
        "elm_ptb_instance", "elm_ptb_regex",
        "elm_1433_instance", "elm_1433_regex",
        "elm_other_lig_instance", "elm_other_lig_regex", "elm_other_lig_regex_span",
        "elm_other_type_instance", "elm_other_type_regex",
        "elm_instance_detail", "elm_instance_partner", "elm_class_domain",
        "nls_15", "nls_15_instance", "nls_15_regex",
        "nes_15", "nes_15_instance", "nes_15_regex",
        "psp_binding", "ptmint_binding", "psp_function", "psp_process",
        "kinase_psp", "enzyme_omnipath_writer",
        "psp_site_known", "psp_lit_LT/MS/CST",
    ]

    n_ok = n_mm = n_noseq = 0
    with open(args.out, "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for sid, s in sites.items():
            acc, res, pos = s["protein"], s["residue"], s["position"]
            seq_acc = acc if acc in seqs else sec.get(acc) if sec.get(acc) in seqs else None
            seq = seqs.get(seq_acc) if seq_acc else None
            row = {c: "" for c in cols}
            row.update({"site_id": sid, "protein": acc, "gene": s["gene"], "residue": res, "position": pos})

            if seq is None:
                row["seq_source"] = ""
                row["residue_check"] = "NO_SEQ"
                n_noseq += 1
            else:
                row["seq_source"] = "swissprot:" + seq_acc + ("" if seq_acc == acc else "(secondary->primary)")
                row["seq_length"] = len(seq)
                if pos > len(seq):
                    row["residue_check"] = "OUT_OF_RANGE"
                    n_mm += 1
                elif seq[pos - 1] != res:
                    row["residue_check"] = "MISMATCH:" + seq[pos - 1]
                    n_mm += 1
                else:
                    row["residue_check"] = "OK"
                    n_ok += 1
                    lo, hi = max(0, pos - 1 - FLANK), min(len(seq), pos + FLANK)
                    fl = seq[lo:pos - 1].rjust(FLANK, "_") + seq[pos - 1].lower() + seq[pos:hi].ljust(FLANK, "_")
                    row["flank_pm7"] = fl

            if row["residue_check"] == "OK":
                # ---------- 正则
                if seq_acc not in scan_cache:
                    scan_cache[seq_acc] = scan_protein(seq, classes)
                hits = scan_cache[seq_acc]
                cat_pres = collections.defaultdict(list)   # 位点正是被捕获的磷酸残基
                cat_span = collections.defaultdict(list)   # 位点仅落在匹配区间内
                nls_rx, nes_rx = [], []
                for elm, hs, he, phos in hits:
                    if elm.startswith("TRG_NLS_") and hs <= pos + WINDOW and he >= pos - WINDOW:
                        nls_rx.append(fmt_hit(elm, hs, he))
                    if elm.startswith("TRG_NES_") and hs <= pos + WINDOW and he >= pos - WINDOW:
                        nes_rx.append(fmt_hit(elm, hs, he))
                    if hs <= pos <= he:
                        cat = category(elm)
                        if pos in phos:
                            cat_pres[cat].append(fmt_hit(elm, hs, he))
                        else:
                            cat_span[cat].append(fmt_hit(elm, hs, he))
                row["elm_sh2_regex"] = join(sorted(set(cat_pres["sh2"])))
                row["elm_ptb_regex"] = join(sorted(set(cat_pres["ptb"])))
                row["elm_1433_regex"] = join(sorted(set(cat_pres["1433"])))
                row["elm_other_lig_regex"] = join(sorted(set(cat_pres["other_lig"])))
                row["elm_other_lig_regex_span"] = join(sorted(set(cat_span["other_lig"])))
                row["elm_other_type_regex"] = join(sorted(set(cat_pres["other_type"])))
                row["nls_15_regex"] = join(sorted(set(nls_rx)))
                row["nes_15_regex"] = join(sorted(set(nes_rx)))

                # ---------- 实验实例
                inst_cat = collections.defaultdict(list)
                inst_detail, partners, cdoms = [], [], []
                nls_in, nes_in = [], []
                for inst in instances.get(acc, []):
                    elm = inst["elm"]
                    overlap_win = inst["start"] <= pos + WINDOW and inst["end"] >= pos - WINDOW
                    if elm.startswith("TRG_NLS_") and overlap_win:
                        nls_in.append(fmt_hit(elm, inst["start"], inst["end"]))
                    if elm.startswith("TRG_NES_") and overlap_win:
                        nes_in.append(fmt_hit(elm, inst["start"], inst["end"]))
                    if not (inst["start"] <= pos <= inst["end"]):
                        continue
                    # 实例坐标来自 ELM 当时的序列版本；检查该区间在当前序列上是否仍匹配类正则
                    c = classes.get(elm)
                    consistent = "?"
                    if c:
                        sub = seq[inst["start"] - 1:inst["end"]]
                        consistent = "ok" if re.fullmatch(c["regex"], sub) else "span_regex_mismatch"
                    iso = "" if "-" not in inst["primary"] else ",isoform:" + inst["primary"]
                    tag = "%s[%d-%d]" % (elm, inst["start"], inst["end"])
                    inst_cat[category(elm)].append(tag)
                    inst_detail.append("%s(%s,PMID:%s,%s%s)" % (tag, inst["id"], inst["pmid"], consistent, iso))
                    for part_acc, pfam, pmid in interactions.get((elm, acc, inst["start"], inst["end"]), []):
                        partners.append("%s->%s[%s]" % (elm, gene_of(part_acc), pfam))
                    for d in class_domains.get(elm, []):
                        cdoms.append("%s->%s" % (elm, d))
                row["elm_sh2_instance"] = join(sorted(set(inst_cat["sh2"])))
                row["elm_ptb_instance"] = join(sorted(set(inst_cat["ptb"])))
                row["elm_1433_instance"] = join(sorted(set(inst_cat["1433"])))
                row["elm_other_lig_instance"] = join(sorted(set(inst_cat["other_lig"])))
                row["elm_other_type_instance"] = join(sorted(set(inst_cat["other_type"])))
                row["elm_instance_detail"] = join(inst_detail)
                row["elm_instance_partner"] = join(sorted(set(partners)))
                row["elm_class_domain"] = join(sorted(set(cdoms)))
                row["nls_15_instance"] = join(sorted(set(nls_in)))
                row["nes_15_instance"] = join(sorted(set(nes_in)))
                row["nls_15"] = "yes" if (nls_rx or nls_in) else "no"
                row["nes_15"] = "yes" if (nes_rx or nes_in) else "no"

            # ---------- 已知结合 / 功能 / 激酶 (与序列版本无关，按 accession+位点查)
            key = "%s_%s%d" % (acc, res, pos)
            row["psp_binding"] = join(sorted(set(psp_inter.get(key, []))))
            row["ptmint_binding"] = join(sorted(set(
                "%s->%s" % (gene_of(t.split("->")[0]), t.split("->")[1]) if "->" in t else t
                for t in ptmint.get(key, []))))
            row["psp_function"] = join(sorted(set(psp_func.get(key, []))))
            row["psp_process"] = join(sorted(set(psp_proc.get(key, []))))
            row["kinase_psp"] = join(sorted(set(psp_enz.get(key, []))))
            row["enzyme_omnipath_writer"] = join(sorted(set(omni_enz.get(key, []))))
            if key in psp_sites:
                lt, ms, cst, _ = psp_sites[key]
                row["psp_site_known"] = "yes"
                row["psp_lit_LT/MS/CST"] = "%s/%s/%s" % (lt or "0", ms or "0", cst or "0")
            else:
                row["psp_site_known"] = "no"

            w.writerow([row[c] for c in cols])

    sys.stderr.write("残基核对: OK %d, 不一致/越界 %d, 无序列 %d\n" % (n_ok, n_mm, n_noseq))
    sys.stderr.write("输出 %s\n" % args.out)


if __name__ == "__main__":
    main()
