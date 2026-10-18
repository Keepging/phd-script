#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
导出 site_meta.tsv 中全部 accession 的 FASTA，按批次拆分，供 14-3-3-Pred / GPS-PBS 等网页服务器手动提交。

用法:
    python make_fasta_batches.py --site-meta site_meta.tsv --fasta data/uniprot/uniprot_sprot.fasta.gz \
        --out-dir fasta_batches --batch-size 50

输出:
    fasta_batches/all_accessions.fasta   全部序列
    fasta_batches/batch_NNN.fasta        每批 --batch-size 条
    fasta_batches/batch_index.tsv        accession -> 批次文件、序列长度、是否存在残基不一致
    --batch-size 1 时每个 accession 单独一个文件 (fasta_batches/single/ACC.fasta)
"""
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_sites import read_fasta_gz, read_sec_ac, read_site_meta


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-meta", default="site_meta.tsv")
    ap.add_argument("--fasta", default="data/uniprot/uniprot_sprot.fasta.gz")
    ap.add_argument("--sec-ac", default="data/uniprot/uniprot_sec_ac.txt.gz")
    ap.add_argument("--annotation", default="site_annotation.tsv",
                    help="annotate_sites.py 的输出；用于在索引里标记残基不一致的蛋白 (可选)")
    ap.add_argument("--out-dir", default="fasta_batches")
    ap.add_argument("--batch-size", type=int, default=50)
    args = ap.parse_args()

    sites = read_site_meta(args.site_meta)
    acc_gene = {}
    for s in sites.values():
        acc_gene.setdefault(s["protein"], s["gene"])
    accs = sorted(acc_gene)
    sec = read_sec_ac(args.sec_ac)
    seqs = read_fasta_gz(args.fasta, set(accs) | {sec[a] for a in accs if a in sec})

    suspect = {}
    if args.annotation and os.path.exists(args.annotation):
        with open(args.annotation, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("residue_check", "OK") != "OK":
                    suspect.setdefault(r["protein"], []).append(r["site_id"] + ":" + r["residue_check"])

    os.makedirs(args.out_dir, exist_ok=True)
    for old in os.listdir(args.out_dir):
        if re.match(r"^(batch_\d+|all_accessions)\.fasta$", old) or old == "batch_index.tsv":
            os.remove(os.path.join(args.out_dir, old))

    records, missing = [], []
    for a in accs:
        src = a if a in seqs else sec.get(a)
        if src not in seqs:
            missing.append(a)
            continue
        records.append((a, acc_gene[a], seqs[src]))

    with open(os.path.join(args.out_dir, "all_accessions.fasta"), "w", encoding="utf-8") as fo:
        for a, g, s in records:
            fo.write(">%s %s len=%d\n%s\n" % (a, g, len(s), wrap(s)))

    index_rows = []
    if args.batch_size <= 1:
        sdir = os.path.join(args.out_dir, "single")
        os.makedirs(sdir, exist_ok=True)
        for a, g, s in records:
            fn = os.path.join("single", a + ".fasta")
            with open(os.path.join(args.out_dir, fn), "w", encoding="utf-8") as fo:
                fo.write(">%s %s len=%d\n%s\n" % (a, g, len(s), wrap(s)))
            index_rows.append((a, g, len(s), fn))
    else:
        nb = (len(records) + args.batch_size - 1) // args.batch_size
        for i in range(nb):
            fn = "batch_%03d.fasta" % (i + 1)
            chunk = records[i * args.batch_size:(i + 1) * args.batch_size]
            with open(os.path.join(args.out_dir, fn), "w", encoding="utf-8") as fo:
                for a, g, s in chunk:
                    fo.write(">%s %s len=%d\n%s\n" % (a, g, len(s), wrap(s)))
                    index_rows.append((a, g, len(s), fn))

    with open(os.path.join(args.out_dir, "batch_index.tsv"), "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "gene", "length", "batch_file", "n_sites", "residue_mismatch_sites"])
        nsites = {}
        for s in sites.values():
            nsites[s["protein"]] = nsites.get(s["protein"], 0) + 1
        for a, g, L, fn in index_rows:
            w.writerow([a, g, L, fn, nsites.get(a, 0), ";".join(suspect.get(a, []))])
        for a in missing:
            w.writerow([a, acc_gene[a], "", "NO_SEQUENCE", nsites.get(a, 0), ""])

    print("accession %d, 有序列 %d, 无序列 %d, 批次文件 %d, 每批 %d" % (
        len(accs), len(records), len(missing),
        len({r[3] for r in index_rows}), args.batch_size))
    if suspect:
        print("残基不一致的蛋白 %d 个 (序列版本可能已更新，见 batch_index.tsv 最后一列)" % len(suspect))


if __name__ == "__main__":
    main()
