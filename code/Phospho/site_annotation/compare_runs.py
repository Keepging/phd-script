#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""比较两次 annotate_sites.py 输出的逐列命中数与逐位点变化。
用法: python compare_runs.py OLD.tsv NEW.tsv > comparison.tsv"""
import csv
import sys

old_p, new_p = sys.argv[1], sys.argv[2]
rd = lambda p: {r["site_id"]: r for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")}
old, new = rd(old_p), rd(new_p)
assert old.keys() == new.keys(), "两次运行的位点集合不同"
skip = {"site_id", "protein", "gene", "residue", "position", "seq_source", "seq_length", "flank_pm7"}
hit = lambda v: bool(v) and v != "no"
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["column", "n_old", "n_new", "delta", "sites_gained", "sites_lost", "sites_value_changed"])
w.writerow(["residue_check!=OK", sum(r["residue_check"] != "OK" for r in old.values()),
            sum(r["residue_check"] != "OK" for r in new.values()),
            sum(r["residue_check"] != "OK" for r in new.values()) - sum(r["residue_check"] != "OK" for r in old.values()),
            "", "", ""])
for c in new[next(iter(new))].keys():
    if c in skip or c == "residue_check":
        continue
    no = sum(hit(r[c]) for r in old.values())
    nn = sum(hit(r[c]) for r in new.values())
    gained = sum(hit(new[s][c]) and not hit(old[s][c]) for s in new)
    lost = sum(hit(old[s][c]) and not hit(new[s][c]) for s in new)
    changed = sum(hit(old[s][c]) and hit(new[s][c]) and old[s][c] != new[s][c] for s in new)
    w.writerow([c, no, nn, nn - no, gained, lost, changed])
