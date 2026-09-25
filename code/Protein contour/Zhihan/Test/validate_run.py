#!/usr/bin/env python3
"""每 run 验收: python3 validate_run.py <report.parquet> <phospho|proteome>

phospho 线判据:
  - report.parquet 可读且非空
  - Modified.Sequence 含 UniMod:21 的行数 > 0
  - 统计 PTM.Site.Confidence >= 0.75 的行数 (occupancy 正式输入的预览, 不作为硬门槛)
  - phosphosites_90/_99 矩阵存在性检查 (报告, 不硬卡)
proteome 线判据:
  - report.parquet 可读且非空
退出码: 0 = PASS, 1 = FAIL, 2 = 环境问题
"""
import sys, os, glob

LOC_CUTOFF = 0.75

def main():
    if len(sys.argv) != 3:
        print("用法: validate_run.py <report.parquet> <phospho|proteome>"); sys.exit(2)
    path, line = sys.argv[1], sys.argv[2]
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        print(f"FAIL: 报告缺失或为空: {path}"); sys.exit(1)
    try:
        import pyarrow.parquet as pq
    except ImportError:
        print("环境问题: 缺 pyarrow (pip install --user pyarrow)"); sys.exit(2)

    pf = pq.ParquetFile(path)
    cols = pf.schema_arrow.names
    seq_col = next((c for c in ("Modified.Sequence", "ModifiedPeptide") if c in cols), None)
    conf_col = next((c for c in cols if "Site.Confidence" in c or "PTM.Site" in c), None)

    total = 0; phos = 0; loc_ok = 0
    read_cols = [c for c in (seq_col, conf_col) if c]
    if not read_cols:
        print(f"FAIL: 未找到序列列。实际列名: {cols}"); sys.exit(1)
    for batch in pf.iter_batches(columns=read_cols, batch_size=100_000):
        total += batch.num_rows
        seqs = batch.column(seq_col).to_pylist() if seq_col else []
        confs = batch.column(conf_col).to_pylist() if conf_col else []
        for i in range(batch.num_rows):
            is_phos = seq_col and seqs[i] and "UniMod:21" in seqs[i]
            if is_phos:
                phos += 1
                if conf_col and confs[i] is not None and confs[i] >= LOC_CUTOFF:
                    loc_ok += 1

    print(f"总前体: {total}")
    if total == 0:
        print("FAIL: 报告为空"); sys.exit(1)

    if line == "phospho":
        print(f"含 UniMod:21: {phos} ({phos/total*100:.1f}%)")
        if conf_col:
            print(f"其中定位置信度 >= {LOC_CUTOFF}: {loc_ok}")
        else:
            print(f"WARN: 未找到定位置信度列 (实际列名: {cols}), 0.75 过滤需人工确认")
        outdir = os.path.dirname(path)
        sites = glob.glob(os.path.join(outdir, "*phosphosites*"))
        print(f"phosphosites 矩阵: {[os.path.basename(s) for s in sites] or '未生成'}")
        if phos == 0:
            print("FAIL: 零磷酸化鉴定"); sys.exit(1)
        print("PASS")
    else:
        if seq_col and phos > 0:
            print(f"WARN: proteome 线出现 {phos} 个 UniMod:21 条目, 库可能用错!")
            sys.exit(1)
        print("PASS")
    sys.exit(0)

if __name__ == "__main__":
    main()
