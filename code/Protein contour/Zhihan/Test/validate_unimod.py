#!/usr/bin/env python3
"""检查 DIA-NN 输出 (谱库 parquet / 主报告 parquet / 矩阵 tsv) 中是否含 UniMod:21 条目。

用法: python3 validate_unimod.py <文件路径>
支持: .parquet (需要 pyarrow 或 pandas+pyarrow) 和 .tsv
退出码: 0 = PASS, 1 = FAIL, 2 = 环境/输入问题
"""
import sys
import os

TARGET = "UniMod:21"


def fail_env(msg):
    print(f"环境问题: {msg}")
    sys.exit(2)


def report(total, hits, label):
    print(f"文件: {label}")
    print(f"总行数: {total}")
    print(f"含 {TARGET} 的行数: {hits}")
    if hits > 0:
        pct = hits / total * 100 if total else 0.0
        print(f"占比: {pct:.1f}%")
        if pct < 5.0:
            print(f"WARN: {TARGET} 占比异常低 (<5%), 建库配置可能仍有问题, 建议人工抽查")
        print(f"PASS — 检测到 {hits} 个含 {TARGET} 的条目")
        sys.exit(0)
    else:
        print(f"FAIL — 未检测到任何 {TARGET} 条目")
        print("可能原因:")
        print(f"  1. --var-mod '{TARGET},79.966331,STY' 未在建库命令中生效")
        print("  2. (仅搜索报告) 输入的 .raw 文件不含目标修饰的信号")
        print("  3. 序列列名与预期不同 — 见上方打印的实际列名")
        sys.exit(1)


def candidate_columns(cols):
    """优先精确匹配, 其次任何含 modified/sequence 的列。"""
    exact = [c for c in cols if c in ("Modified.Sequence", "ModifiedPeptide", "ModifiedSequence")]
    if exact:
        return exact
    return [c for c in cols if "modified" in c.lower() or "sequence" in c.lower()]


def check_parquet(path):
    try:
        import pyarrow.parquet as pq
    except ImportError:
        fail_env("缺少 pyarrow。请执行: pip install --user pyarrow  (或改用 tsv 矩阵作为验证对象)")
    pf = pq.ParquetFile(path)
    cols = pf.schema_arrow.names
    print(f"实际列名 ({len(cols)}): {cols}")
    targets = candidate_columns(cols)
    if not targets:
        print("FAIL: 未找到任何序列相关列, 无法验证")
        sys.exit(1)
    print(f"检查列: {targets}")
    total = 0
    hits = 0
    for batch in pf.iter_batches(columns=targets, batch_size=100_000):
        n = batch.num_rows
        total += n
        row_hit = [False] * n
        for col in targets:
            vals = batch.column(col).to_pylist()
            for i, v in enumerate(vals):
                if v and TARGET in str(v):
                    row_hit[i] = True
        hits += sum(row_hit)
    report(total, hits, path)


def check_tsv(path):
    import csv
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames or []
        print(f"实际列名 ({len(cols)}): {cols}")
        targets = candidate_columns(cols)
        if not targets:
            print("FAIL: 未找到任何序列相关列, 无法验证")
            sys.exit(1)
        print(f"检查列: {targets}")
        total = 0
        hits = 0
        for row in reader:
            total += 1
            if any(TARGET in (row.get(c) or "") for c in targets):
                hits += 1
    report(total, hits, path)


def main():
    if len(sys.argv) != 2:
        fail_env("用法: python3 validate_unimod.py <report.parquet | lib.parquet | *.tsv>")
    path = sys.argv[1]
    if not os.path.isfile(path):
        fail_env(f"文件不存在: {path}")
    if os.path.getsize(path) == 0:
        print(f"FAIL: 文件为空: {path}")
        sys.exit(1)
    if path.endswith(".parquet"):
        check_parquet(path)
    elif path.endswith((".tsv", ".txt", ".csv")):
        check_tsv(path)
    else:
        fail_env(f"不认识的扩展名: {path} (支持 .parquet / .tsv)")


if __name__ == "__main__":
    main()
