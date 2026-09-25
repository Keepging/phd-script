#!/bin/bash
# ============================================================
# 列出 PXD023690 (Martinez-Val 2021) 的可下载 .raw 文件
# 用 PRIDE REST API, 不猜 FTP 年月路径
# 交互执行, 从输出中人工挑 1 个 phospho-enriched 的 run,
# 把其 URL 作为参数传给 stage2: sbatch stage2_full_search.slurm <RAW_URL>
# ============================================================
set -euo pipefail
ACC="${1:-PXD023690}"
API="https://www.ebi.ac.uk/pride/ws/archive/v2/projects/${ACC}/files?pageSize=500"

curl -sf "$API" | python3 -c '
import json, sys
data = json.load(sys.stdin)
files = data if isinstance(data, list) else data.get("_embedded", {}).get("files", [])
rows = []
for f in files:
    name = f.get("fileName", "")
    if not name.lower().endswith(".raw"):
        continue
    size_mb = (f.get("fileSizeBytes") or 0) / 1e6
    url = ""
    for loc in f.get("publicFileLocations", []):
        if "ftp" in (loc.get("name") or "").lower() or (loc.get("value") or "").startswith(("ftp", "http")):
            url = loc.get("value", "")
            break
    rows.append((size_mb, name, url))
rows.sort()
print(f"共 {len(rows)} 个 .raw 文件 (按大小升序):\n")
for size_mb, name, url in rows:
    tag = " <-- 疑似 phospho" if any(k in name.lower() for k in ("phos", "psty", "enrich")) else ""
    print(f"{size_mb:9.1f} MB  {name}{tag}")
    print(f"             {url}")
'
