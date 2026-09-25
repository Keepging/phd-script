#!/bin/bash
# ============================================================
# 01_prepare_inputs.sh — login 节点交互执行, 约 10 分钟
# 产出:
#   $VSC_SCRATCH/rerun/fasta/human_<日期>_UP<release>.fasta
#   $VSC_SCRATCH/rerun/lists/phospho_files.list   (120 个 URL)
#   $VSC_SCRATCH/rerun/lists/proteome_files.list  (120 个 URL)
# ============================================================
set -euo pipefail
BASE="$VSC_SCRATCH/rerun"
mkdir -p "$BASE"/{fasta,lists}

# ---- T2: FASTA (UniProt 最新 release, human reviewed canonical) ----
echo "== 下载 FASTA =="
UURL='https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29'
curl -sf -D "$BASE/fasta/uniprot_headers.txt" -o "$BASE/fasta/human_raw.fasta" "$UURL"
REL=$(grep -i '^x-uniprot-release:' "$BASE/fasta/uniprot_headers.txt" | tr -d '\r' | awk '{print $2}')
[[ -n "$REL" ]] || { echo "FAIL: 未取到 UniProt release 号, 见 uniprot_headers.txt"; exit 1; }
FASTA="$BASE/fasta/human_$(date +%F)_UP${REL}.fasta"
mv "$BASE/fasta/human_raw.fasta" "$FASTA"
N=$(grep -c '>' "$FASTA")
echo "FASTA: $FASTA"
echo "Release: $REL | 蛋白数: $N (预期 ~20,431)"
(( N > 20000 && N < 21500 )) || { echo "FAIL: 蛋白数异常, 检查下载"; exit 1; }
grep -c '|' <(head -1 "$FASTA") >/dev/null || true

# ---- T1: PRIDE 文件清单 ----
echo "== 拉取 PXD023690 文件清单 (分页) =="
python3 - "$BASE/lists" << 'PYEOF'
import json, sys, urllib.request

outdir = sys.argv[1]
acc = "PXD023690"
all_files = []
page = 0
while True:
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files?page={page}&pageSize=100"
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    files = data if isinstance(data, list) else data.get("_embedded", {}).get("files", [])
    if not files:
        break
    all_files.extend(files)
    page += 1
    if page > 30:
        break
print(f"项目文件总数: {len(all_files)}")

def url_of(f):
    for loc in f.get("publicFileLocations", []):
        v = loc.get("value", "")
        if v.startswith(("http", "ftp")):
            return v.replace("ftp://", "https://", 1) if v.startswith("ftp://ftp.pride") else v
    return ""

def pick(tag):
    sel = []
    for f in all_files:
        name = f.get("fileName", "")
        if name.endswith(".raw") and f"DIA_{tag}_HeLa_SUBCELL_EGF" in name:
            u = url_of(f)
            if u:
                sel.append((name, u))
    sel.sort()
    return sel

for tag, outname in (("Phos", "phospho_files.list"), ("Prot", "proteome_files.list")):
    sel = pick(tag)
    with open(f"{outdir}/{outname}", "w") as fh:
        for _, u in sel:
            fh.write(u + "\n")
    print(f"{tag}: {len(sel)} 个文件 -> {outname}")
    if len(sel) != 120:
        print(f"!! 警告: {tag} 数量 != 120 (5 时间点 x 6 组分 x 4 重复)。")
        print("!! 按文档要求: 对不上立即停下回报, 先人工核对命名规则再继续。")
PYEOF

echo "== 完成。核对上方两个数量均为 120 后再进入建库步骤 =="
# 顺手确保 pyarrow 可用 (验证脚本需要)
module load Python/3.13.1-GCCcore-14.2.0 2>/dev/null || true
python3 -c "import pyarrow" 2>/dev/null || pip install --user --quiet pyarrow && echo "pyarrow OK"
