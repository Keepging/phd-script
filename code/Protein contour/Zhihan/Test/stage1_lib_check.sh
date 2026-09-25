#!/bin/bash
# ============================================================
# 阶段 1: 小 FASTA 建库 + 验证 UniMod:21 是否进库
# 目的: 隔离验证 "--var-mod UniMod:21 修复是否生效" (纯建库端问题,
#       与样品表达谱无关, 因此小 FASTA 是合法且最快的验证方式)
# 运行: login 节点交互执行即可, 预期 < 30 min
# ============================================================
set -euo pipefail

# module 在非 login shell 中可能不可用, 兜底
source /etc/profile.d/modules.sh 2>/dev/null || true

# ---- 配置区 ----
DIANN="$VSC_DATA/software/diann/diann-2.2.0/diann-linux"
WORKDIR="$VSC_SCRATCH/phospho_test"
FASTA_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
N_PROTEINS="${N_PROTEINS:-150}"   # 子集蛋白数, 可用环境变量覆盖
THREADS="${THREADS:-$(nproc)}"

module load GCCcore/14.2.0
mkdir -p "$WORKDIR"/{input,output}
cd "$WORKDIR"

# ---- S1: 下载完整 FASTA (阶段 2 复用, 只下一次) ----
if [[ ! -s input/human.fasta ]]; then
    echo "[S1] 下载 FASTA ..."
    wget -q -O input/human.fasta.gz "$FASTA_URL"
    gunzip -f input/human.fasta.gz
fi
[[ -s input/human.fasta ]] || { echo "FAIL: FASTA 为空"; exit 1; }
head -c 1 input/human.fasta | grep -q '>' || { echo "FAIL: FASTA 首字符非 '>'"; exit 1; }
echo "[S1] OK: $(grep -c '>' input/human.fasta) 条蛋白"

# ---- S1b: 截取前 N 条蛋白作为测试子集 ----
awk -v n="$N_PROTEINS" '/^>/{c++} c<=n' input/human.fasta > input/subset.fasta
echo "[S1b] 子集: $(grep -c '>' input/subset.fasta) 条蛋白 -> input/subset.fasta"

# ---- S3-mini: 用与生产完全相同的修饰配置建库 ----
# 注意: --var-mod-max 保持 3, 与生产配置一致 —— 阶段 1 验证的就是这套配置本身
# --out-lib 直接指定 .parquet, DIA-NN 2.x 支持自家 parquet 库格式, 便于后续用 pandas 检查
echo "[S3-mini] 建库开始: $(date)"
"$DIANN" \
  --fasta "$WORKDIR/input/subset.fasta" \
  --fasta-search \
  --predictor \
  --gen-spec-lib \
  --cut "K*,R*,!*P" \
  --missed-cleavages 2 \
  --min-pep-len 7 \
  --max-pep-len 30 \
  --min-pr-charge 2 \
  --max-pr-charge 4 \
  --min-pr-mz 300 \
  --max-pr-mz 1800 \
  --min-fr-mz 200 \
  --max-fr-mz 1800 \
  --fixed-mod "UniMod:4,57.021464,C" \
  --var-mod "UniMod:35,15.994915,M" \
  --var-mod "UniMod:1,42.010565,*n" \
  --var-mod "UniMod:21,79.966331,STY" \
  --var-mod-max 3 \
  --threads "$THREADS" \
  --out-lib "$WORKDIR/output/stage1_lib.parquet" \
  --out "$WORKDIR/output/stage1_lib_gen.tsv" \
  2>&1 | tee "$WORKDIR/output/stage1_lib_gen.log"
echo "[S3-mini] 建库结束: $(date)"

# ---- 验收 1: log 检查 ----
grep -qi "precursors generated" output/stage1_lib_gen.log \
  || echo "WARN: log 中未见 'precursors generated' 字样, 请人工检查 log"
# 排除正常输出中的 'mass error'/'Mass Error' 等误报
if grep -i "error" output/stage1_lib_gen.log | grep -vi "mass error" | grep -vi "mass accuracy" | grep -q .; then
    echo "WARN: log 中存在疑似 error 行:"
    grep -i "error" output/stage1_lib_gen.log | grep -vi "mass error" | grep -vi "mass accuracy"
fi

# ---- 验收 2: 定位实际生成的谱库文件 ----
# DIA-NN 可能按指定名落盘, 也可能改扩展名 (.speclib / .predicted.speclib)
LIB_FILE=""
for f in output/stage1_lib.parquet output/stage1_lib.speclib output/stage1_lib.predicted.speclib; do
    [[ -s "$f" ]] && LIB_FILE="$f" && break
done
if [[ -z "$LIB_FILE" ]]; then
    echo "FAIL: 未找到生成的谱库文件。output/ 实际内容:"
    ls -la output/
    exit 1
fi
echo "[验收] 谱库文件: $LIB_FILE ($(du -h "$LIB_FILE" | cut -f1))"

# 若只有二进制 .speclib, 转换为 parquet 以便检查
if [[ "$LIB_FILE" != *.parquet ]]; then
    echo "[转换] .speclib -> .parquet ..."
    "$DIANN" --lib "$LIB_FILE" --out-lib "$WORKDIR/output/stage1_lib.parquet" --threads "$THREADS" \
      2>&1 | tee -a "$WORKDIR/output/stage1_lib_gen.log"
    LIB_FILE="output/stage1_lib.parquet"
    [[ -s "$LIB_FILE" ]] || { echo "FAIL: 转换 parquet 失败"; exit 1; }
fi

# ---- S5-mini: 检查库中是否含 UniMod:21 前体 ----
module load Python/3.13.1-GCCcore-14.2.0
python3 "$(dirname "$0")/validate_unimod.py" "$WORKDIR/$LIB_FILE"
