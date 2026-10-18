#!/usr/bin/env bash
# 重新获取 annotate_sites.py 所需的数据文件。
#   bash fetch_data.sh          # 从本次使用的镜像重新下载 (与 README 中记录的版本一致)
#   bash fetch_data.sh --live   # 另外从官方站点抓取最新 ELM 类/实例/互作表 和 UniProt 序列 (需要能访问 elm.eu.org / uniprot.org)
#
# 大文件 (Swiss-Prot FASTA、PSP 位点表) 不进 git，见 .gitignore。
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/uniprot data/psp data/elm data/ptm_pose_annotations

# ---- 1. UniProt / PhosphoSitePlus: protmapper (Sorger lab, INDRA) 维护的公开 S3 镜像 ----
#   https://github.com/indralab/protmapper  resources.py -> bucket "bigmech", key protmapper/<version>/<file>
PM_VER=0.0.29
S3=https://bigmech.s3.amazonaws.com/protmapper/$PM_VER
for f in uniprot_sprot.fasta.gz uniprot_sprot_varsplic.fasta.gz uniprot_sec_ac.txt.gz uniprot_entries.tsv.gz; do
  [ -s data/uniprot/$f ] || curl -fsSL -o data/uniprot/$f "$S3/$f"
done
[ -s data/psp/Phosphorylation_site_dataset.tsv.gz ] || curl -fsSL -o data/psp/Phosphorylation_site_dataset.tsv.gz "$S3/Phosphorylation_site_dataset.tsv.gz"

# ---- 2. ELM: SLiMSuite / SLiMEnrich (Edwards lab) 仓库里随软件分发的 ELM 快照 ----
if [ ! -s data/elm/elm2019.motifs ]; then
  tmp=$(mktemp -d)
  git clone -q --depth 1 https://github.com/slimsuite/slimsuite "$tmp/slimsuite"
  cp "$tmp/slimsuite/data/elm2019.motifs" data/elm/
  git clone -q --depth 1 https://github.com/slimsuite/SLiMEnrich "$tmp/SLiMEnrich"
  cp "$tmp/SLiMEnrich/example/elm_classes.2018-07-17.tsv" data/elm/
  cp "$tmp/SLiMEnrich/data/elm_instances.tsv"            data/elm/elm_instances.2018-07-17.tsv
  cp "$tmp/SLiMEnrich/data/elm_interactions.tsv"         data/elm/elm_interactions.2018.tsv
  cp "$tmp/SLiMEnrich/data/elm_interaction_domains.tsv"  data/elm/elm_interaction_domains.2018.tsv
  rm -rf "$tmp"
fi

# ---- 3. PhosphoSitePlus / PTMInt / OmniPath 位点注释: PTM-POSE (Naegle lab) PyPI 包内再分发的 GMT ----
if [ ! -s data/ptm_pose_annotations/PhosphoSitePlus_Interactions.gmt ]; then
  tmp=$(mktemp -d)
  pip download -q --no-deps --ignore-requires-python "ptm-pose==0.4.0" -d "$tmp"
  python3 - "$tmp" <<'EOF'
import glob, gzip, os, sys, zipfile
whl = glob.glob(os.path.join(sys.argv[1], "ptm_pose-*.whl"))[0]
z = zipfile.ZipFile(whl)
for name in z.namelist():
    if name.startswith("ptm_pose/Resource_Files/Annotations/") and name.endswith(".gmt.gz"):
        src, base = name.split("/")[-2:]
        if src in ("PhosphoSitePlus", "PTMInt", "OmniPath"):
            out = os.path.join("data/ptm_pose_annotations", "%s_%s" % (src, base[:-3]))
            open(out, "wb").write(gzip.decompress(z.read(name)))
open("data/ptm_pose_annotations/last_updated.txt", "wb").write(z.read("ptm_pose/Resource_Files/last_updated.txt"))
EOF
  rm -rf "$tmp"
fi

# ---- 4. 可选: 官方最新版 (本会话网络策略不可达，未在此验证) ----
if [ "${1:-}" = "--live" ]; then
  mkdir -p data/elm/live data/uniprot/live
  curl -fsSL -o data/elm/live/elms_index.tsv      "http://elm.eu.org/elms/elms_index.tsv"
  curl -fsSL -o data/elm/live/elm_instances.tsv   "http://elm.eu.org/instances.tsv?q=*&taxon=Homo%20sapiens&instance_logic=true%20positive"
  curl -fsSL -o data/elm/live/elm_interactions.tsv "http://elm.eu.org/interactions/as_tsv"
  curl -fsSL -o data/elm/live/elm_interaction_domains.tsv "http://elm.eu.org/interactiondomains.tsv"
  # 当前 Swiss-Prot 人源全集 (含本表全部 accession)
  curl -fsSL -o data/uniprot/live/uniprot_sprot_human.fasta.gz \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29"
  echo "最新数据已放在 data/*/live/，用法示例:"
  echo "  python annotate_sites.py --fasta data/uniprot/live/uniprot_sprot_human.fasta.gz \\"
  echo "      --elm-motifs '' --elm-classes data/elm/live/elms_index.tsv --elm-instances data/elm/live/elm_instances.tsv \\"
  echo "      --elm-interactions data/elm/live/elm_interactions.tsv --elm-interaction-domains data/elm/live/elm_interaction_domains.tsv"
fi
echo done
