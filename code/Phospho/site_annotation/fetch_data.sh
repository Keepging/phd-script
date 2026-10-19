#!/usr/bin/env bash
# 重新获取 annotate_sites.py 所需的数据文件。
#   bash fetch_data.sh          # 从本次使用的镜像重新下载 (与 README 中记录的版本一致)
#   bash fetch_data.sh --live   # 另外抓官方当前版本: ELM (elm.eu.org)、UniProt 当前 release、PyPI 最新 PTM-POSE -> data/*/live/，版本记入 data/LIVE_VERSIONS.tsv
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

# ---- 4. --live: 官方当前版本，每个资源的版本 / 下载日期写入 data/LIVE_VERSIONS.tsv ----
if [ "${1:-}" = "--live" ]; then
  set +e
  TODAY=$(date -u +%Y-%m-%d)
  V=data/LIVE_VERSIONS.tsv
  printf 'resource\tfile\tversion\tdownload_date_utc\tsource_url\tstatus\n' > $V
  log() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$TODAY" "$4" "$5" >> $V; }
  mkdir -p data/elm/live data/uniprot/live data/ptm_pose_annotations/live

  # 4a. ELM (elm.eu.org; 只走 https，代理不支持明文 http)
  ELM=https://elm.eu.org
  for pair in \
      "elms_index.tsv|$ELM/elms/elms_index.tsv" \
      "elm_instances.tsv|$ELM/instances.tsv?q=*&taxon=Homo%20sapiens&instance_logic=true%20positive" \
      "elm_interactions.tsv|$ELM/interactions/as_tsv" \
      "elm_interaction_domains.tsv|$ELM/interactiondomains.tsv"; do
    f=${pair%%|*}; u=${pair#*|}
    if curl -fsSL --retry 2 -m 300 -o data/elm/live/$f.part "$u" && [ -s data/elm/live/$f.part ]; then
      mv data/elm/live/$f.part data/elm/live/$f
      ver=$(grep -m2 -E '^#ELM_.*(Version|Date)' data/elm/live/$f | tr -d '\r' | paste -sd ' ' -)
      log ELM "data/elm/live/$f" "${ver:-see file header}" "$u" ok
    else
      rm -f data/elm/live/$f.part
      log ELM "data/elm/live/$f" "" "$u" "FAILED: $(curl -sS -m 30 -o /dev/null "$u" 2>&1 | head -1)"
      echo "ELM 下载失败: $u (annotate_sites.py --live 会回退到 data/elm/ 下的 2018/2019 快照)" >&2
    fi
  done

  # 4b. UniProt 当前 release: 人源 Swiss-Prot 序列 + 基因名表 + sec_ac
  Q='%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29'
  U1="https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=$Q"
  U2="https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&fields=accession,gene_primary,length,reviewed&query=$Q"
  U3="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt"
  U4="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"
  curl -fsSL -m 900 -D data/uniprot/live/fasta.headers -o data/uniprot/live/uniprot_sprot_human.fasta.gz "$U1" \
    && curl -fsSL -m 900 -o data/uniprot/live/uniprot_entries_human.tsv.gz "$U2" \
    && curl -fsSL -m 900 "$U3" | gzip -c > data/uniprot/live/sec_ac.txt.gz \
    && curl -fsSL -m 60 -o data/uniprot/live/reldate.txt "$U4"
  st=$?
  rel=$(grep -i -m1 '^x-uniprot-release:' data/uniprot/live/fasta.headers | cut -d' ' -f2 | tr -d '\r')
  rdate=$(grep -i -m1 '^x-uniprot-release-date:' data/uniprot/live/fasta.headers | cut -d' ' -f2- | tr -d '\r')
  [ $st -eq 0 ] && s=ok || s=FAILED
  log UniProt data/uniprot/live/uniprot_sprot_human.fasta.gz "release $rel ($rdate); $(zcat data/uniprot/live/uniprot_sprot_human.fasta.gz | grep -c '^>') entries" "$U1" $s
  log UniProt data/uniprot/live/uniprot_entries_human.tsv.gz "release $rel" "$U2" $s
  log UniProt data/uniprot/live/sec_ac.txt.gz "$(grep -m1 Swiss-Prot data/uniprot/live/reldate.txt | tr -d '\r')" "$U3" $s

  # 4c. PTM-POSE: PyPI 最新版 (不锁版本)
  tmp=$(mktemp -d)
  if pip download -q --no-deps --ignore-requires-python ptm-pose -d "$tmp"; then
    whl=$(ls "$tmp"/ptm_pose-*.whl | head -1)
    python3 - "$whl" <<'EOF'
import gzip, os, sys, zipfile
z = zipfile.ZipFile(sys.argv[1])
for name in z.namelist():
    if name.startswith("ptm_pose/Resource_Files/Annotations/") and name.endswith(".gmt.gz"):
        src, base = name.split("/")[-2:]
        if src in ("PhosphoSitePlus", "PTMInt", "OmniPath"):
            open(os.path.join("data/ptm_pose_annotations/live", "%s_%s" % (src, base[:-3])), "wb").write(gzip.decompress(z.read(name)))
open("data/ptm_pose_annotations/live/last_updated.txt", "wb").write(z.read("ptm_pose/Resource_Files/last_updated.txt"))
EOF
    pv=$(basename "$whl" | cut -d- -f2)
    log PTM-POSE data/ptm_pose_annotations/live/ "ptm-pose $pv; wheel sha256 $(sha256sum "$whl" | cut -c1-16)…; Resource_Files/last_updated $(cat data/ptm_pose_annotations/live/last_updated.txt)" "https://pypi.org/project/ptm-pose/$pv/" ok
  else
    log PTM-POSE data/ptm_pose_annotations/live/ "" "https://pypi.org/project/ptm-pose/" FAILED
  fi
  rm -rf "$tmp"

  # 4d. PhosphoSitePlus 位点全表: 无官方匿名下载，仍用 protmapper 镜像 (protmapper 最新版仍为 $PM_VER)
  log PhosphoSitePlus data/psp/Phosphorylation_site_dataset.tsv.gz "protmapper $PM_VER S3 mirror, file dated 2024-01-19" "$S3/Phosphorylation_site_dataset.tsv.gz" "unchanged (no newer redistributed copy)"

  set -e
  cat $V
  echo "运行: python annotate_sites.py --live --out site_annotation.tsv"
fi
echo done
