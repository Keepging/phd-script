# 磷酸化位点 ELM / NLS-NES / 已知结合注释

对 `site_meta.tsv` 的每个 `site_id`（同一 site_id 的多个 multiplicity 只算一次）做一行注释。
输入 11,826 行 → 10,873 个位点，3,768 个蛋白（全部为 Swiss-Prot accession）。

## 产物

| 文件 | 内容 |
|---|---|
| `site_annotation.tsv` | 每位点一行，列定义见下 |
| `fasta_batches/all_accessions.fasta` | 全部 3,768 条序列 |
| `fasta_batches/batch_001.fasta` … `batch_076.fasta` | 每批 50 条，供网页服务器手动提交 |
| `fasta_batches/batch_index.tsv` | accession → 批次文件、长度、位点数、残基不一致的位点 |
| `annotate_sites.py` | 注释脚本（标准库，Python 3.8+） |
| `make_fasta_batches.py` | FASTA 导出/分批脚本，`--batch-size N` 改批量，`--batch-size 1` 每个 accession 一个文件 |
| `fetch_data.sh` | 重新下载数据；`--live` 额外抓官方最新版 |
| `data/` | 数据快照（大文件不进 git，由 `fetch_data.sh` 下载） |

重跑：

```bash
bash fetch_data.sh                 # 只在 data/uniprot、data/psp 缺文件时需要
python annotate_sites.py           # -> site_annotation.tsv
python make_fasta_batches.py --batch-size 50
```

## 数据来源与版本

本会话的云环境网络策略只放行包管理器、GitHub 和 S3/GCS，**elm.eu.org、uniprot.org、omnipathdb.org、ebi.ac.uk、iptmnet、dbPTM、EPSD、NLSdb(rostlab.org) 全部被拦截**。
因此全部数据来自这些资源被合法再分发的副本，版本比官方当前版本旧，见下表。想换成最新版：在能上网的机器跑 `bash fetch_data.sh --live`，再按脚本末尾打印的命令重跑（`--live` 部分的 URL 未在本会话验证）。

| 数据 | 用途 | 再分发来源 | 版本 | 许可 |
|---|---|---|---|---|
| UniProtKB/Swiss-Prot 全库 FASTA、`sec_ac.txt`、人源条目表 | 序列、残基核对、次级 accession、伙伴基因名 | protmapper 0.0.29 (INDRA / Sorger lab) 的公开 S3 镜像 `s3://bigmech/protmapper/0.0.29/` | 镜像文件时间戳 2024-01-19；Swiss-Prot 共 570,420 条，与 UniProt release **2023_05** 的条目数一致（仅按条目数推断） | CC BY 4.0 |
| ELM 类正则 (270 类) | 正则匹配 | SLiMSuite 仓库 `data/elm2019.motifs`（ELM 2019 快照，269 类）+ SLiMEnrich 仓库 `example/elm_classes.2018-07-17.tsv` 补 1 类（DOC_CYCLIN_1） | ELM 2019 / 2018-07-17 | ELM: CC BY-NC-SA 4.0 |
| ELM 实验实例 (人源 true positive 1,716 条) | “实验实例命中” | SLiMEnrich 仓库 `data/elm_instances.tsv` | ELM 下载日期 2018-07-17 | 同上 |
| ELM motif–domain 互作表、类→结构域表 | 实例对应的结合伙伴 / 结构域 | SLiMEnrich 仓库 `data/elm_interactions.tsv`, `elm_interaction_domains.tsv` | 2018 | 同上 |
| PhosphoSitePlus 调控位点 (结合伙伴 / 功能 / 过程 / 激酶) | `psp_*`, `kinase_psp` | PTM-POSE 0.4.0 (Naegle lab, PyPI) 包内 `Resource_Files/Annotations/PhosphoSitePlus/*.gmt` | PTM-POSE 资源日期 2025-03-30 | PSP: CC BY-NC-SA 3.0，非商业使用，须引用 Hornbeck et al. 2015 |
| PhosphoSitePlus 磷酸化位点全表 | `psp_site_known`, 文献计数 | protmapper S3 镜像 `Phosphorylation_site_dataset.tsv.gz` | 2024-01-19 快照 | 同上 |
| PTMInt (PTM 介导的互作) | `ptmint_binding` | PTM-POSE 0.4.0 `PTMInt/Interactions.gmt` | 2025-03-30 | PTMInt 免费学术使用 |
| OmniPath 酶–底物 (writer) | `enzyme_omnipath_writer` | PTM-POSE 0.4.0 `OmniPath/Writer_Enzyme.gmt` | 2025-03-30 | 各子库许可不一，学术用 |

**没找到的**：

* **NLSdb**：只在 rostlab.org 提供，本环境不可达，也没有可再分发副本。`nls_15` 只基于 ELM `TRG_NLS_*`（Bipartite_1、MonoCore_2、MonoExtC_3、MonoExtN_4），`nes_15` 只基于 `TRG_NES_CRM1_1`。
* **iPTMnet、dbPTM、EPSD、UniProt 位点特征**：站点被拦截，且没有找到合法镜像；未使用。
* **OmniPath 直接下载**：被拦截；只用了 PTM-POSE 再分发的 writer-enzyme 表（里面混有非激酶的上游调控者，如细胞因子，见列说明）。
* ELM 的 SH2 类在 2019 之后重新整理扩充过（据记忆新增了 `LIG_SH2_NCK_1`、`LIG_SH2_SFK_2`、`LIG_SH2_GRB2like` 等），本快照只有 `LIG_SH2_GRB2 / PTP2 / SRC / STAT3 / STAT5 / STAT6`；`LIG_14-3-3_*` 和 `LIG_PTB_*` 的 2019 版正则与当前版本基本一致，但建议最终投稿前用 `--live` 重跑一次。

## 序列版本核对

每个位点都核对了 Swiss-Prot 快照中该位置的残基：

| residue_check | 位点数 |
|---|---|
| OK | 10,797 |
| MISMATCH（残基不同） | 60 |
| OUT_OF_RANGE（超出序列长度） | 16 |

不一致的 76 个位点来自 24 个蛋白，说明你的 UniProt 版本比 2023_05 新、这些序列已更新：
A9YTQ3, P08195, P49761, P78312, P78549, Q13636, Q14160, Q14241, Q14669, Q2KHR2, Q52LA3, Q5SNV9, Q8WX92, Q92538, Q96GQ7, Q96JG9, Q96JP2, Q96RK0, Q9BWU0, Q9BYK8, Q9HB09, Q9NVI7, Q9NYZ3, Q9ULE0。
这些位点所有依赖序列的列（ELM 正则/实例、NLS/NES、flank）留空，只填按 accession+位点查到的 PSP/PTMInt 列；FASTA 批次里这 24 条仍是旧序列（`batch_index.tsv` 最后一列标出）。用 `--fasta` 指向你自己的 UniProt 导出即可修正。

## `site_annotation.tsv` 列定义

通用：`site_id, protein, gene, residue, position, seq_source, seq_length, residue_check, flank_pm7`（±7 残基，中心小写，`_` 补齐）。

ELM 命中格式一律 `类名[起-止]`（1-based 闭区间），多个用 `;` 分隔。

| 列 | 含义 |
|---|---|
| `elm_sh2_instance` / `elm_sh2_regex` | `LIG_SH2_*`：ELM 人源实验实例区间包含该位点 / 正则匹配且**该位点正是正则括号标记的 pY** |
| `elm_ptb_instance` / `elm_ptb_regex` | `LIG_PTB_*`，同上（`LIG_PTB_Phospho_1` 的 pY；`LIG_PTB_Apo_2` 无磷酸残基标记，只会出现在 span 列） |
| `elm_1433_instance` / `elm_1433_regex` | `LIG_14-3-3_*`，同上（pS/pT） |
| `elm_other_lig_instance` | 其他 `LIG_*` 实验实例区间包含该位点 |
| `elm_other_lig_regex` | 其他 `LIG_*` 正则匹配且该位点是括号标记的 S/T/Y（磷酸依赖型，如 `LIG_FHA_1/2`、`LIG_BRCT_BRCA1_1`、`LIG_WW_3` 等） |
| `elm_other_lig_regex_span` | 其他 `LIG_*` 正则匹配区间包含该位点，但该位点不是（或该类没有）标记的磷酸残基（如 `LIG_SH3_*`、`LIG_PDZ_*` 内恰好有个 S/T）——即“包含该磷酸位点的类名”的宽松版本 |
| `elm_other_type_instance` / `elm_other_type_regex` | 非 LIG 类（DEG/DOC/MOD/TRG/CLV）的实例包含位点 / 正则且位点是标记的 S/T/Y（如 `MOD_CDK_SPK_2`、`DEG_SCF_TRCP1_1`、`DOC_WW_Pin1_4`）。`MOD_*` 几乎每个位点都会中，供参考 |
| `elm_instance_detail` | 命中的全部实例：`类[区间](ELMI编号,PMID:…,ok或span_regex_mismatch[,isoform:…])`。`ok` = 该区间在当前序列上仍能被该类正则整段匹配；`span_regex_mismatch` 提示 2018 坐标可能已漂移 |
| `elm_instance_partner` | ELM 互作表里该实例实验验证的结合蛋白：`类->基因(accession)[Pfam]`，仅人–人 |
| `elm_class_domain` | 命中实例所属类的结合结构域（类级别，Pfam id:名） |
| `nls_15`, `nls_15_instance`, `nls_15_regex` | 位点 ±15 残基窗口内是否有 `TRG_NLS_*`（实例或正则区间与窗口重叠）及明细 |
| `nes_15`, `nes_15_instance`, `nes_15_regex` | 同上，`TRG_NES_*` |
| `psp_binding` | PhosphoSitePlus 调控位点 ON_PROT_INTERACT：`伙伴(INDUCES/DISRUPTS)` |
| `ptmint_binding` | PTMInt：`伙伴基因(accession)->Enhance/Inhibit` |
| `psp_function`, `psp_process` | PSP 调控位点的功能 / 生物过程注释 |
| `kinase_psp` | PSP 激酶–底物表里该位点的激酶 |
| `enzyme_omnipath_writer` | OmniPath 酶–底物（writer）表的上游酶；**包含文献挖掘来的非激酶调控者**（如 IFNG、IL1B），需自行甄别 |
| `psp_site_known`, `psp_lit_LT/MS/CST` | 该位点是否在 PSP 位点全表（2024-01 快照）中；低通量文献 / 质谱文献 / CST 内部质谱计数 |

命中数概览（10,873 位点）：SH2 实例 1、正则 88（全为 Y）；PTB 实例 5、正则 5；14-3-3 实例 17、正则 915；其他 LIG 磷酸依赖正则 502、宽松 span 5,854；NLS ±15 1,129；NES ±15 192；PSP 结合伙伴 479；PTMInt 263；PSP 已知位点 9,992。

## 方法细节

* 正则用 Python `re`，以 `(?=(…))` 前瞻取**所有重叠匹配**；ELM 正则里的括号即修饰残基标记，取捕获长度为 1 且为 S/T/Y 的组作为磷酸残基位置。
* 实例按 `Primary_Acc` 及 `Accessions` 列的全部 accession（去 isoform 后缀）索引；只用 `InstanceLogic = true positive`、`Organism = Homo sapiens`。标注在 isoform 上的实例（如 `Q13480-2`）仍按坐标匹配，但在 `elm_instance_detail` 里标 `isoform:`。
* PSP/PTMInt/OmniPath GMT 以 `ACC_S123` 查；GMT 里 isoform 特异条目（如 `O94875-2_Y59`）不会匹配到 canonical 位点。
* 同一 `site_id` 的多条 multiplicity 行只取第一条的 protein/gene。

## 网页服务器输入（14-3-3-Pred、GPS-PBS）

未爬取。`fasta_batches/` 默认每批 50 条、共 76 个文件，头行 `>ACC GENE len=N`。
两个服务器的单次上限本会话无法核对（站点不可达）：14-3-3-Pred（compbio.dundee.ac.uk/1433pred）网页版据记忆一次只收一条序列，批量请用其 standalone 脚本或 `--batch-size 1`（生成 `fasta_batches/single/ACC.fasta`）；GPS-PBS（Cuckoo 组）网页版请按页面提示的条数上限重跑 `make_fasta_batches.py --batch-size N`。最长序列 8,797 aa，注意服务器的长度上限（见 `batch_index.tsv` 的 length 列）。

## 引用

* ELM: Kumar M. et al. *NAR* 2024 (ELM 2024)；本快照为 ELM 2018/2019。
* SLiMSuite / SLiMEnrich: Edwards R.J. 等，https://github.com/slimsuite 。
* protmapper: Bachman J.A., Sorger P.K., Gyori B.M. *bioRxiv* 2019 / https://github.com/indralab/protmapper 。
* PTM-POSE: Crowl S. et al., Naegle lab, https://github.com/NaegleLab/PTM-POSE 。
* PhosphoSitePlus: Hornbeck P.V. et al. *NAR* 2015, 43:D512–D520（CC BY-NC-SA 3.0）。
* PTMInt: Hong X. et al. *Bioinformatics* 2023。
* OmniPath: Türei D. et al. *Mol Syst Biol* 2021。
