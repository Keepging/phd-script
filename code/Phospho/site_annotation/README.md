# 磷酸化位点 ELM / NLS-NES / 已知结合注释

对 `site_meta.tsv` 的每个 `site_id`（同一 site_id 的多个 multiplicity 只算一次）做一行注释。
输入 11,826 行 → 10,873 个位点，3,768 个蛋白（全部为人源 Swiss-Prot accession，在 2026_03 中仍全部是主 accession）。

## 产物

| 文件 | 内容 |
|---|---|
| `site_annotation.tsv` | **当前版**（2026-10-19，`--live`）：UniProt 2026_03 序列 + PyPI 最新 PTM-POSE；ELM 仍为 2018/2019 快照（原因见下） |
| `site_annotation.sources.tsv` | 当前版实际读取的每个数据文件路径与 md5 |
| `site_annotation_snapshot2019.tsv` | 首版（2026-10-18）：UniProt 2023_05 镜像序列，其他同上；保留作对照 |
| `site_annotation_comparison.tsv` | 新旧两版逐列命中数对比（`compare_runs.py` 生成） |
| `data/LIVE_VERSIONS.tsv` | `fetch_data.sh --live` 每个资源的版本、下载日期、URL、成功/失败 |
| `fasta_batches/` | 首版导出的 FASTA 批次（**未随本次更新重新生成**，序列仍是 2023_05 版） |
| `annotate_sites.py` | 注释脚本（标准库，Python 3.8+）；`--live` 使用 `data/*/live/`，缺失时逐文件回退到快照 |
| `make_fasta_batches.py` | FASTA 导出/分批脚本 |
| `compare_runs.py` | 两次输出的逐列对比 |
| `fetch_data.sh` | 下载快照数据；`--live` 额外抓官方当前版本 |

重跑：

```bash
bash fetch_data.sh --live                              # 快照 + 当前版本，写 data/LIVE_VERSIONS.tsv
python annotate_sites.py --live --out site_annotation.tsv
python annotate_sites.py --out site_annotation_snapshot2019.tsv   # 复现首版（已核对逐字节一致）
python compare_runs.py site_annotation_snapshot2019.tsv site_annotation.tsv > site_annotation_comparison.tsv
```

## 数据来源与版本（当前版 `site_annotation.tsv`）

本次网络策略为 Full，`fetch_data.sh --live` 于 **2026-10-19** 运行，结果见 `data/LIVE_VERSIONS.tsv`：

| 数据 | 用途 | 来源 | 版本 / 日期 | 状态 |
|---|---|---|---|---|
| UniProtKB/Swiss-Prot 人源序列（20,431 条） | 序列、残基核对、正则扫描 | rest.uniprot.org stream `reviewed:true AND organism_id:9606` | **Release 2026_03（2026-09-02）** | 已更新 |
| UniProt `sec_ac.txt` | 次级→主 accession | ftp.uniprot.org current_release | 2026_03 | 已更新 |
| UniProt 人源基因名表 | 伙伴 accession→基因名 | rest.uniprot.org tsv | 2026_03（与快照表合并，当前版优先） | 已更新 |
| PTM-POSE 注释（PSP 结合/功能/过程/激酶、PTMInt、OmniPath writer） | `psp_*`、`kinase_psp`、`ptmint_binding`、`enzyme_omnipath_writer` | PyPI `ptm-pose` 最新版 | **0.4.0**（2026-05-19 上传，wheel sha256 `bc214ce0…`）；包内数据 `last_updated` 2025-03-30 | PyPI 无更新版本；与首版用的 wheel 完全相同，所用 GMT 逐字节一致 |
| ELM 类 / 实例 / 互作 / 类→结构域 | ELM 各列、NLS/NES | elm.eu.org | — | **下载失败**，回退到快照（见下） |
| PhosphoSitePlus 位点全表 | `psp_site_known`、文献计数 | protmapper 0.0.29 S3 镜像 | 2024-01-19 | PSP 无匿名下载；protmapper 仍为 0.0.29，无更新镜像 |

**ELM 未能更新。** elm.eu.org（194.94.44.88）对本容器出口在 TLS 握手阶段直接重置连接（`curl: (35) Recv failure: Connection reset by peer`），代理已放行（CONNECT 200），同属 EMBL 的 www.embl.org、www.ebi.ac.uk、www.huber.embl.de 均可访问，web.archive.org 同样被重置；明文 http 走不了代理。也没有找到比 2019 更新的 ELM 公开副本（OmniPath 只在互作表里用 ELM，不含正则和实例）。因此当前版的 ELM 部分仍是：

| ELM 文件 | 来源 | 版本 |
|---|---|---|
| 类正则 270 类 | SLiMSuite `data/elm2019.motifs`（269 类）+ SLiMEnrich `elm_classes.2018-07-17.tsv` 补 DOC_CYCLIN_1 | ELM 2019 / 2018-07-17 |
| 实验实例（人源 true positive 1,716 条） | SLiMEnrich `data/elm_instances.tsv` | 2018-07-17 |
| 互作表、类→结构域表 | SLiMEnrich `data/elm_interactions.tsv`、`elm_interaction_domains.tsv` | 2018 |

ELM 2024 版有 356 类、4,283 个实例；2019 之后 SH2 类重新整理过（新增如 `LIG_SH2_NCK_1`、`LIG_SH2_SFK_2`、`LIG_SH2_GRB2like`，据记忆），本快照只有 `LIG_SH2_GRB2/PTP2/SRC/STAT3/STAT5/STAT6`。在能访问 elm.eu.org 的机器上跑 `bash fetch_data.sh --live` 后再 `python annotate_sites.py --live`，脚本会自动改用 `data/elm/live/` 下的官方文件（此时只用官方 `elms_index.tsv` 的正则，不再用 SLiMSuite 快照）。

ELM 数据的许可是 ELM 学术许可协议（非商业学术使用，见 elm.eu.org/media/Elm_academic_license.pdf）；首版 README 写成 CC BY-NC-SA 4.0 有误，已更正。

**仍未使用的资源**：NLSdb、iPTMnet、dbPTM、EPSD、UniProt 位点特征。本次任务只要求更新 ELM、UniProt、PTM-POSE，这些未重新尝试。`nls_15` 只基于 ELM `TRG_NLS_*`，`nes_15` 只基于 `TRG_NES_CRM1_1`。

## 首版 `site_annotation_snapshot2019.tsv` 的来源

UniProtKB/Swiss-Prot 全库 FASTA、`sec_ac`、条目表来自 protmapper 0.0.29 S3 镜像（文件日期 2024-01-19，570,420 条，与 Release 2023_05 条目数一致）；ELM、PTM-POSE、PSP 位点表同上表。

## 序列版本核对：76 → 1

| residue_check | 首版（2023_05） | 当前版（2026_03） |
|---|---|---|
| OK | 10,797 | 10,872 |
| MISMATCH | 60 | 1 |
| OUT_OF_RANGE | 16 | 0 |

* 首版 24 个蛋白上的 76 个不匹配位点，在 2026_03 序列上全部匹配。
* 新出现 1 个：**Q6IQ22_S106**（RAB12）。RAB12 在 2026_03 升到序列版本 4，N 端延长 96 个残基（244 → 340 aa），你的数据用的是旧版本，旧 S106 对应新序列 S202。该位点序列相关列留空，PSP/PTMInt 列照常按 accession+位点填写。
* **Q2KHR2_S1225** 在两版里残基都是 S，但首版序列已过时、上下文完全不同（旧 flank `SYLTKSNsTGQINFS`，新 `NIAQRSQsVPLTVMM`），首版的 ELM 结果对它是错的，当前版已更正（新增 14-3-3 CanoR 正则命中）。这类“残基巧合相同”的情况残基核对查不出，只有换对版本才能消除。

## 新旧命中数对比

完整表见 `site_annotation_comparison.tsv`。下表只列有变化的列，其余列（全部 ELM 实例列、ELM 互作伙伴、NES、全部 PSP/PTMInt/OmniPath 列）新旧完全相同。

| 列 | 首版 | 当前版 | 新增位点 | 失去位点 |
|---|---|---|---|---|
| `elm_1433_regex` | 915 | 922 | 7 | 0 |
| `elm_other_lig_regex` | 502 | 507 | 5 | 0 |
| `elm_other_lig_regex_span` | 5,854 | 5,902 | 50 | 2 |
| `elm_other_type_regex` | 9,096 | 9,166 | 72 | 2 |
| `nls_15` | 1,129 | 1,144 | 16 | 1 |

变化几乎全部来自首版因残基不匹配而留空的 76 个位点，现在可以扫描了；另有 Q2KHR2_S1225 的更正，以及 Q6IQ22_S106 因新不匹配而清空。SH2、PTB 两类和所有实验实例列没有变化。由于 ELM 本身没有更新，这次对比只反映序列版本的影响。

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

命中数概览（当前版，10,873 位点）：SH2 实例 1、正则 88（全为 Y）；PTB 实例 5、正则 5；14-3-3 实例 17、正则 922；其他 LIG 磷酸依赖正则 507、宽松 span 5,902；NLS ±15 1,144；NES ±15 192；PSP 结合伙伴 479；PTMInt 263；PSP 已知位点 9,992。

## 方法细节

* 正则用 Python `re`，以 `(?=(…))` 前瞻取**所有重叠匹配**；ELM 正则里的括号即修饰残基标记，取捕获长度为 1 且为 S/T/Y 的组作为磷酸残基位置。
* 实例按 `Primary_Acc` 及 `Accessions` 列的全部 accession（去 isoform 后缀）索引；只用 `InstanceLogic = true positive`、`Organism = Homo sapiens`。标注在 isoform 上的实例（如 `Q13480-2`）仍按坐标匹配，但在 `elm_instance_detail` 里标 `isoform:`。
* PSP/PTMInt/OmniPath GMT 以 `ACC_S123` 查；GMT 里 isoform 特异条目（如 `O94875-2_Y59`）不会匹配到 canonical 位点。
* 同一 `site_id` 的多条 multiplicity 行只取第一条的 protein/gene。

## 网页服务器输入（14-3-3-Pred、GPS-PBS）

未爬取。`fasta_batches/` 是首版导出的（本次按要求未重新生成，序列为 2023_05 版，其中 24 个蛋白已过时）；需要新版时运行 `python make_fasta_batches.py --fasta data/uniprot/live/uniprot_sprot_human.fasta.gz --sec-ac data/uniprot/live/sec_ac.txt.gz`。默认每批 50 条、共 76 个文件，头行 `>ACC GENE len=N`。
两个服务器的单次上限本会话无法核对（站点不可达）：14-3-3-Pred（compbio.dundee.ac.uk/1433pred）网页版据记忆一次只收一条序列，批量请用其 standalone 脚本或 `--batch-size 1`（生成 `fasta_batches/single/ACC.fasta`）；GPS-PBS（Cuckoo 组）网页版请按页面提示的条数上限重跑 `make_fasta_batches.py --batch-size N`。最长序列 8,797 aa，注意服务器的长度上限（见 `batch_index.tsv` 的 length 列）。

## 引用

* ELM: Kumar M. et al. *NAR* 2024 (ELM 2024)；本分析实际使用 ELM 2018/2019 快照。
* UniProt: The UniProt Consortium, *NAR* 2025；序列为 Release 2026_03。
* SLiMSuite / SLiMEnrich: Edwards R.J. 等，https://github.com/slimsuite 。
* protmapper: Bachman J.A., Sorger P.K., Gyori B.M. *bioRxiv* 2019 / https://github.com/indralab/protmapper 。
* PTM-POSE: Crowl S. et al., Naegle lab, https://github.com/NaegleLab/PTM-POSE 。
* PhosphoSitePlus: Hornbeck P.V. et al. *NAR* 2015, 43:D512–D520（CC BY-NC-SA 3.0）。
* PTMInt: Hong X. et al. *Bioinformatics* 2023。
* OmniPath: Türei D. et al. *Mol Syst Biol* 2021。
