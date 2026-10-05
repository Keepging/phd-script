# 04 脚本重复：逐对比较

比较对象：`code/` 里的副本（`code/<X>` == `D:\博士\<X>`）。字节比较用 `cmp`，行级比较用 Python `difflib.SequenceMatcher`（按行，不忽略空白）。
「大小」「修改时间」取自 `inventory_2026-09-25_gpt.xlsx`（本机原文件；`code/` 副本的 .py 因换行符 CRLF→LF 字节数略小，内容不变）。
Notebook 只比较 code cell 的 source，不比较 outputs / metadata。

## 表 A：`code/Phospho` 与 `code/Protein contour/Phospho` 同名文件

两个目录只有一个同名文件。

| 文件名 | 甲：`D:\博士\Phospho\` | 乙：`D:\博士\Protein contour\Phospho\` | 结论 | 差异大小 | 哪个更新 |
|---|---|---|---|---|---|
| 0210trying.ipynb | 20,681,513 B，2026-04-02 17:05；117 个 code cell，14,475 行代码 | 8,542,591 B，2026-02-28 20:51；78 个 code cell，9,898 行代码 | 有差异 | 乙的 9,898 行里 9,854 行在甲中原样存在（99.6%）；甲比乙多 4,600 行左右（见下表） | 甲（`D:\博士\Phospho\`）修改时间晚 33 天，cell 数多 39 个 |

0210trying.ipynb 甲 vs 乙 的 cell 级差异（cell 序号从 0 起，只列 code cell）：

| 差异类型 | 甲（Phospho）cell | 乙（Protein contour\Phospho）cell | 内容 |
|---|---|---|---|
| 替换 | 0–3（192 行） | 0–2（158 行） | 开头 cell 顺序不同；甲多一个 8 行 openpyxl 读 `Full_data.xlsx` 的 cell；「输入文件路径」cell 里 `JSON_FOLDER` 甲写 `D:\博士\Protein contour\Phospho\biophys_json`，乙写 `D:\博士\Phospho\biophys_json`（下同：乙里所有 `D:\博士\Phospho\...` 的数据路径在甲里都改成 `D:\博士\Protein contour\Phospho\...`） |
| 乙多一个 | （无） | 26（32 行） | 「快速启动：加载所有已保存的数据」cell，甲把它挪到了 cell 0 |
| 替换 | 44（128 行） | 44（128 行） | 只差 1 行：`Full_data.xlsx` 的路径前缀（同上） |
| 甲多两个 | 71–72（351 行） | （无） | 「Cell B-slide: Slide 7 专用 Within-Group High vs Low」「Cell C-slide: Slide 8 专用 Between-Group Heatmap」 |
| 替换 | 74–75（508 行） | 72（364 行） | 甲在「Cell D 非磷酸化对照组」前多一个 144 行的 `# %%` 绘图 cell；Cell D 本身相同 |
| 替换 | 80–116（4,081 行） | 77（1 行，空） | 甲末尾多 37 个 cell：Fig3-Proposal v4、0318 新图（BH 校正、合并 violin）、BANDLE 数据导出、±5 氨基酸窗口分析 Cell 1–10、Sequence Logo、gprofiler GO 富集、Top30 热图、whole-protein biophysical features、whole proteome 丰度、Kinase Analysis Step 1c 等 |

## 表 B：`code/protein_contour_wiki/raw/code` 与 `code/Protein contour/Phospho` 同名脚本

| 文件名 | `D:\博士\Protein contour\Phospho\` 大小 / 修改时间 | `D:\博士\protein_contour_wiki\raw\code\` 大小 / 修改时间 | 结论 | SHA-256 前 12 位（两份相同） |
|---|---|---|---|---|
| kinase_analysis.py | 15,249 B，2026-04-09 18:48:28 | 15,249 B，2026-04-09 18:48:28 | 完全相同（字节相同，修改时间也相同） | 62481cad0157 |
| kinase_analysis_v2.py | 21,894 B，2026-04-10 00:20:14 | 21,894 B，2026-04-10 00:20:14 | 完全相同 | 16989808a313 |
| kinase_correlation_analysis.py | 23,680 B，2026-04-10 16:34:17 | 23,680 B，2026-04-10 16:34:17 | 完全相同 | da84d9e29fc5 |
| length_matched_temporal_analysis.py | 19,175 B，2026-05-06 00:35:45 | 19,175 B，2026-05-06 00:35:45 | 完全相同 | ae8257ffaee0 |
| null_and_nonself_analysis.py | 11,159 B，2026-04-10 16:47:06 | 11,159 B，2026-04-10 16:47:06 | 完全相同 | 30a7805c2a4c |
| task6_s5c_detection.py | 14,557 B，2026-04-10 00:32:21 | 14,557 B，2026-04-10 00:32:21 | 完全相同 | cb860742aac3 |
| three_analyses.py | 12,904 B，2026-04-28 16:28:42 | 12,904 B，2026-04-28 16:28:42 | 完全相同 | 10eff6b5fa85 |
| top20_data_driven_kinase_heatmap.py | 9,157 B，2026-05-06 11:06:36 | 9,157 B，2026-05-06 11:06:36 | 完全相同 | bcb5d39809c2 |
| top20_qc_filter.py | 12,776 B，2026-05-06 11:21:04 | 12,776 B，2026-05-06 11:21:04 | 完全相同 | a7b8f7309170 |

`protein_contour_wiki/raw/code/SDRF/` 下的 3 个脚本（annotate_sdrf.py、fetch_metadata.py、pattern_scan.py）和 `protein_contour_wiki/wiki/figures/fig01_overview/fig01_overview.py` 在 `Protein contour/Phospho` 下没有同名文件。fig01_overview.py 与 top20_qc_filter.py 行数相同（323 行）但内容不同（323 行里只有 38 行相同）。

## 表 C：其他同名 / 同源脚本（不在任务指定的两组里，顺带列出）

| 甲 | 乙 | 结论 | 差异大小 | 哪个更新 |
|---|---|---|---|---|
| `D:\博士\Protein contour\Phospho\Test\test.py`（226,173 B，2026-02-18 14:41；6,395 行） | `D:\博士\Protein contour\Phospho\test.py`（282,853 B，2026-02-19 00:13；7,927 行） | 有差异 | 甲的 6,395 行里 5,937 行在乙里原样存在；乙在开头多 139 行（`# %%` 启动 cell）、中间多 4 处小插入（34+5+12 行），并把甲末尾的「Cell 26a 方案A / Cell 26b 方案B」约 540 行改写成约 1,880 行（终极修正版、5 组箱线图、Cell 30 S/T/Y Mobility Rate 等） | 乙（`Phospho\test.py`）晚 10 小时 |
| `D:\博士\Protein contour\Phospho\test.py`（7,927 行） | `D:\博士\Protein contour\Phospho\0220.py`（300,568 B，2026-02-20 13:39；8,366 行） | 有差异，同源 | 甲的 7,927 行全部在乙中原样存在；乙末尾多 2 个 cell 共 439 行 | 乙（0220.py） |
| `D:\博士\Protein contour\Phospho\0220.py`（8,366 行） | `D:\博士\Protein contour\Phospho\0228.py`（365,442 B，2026-02-28 20:30；10,089 行） | 有差异，同源 | 乙删了甲的 22 行（一个「看看目前 notebook 里有哪些 dataframe 变量」cell 和 3 条分隔注释），末尾多 2 个 cell 共 1,745 行 | 乙（0228.py） |
| `D:\博士\Protein contour\Phospho\0228.py`（10,089 行） | `D:\博士\Protein contour\Phospho\0210trying.ipynb`（9,898 行代码） | 有差异，同源 | 乙的 9,898 行里 9,896 行在甲中原样存在；甲比乙多约 190 行（0228.py 是这份 notebook 加上末尾几个 cell 的导出） | 甲（0228.py，20:30）与乙（notebook，20:51）同一天；notebook 晚 21 分钟 |
| `D:\博士\Protein contour\Phospho\0228.py`（10,089 行） | `D:\博士\Phospho\0331.py`（513,898 B，2026-03-31 15:24；13,989 行） | 有差异，同源 | 甲的 10,089 行里 10,005 行在乙中原样存在；差异：乙把 `D:\博士\Phospho\` 数据路径改成 `D:\博士\Protein contour\Phospho\`（BASE_PATH、JSON_FOLDER、OUTPUT_FILE、Full_data.xlsx 共 6 处），删了「快速启动」cell（34 行）和「Cell E 汇总检查」（36 行），末尾多约 3,900 行 | 乙（0331.py） |
| `D:\博士\Phospho\0331.py`（13,989 行） | `D:\博士\Phospho\0210trying.ipynb`（14,475 行代码） | 有差异，同源 | 甲的 13,989 行里 13,707 行在乙中原样存在（98%） | 乙（notebook，2026-04-02） |

同源链（按修改时间）：`Test\test.py`（02-18）→ `test.py`（02-19）→ `0220.py`（02-20）→ `0228.py` / `Protein contour\Phospho\0210trying.ipynb`（02-28）→ `Phospho\0331.py`（03-31）→ `Phospho\0210trying.ipynb`（04-02）。每一版基本是上一版加上末尾新 cell；仅 0228 → 0331 这一步把数据根目录从 `D:\博士\Phospho` 改为 `D:\博士\Protein contour\Phospho`。
