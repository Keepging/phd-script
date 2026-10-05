# 01 规模表：按目录聚合（到第 3 层）

来源：`inventory_2026-09-25_gpt.xlsx` 工作表「全部文件」，共 61044 行（文件）。

口径：
- 层级：第 1 层 = `D:\博士\<文件夹>`（4 个根），第 2 层 = 根下一级子目录，第 3 层 = 再下一级。
- 「文件数」= 该目录及其所有子目录内的文件数（递归）；「直接文件数」= 直接位于该目录、不在子目录内的文件数。父目录的数字包含子目录。
- 「总大小」按 xlsx 的字节数求和。
- 「主要类型」= xlsx「类型」列前 3 名（含数量）；「主要扩展名」= 路径扩展名前 3 名。
- 活 / 死 / 不确定 = xlsx「状态」列原值（本表未做任何重判，重判结果见 02）。
- 按文件数降序；同数时按层级、路径排。只列出 xlsx 里实际有文件的目录。

全部合计：61044 个文件，63.0 GB；活 385，死 38875，不确定 21784。

| 目录 | 层 | 文件数 | 直接文件数 | 总大小 | 主要类型 | 主要扩展名 | 活 | 死 | 不确定 |
|---|---|---:|---:|---:|---|---|---:|---:|---:|
| D:\博士\Protein contour | 1 | 60352 | 21 | 61.8 GB | 数据表 58359、文档 958、图 542 | .tsv 24950、.json 22965、.txt 5172 | 366 | 38518 | 21468 |
| D:\博士\Protein contour\Zhihan | 2 | 37519 | 1 | 47.8 GB | 数据表 37474、脚本 29、文档 8 | .tsv 24945、.parquet 5042、.txt 4990 | 29 | 37446 | 44 |
| D:\博士\Protein contour\Phospho | 2 | 22477 | 253 | 3.8 GB | 数据表 20653、文档 913、图 470 | .json 20461、.md 861、.png 399 | 326 | 842 | 21309 |
| D:\博士\Protein contour\Phospho\biophys_json | 3 | 20421 | 20421 | 2.6 GB | 数据表 20421 | .json 20420、.csv 1 | 0 | 0 | 20421 |
| D:\博士\Protein contour\Zhihan\DIANN_v2.2.0_phospho | 3 | 14129 | 0 | 15.9 GB | 数据表 14128、其他 1 | .tsv 9972、.txt 1662、.parquet 1662 | 0 | 14125 | 4 |
| D:\博士\Protein contour\Zhihan\DIANN_v2.2.0_acetyl | 3 | 12467 | 0 | 15.8 GB | 数据表 12466、其他 1 | .tsv 8310、.txt 1662、.parquet 1662 | 0 | 12467 | 0 |
| D:\博士\Protein contour\Zhihan\DIANN_v2.2.0 | 3 | 10805 | 0 | 15.7 GB | 数据表 10804、其他 1 | .tsv 6648、.txt 1662、.parquet 1662 | 0 | 10796 | 9 |
| D:\博士\Protein contour\Phospho\.agents | 3 | 1409 | 0 | 173.6 MB | 文档 867、图 202、脚本 185 | .md 843、.png 202、.py 180 | 161 | 807 | 441 |
| D:\博士\protein_contour_wiki | 1 | 507 | 8 | 341.4 MB | 其他 320、文档 125、数据表 44 | (无扩展名) 304、.md 94、.pdf 29 | 16 | 335 | 156 |
| D:\博士\protein_contour_wiki\.git | 2 | 318 | 7 | 147.8 MB | 其他 318 | (无扩展名) 304、.sample 14 | 0 | 318 | 0 |
| D:\博士\protein_contour_wiki\.git\objects | 3 | 289 | 0 | 147.7 MB | 其他 289 | (无扩展名) 289 | 0 | 289 | 0 |
| D:\博士\Phospho | 1 | 163 | 141 | 434.7 MB | 图 131、数据表 29、脚本 3 | .tiff 66、.png 65、.csv 29 | 3 | 22 | 138 |
| D:\博士\Protein contour\scope3p data | 2 | 153 | 0 | 1.9 GB | 数据表 137、图 7、文档 6 | .txt 129、.csv 8、.png 7 | 3 | 129 | 21 |
| D:\博士\Protein contour\Zhihan\Test | 3 | 117 | 51 | 387.6 MB | 数据表 75、脚本 29、文档 8 | .parquet 56、.py 16、.tsv 15 | 29 | 58 | 30 |
| D:\博士\Protein contour\scope3p data\116 Case raw peptide data | 3 | 116 | 116 | 1.7 GB | 数据表 116 | .txt 116 | 0 | 116 | 0 |
| D:\博士\protein_contour_wiki\wiki | 2 | 91 | 6 | 1.5 MB | 文档 88、图 2、脚本 1 | .md 88、.png 1、.py 1 | 3 | 3 | 85 |
| D:\博士\Protein contour\Phospho\SpatialProteoDynamics.github.io-v1.0 | 3 | 90 | 0 | 57.5 MB | 脚本 71、文档 8、图 4 | .js 68、.html 6、.png 3 | 3 | 8 | 79 |
| D:\博士\Protein contour\Phospho\multi_kinase_analysis | 3 | 90 | 90 | 54.6 MB | 图 84、数据表 4、脚本 2 | .png 42、.tiff 42、.csv 4 | 8 | 0 | 82 |
| D:\博士\Protein contour\TCA_Succciny | 2 | 73 | 0 | 294.7 MB | 数据表 49、图 17、脚本 5 | .csv 26、.png 11、.rds 10 | 5 | 26 | 42 |
| D:\博士\protein_contour_wiki\raw | 2 | 72 | 0 | 189.8 MB | 数据表 31、文档 29、脚本 12 | .pdf 28、.csv 27、.py 12 | 12 | 3 | 57 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803 | 2 | 71 | 3 | 12.3 MB | 图 46、文档 13、数据表 10 | .svg 25、.png 21、.md 10 | 0 | 38 | 33 |
| D:\博士\Protein contour\TCA_Succciny\代码 | 3 | 49 | 32 | 126.3 MB | 数据表 31、图 11、脚本 5 | .csv 19、.png 9、.rds 9 | 5 | 14 | 30 |
| D:\博士\Protein contour\Phospho\martinez_network_check | 3 | 44 | 41 | 25.4 MB | 数据表 16、脚本 14、图 10 | .csv 16、.py 14、.png 9 | 32 | 1 | 11 |
| D:\博士\Protein contour\Phospho\bandle_input | 3 | 32 | 32 | 347.8 MB | 数据表 32 | .csv 26、.rdata 5、.txt 1 | 3 | 5 | 24 |
| D:\博士\protein_contour_wiki\raw\data | 3 | 31 | 15 | 60.4 MB | 数据表 31 | .csv 27、.xlsx 4 | 0 | 2 | 29 |
| D:\博士\protein_contour_wiki\raw\papers | 3 | 28 | 28 | 128.3 MB | 文档 28 | .pdf 27、.pptx 1 | 0 | 1 | 27 |
| D:\博士\protein_contour_wiki\wiki\literature | 3 | 25 | 25 | 51.9 KB | 文档 25 | .md 25 | 0 | 0 | 25 |
| D:\博士\Protein contour\Phospho\Result | 3 | 23 | 23 | 21.0 MB | 数据表 12、图 11 | .csv 12、.png 11 | 0 | 3 | 20 |
| D:\博士\博后_大导开会会议记录 | 1 | 22 | 22 | 440.3 MB | 会议记录 22 | .txt 15、.wav 5、.pdf 1 | 0 | 0 | 22 |
| D:\博士\Phospho\bandle_input | 2 | 22 | 22 | 14.0 MB | 数据表 22 | .csv 22 | 0 | 22 | 0 |
| D:\博士\Protein contour\Phospho\kinase_correlation | 3 | 22 | 19 | 11.0 MB | 图 13、数据表 9 | .png 13、.csv 9 | 19 | 0 | 3 |
| D:\博士\Protein contour\Phospho\layerA_full | 3 | 22 | 22 | 68.5 MB | 数据表 19、文档 1、脚本 1 | .csv 12、.json 5、.md 1 | 20 | 0 | 2 |
| D:\博士\protein_contour_wiki\wiki\analyses | 3 | 21 | 21 | 48.6 KB | 文档 21 | .md 21 | 0 | 0 | 21 |
| D:\博士\Protein contour\Phospho_HeLa_CellCycle | 2 | 19 | 19 | 322.5 MB | 数据表 18、文档 1 | .txt 18、.pdf 1 | 3 | 14 | 2 |
| D:\博士\Protein contour\Proteomeandfractionatedproteome | 2 | 19 | 0 | 5.6 GB | 数据表 18、文档 1 | .txt 18、.pdf 1 | 0 | 19 | 0 |
| D:\博士\Protein contour\Proteomeandfractionatedproteome\txt | 3 | 19 | 19 | 5.6 GB | 数据表 18、文档 1 | .txt 18、.pdf 1 | 0 | 19 | 0 |
| D:\博士\Protein contour\scope3p data\Scop3p_data | 3 | 19 | 13 | 67.2 MB | 数据表 13、文档 4、脚本 2 | .txt 13、.pdf 4、.r 1 | 2 | 13 | 4 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\exports | 3 | 18 | 3 | 7.6 MB | 图 15、文档 3 | .png 15、.pptx 3 | 0 | 10 | 8 |
| D:\博士\Protein contour\TCA_Succciny\Result | 3 | 15 | 15 | 74.2 MB | 数据表 9、图 6 | .csv 6、.tiff 4、.xlsx 3 | 0 | 6 | 9 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\backup | 3 | 15 | 0 | 50.2 KB | 图 15 | .svg 15 | 0 | 15 | 0 |
| D:\博士\protein_contour_wiki\.git\hooks | 3 | 14 | 14 | 25.2 KB | 其他 14 | .sample 14 | 0 | 14 | 0 |
| D:\博士\protein_contour_wiki\raw\code | 3 | 12 | 9 | 158.9 KB | 脚本 12 | .py 12 | 12 | 0 | 0 |
| D:\博士\protein_contour_wiki\.obsidian | 2 | 11 | 6 | 748.1 KB | 数据表 8、脚本 2、其他 1 | .json 8、.js 1、.sh 1 | 1 | 1 | 9 |
| D:\博士\protein_contour_wiki\wiki\findings | 3 | 11 | 11 | 30.3 KB | 文档 11 | .md 11 | 0 | 0 | 11 |
| D:\博士\Protein contour\Phospho\scaffold_check | 3 | 10 | 10 | 1.1 MB | 数据表 7、脚本 2、文档 1 | .csv 5、.py 2、.json 2 | 4 | 0 | 6 |
| D:\博士\Protein contour\Phospho\anomaly_first_codex | 3 | 9 | 9 | 51.4 KB | 数据表 6、文档 2、脚本 1 | .csv 5、.md 2、.txt 1 | 7 | 0 | 2 |
| D:\博士\Protein contour\Phospho\kinase_results | 3 | 9 | 0 | 2.9 MB | 图 7、数据表 2 | .png 7、.csv 2 | 9 | 0 | 0 |
| D:\博士\Protein contour\scope3p data\116_integration | 3 | 9 | 1 | 5.5 MB | 图 4、数据表 4、文档 1 | .png 4、.csv 4、.md 1 | 0 | 0 | 9 |
| D:\博士\Protein contour\scope3p data\observed_vs_unobserved | 3 | 9 | 2 | 114.9 MB | 数据表 4、图 3、文档 1 | .csv 4、.png 3、.md 1 | 1 | 0 | 8 |
| D:\博士\Protein contour\Phospho\cellcycle_analysis | 3 | 8 | 8 | 1.1 MB | 数据表 4、图 2、脚本 2 | .csv 4、.py 2、.png 1 | 7 | 0 | 1 |
| D:\博士\protein_contour_wiki\paper | 2 | 7 | 2 | 6.0 KB | 数据表 5、文档 2 | .json 5、.canvas 1、.md 1 | 0 | 6 | 1 |
| D:\博士\Protein contour\Phospho\pilot_layerA | 3 | 7 | 6 | 4.4 MB | 数据表 5、脚本 1、其他 1 | .tsv 2、.csv 2、.py 1 | 5 | 1 | 1 |
| D:\博士\Protein contour\TCA_Succciny\Supplment data | 3 | 7 | 7 | 36.1 MB | 数据表 7 | .fasta 2、.tsv 2、.txt 1 | 0 | 4 | 3 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\confirm_ui | 3 | 6 | 6 | 30.7 KB | 数据表 5、其他 1 | .json 5、.log 1 | 0 | 1 | 5 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\notes | 3 | 6 | 6 | 16.5 KB | 文档 6 | .md 6 | 0 | 0 | 6 |
| D:\博士\protein_contour_wiki\wiki\gaps | 3 | 6 | 6 | 6.0 KB | 文档 6 | .md 6 | 0 | 0 | 6 |
| D:\博士\Protein contour\Phospho\Test | 3 | 5 | 5 | 244.2 KB | 脚本 5 | .r 4、.py 1 | 4 | 0 | 1 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\sources | 3 | 5 | 5 | 1.6 MB | 图 3、文档 1、数据表 1 | .png 3、.md 1、.txt 1 | 0 | 1 | 4 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\svg_final | 3 | 5 | 5 | 1.3 MB | 图 5 | .svg 5 | 0 | 5 | 0 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\svg_output | 3 | 5 | 5 | 16.7 KB | 图 5 | .svg 5 | 0 | 5 | 0 |
| D:\博士\protein_contour_wiki\.obsidian\plugins | 3 | 5 | 0 | 739.7 KB | 数据表 2、脚本 2、其他 1 | .json 2、.js 1、.sh 1 | 1 | 1 | 3 |
| D:\博士\protein_contour_wiki\paper\.obsidian | 3 | 5 | 5 | 5.8 KB | 数据表 5 | .json 5 | 0 | 5 | 0 |
| D:\博士\protein_contour_wiki\wiki\concepts | 3 | 5 | 5 | 9.6 KB | 文档 5 | .md 5 | 0 | 0 | 5 |
| D:\博士\protein_contour_wiki\wiki\figures | 3 | 5 | 1 | 1.2 MB | 文档 2、图 2、脚本 1 | .md 2、.png 1、.py 1 | 3 | 1 | 1 |
| D:\博士\Protein contour\Phospho\length_matched_temporal | 3 | 4 | 4 | 1.2 MB | 数据表 2、图 2 | .csv 2、.png 2 | 4 | 0 | 0 |
| D:\博士\Protein contour\Phospho\pilot_layerA_codex | 3 | 4 | 4 | 28.9 MB | 数据表 3、脚本 1 | .csv 2、.py 1、.json 1 | 3 | 1 | 0 |
| D:\博士\Protein contour\Phospho\three_analyses | 3 | 4 | 4 | 881.8 KB | 图 3、数据表 1 | .png 3、.csv 1 | 4 | 0 | 0 |
| D:\博士\Protein contour\Phospho\wincrashreport-x64 | 3 | 4 | 4 | 363.4 KB | 其他 2、数据表 1、文档 1 | .txt 1、.cfg 1、.chm 1 | 0 | 3 | 1 |
| D:\博士\protein_contour_wiki\.git\logs | 3 | 4 | 1 | 7.4 KB | 其他 4 | (无扩展名) 4 | 0 | 4 | 0 |
| D:\博士\protein_contour_wiki\wiki\hypotheses | 3 | 4 | 4 | 7.0 KB | 文档 4 | .md 4 | 0 | 0 | 4 |
| D:\博士\protein_contour_wiki\wiki\methods | 3 | 4 | 4 | 8.3 KB | 文档 4 | .md 4 | 0 | 0 | 4 |
| D:\博士\Protein contour\Phospho\anomaly_probe | 3 | 3 | 3 | 7.5 KB | 数据表 2、脚本 1 | .csv 2、.py 1 | 3 | 0 | 0 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\images | 3 | 3 | 3 | 1.6 MB | 图 3 | .png 3 | 0 | 0 | 3 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\validation | 3 | 3 | 3 | 9.4 KB | 数据表 3 | .json 3 | 0 | 0 | 3 |
| D:\博士\protein_contour_wiki\.git\refs | 3 | 3 | 0 | 112 B | 其他 3 | (无扩展名) 3 | 0 | 3 | 0 |
| D:\博士\Protein contour\Phospho\outputs | 3 | 2 | 2 | 480 B | 数据表 2 | .csv 2 | 0 | 0 | 2 |
| D:\博士\Protein contour\TCA_Succciny\Core data | 3 | 2 | 2 | 58.1 MB | 数据表 2 | .txt 2 | 0 | 2 | 0 |
| D:\博士\protein_contour_wiki\wiki\meetings | 3 | 2 | 2 | 4.8 KB | 文档 2 | .md 2 | 0 | 0 | 2 |
| D:\博士\Protein contour\Phospho\.claude | 3 | 1 | 0 | 8.0 KB | 文档 1 | .md 1 | 0 | 1 | 0 |
| D:\博士\Protein contour\Phospho\figures | 3 | 1 | 1 | 124.3 KB | 图 1 | .png 1 | 0 | 0 | 1 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\analysis | 3 | 1 | 1 | 529 B | 数据表 1 | .csv 1 | 0 | 0 | 1 |
| D:\博士\Protein contour\meeting_2026-08-03_pptmaster_ppt169_20260803\live_preview | 3 | 1 | 1 | 56.5 KB | 其他 1 | .log 1 | 0 | 1 | 0 |
| D:\博士\protein_contour_wiki\.git\info | 3 | 1 | 1 | 291 B | 其他 1 | (无扩展名) 1 | 0 | 1 | 0 |
| D:\博士\protein_contour_wiki\raw\meetings | 3 | 1 | 1 | 1.0 MB | 文档 1 | .pdf 1 | 0 | 0 | 1 |
| D:\博士\protein_contour_wiki\wiki\datasets | 3 | 1 | 1 | 3.7 KB | 文档 1 | .md 1 | 0 | 0 | 1 |
| D:\博士\protein_contour_wiki\wiki\proteins | 3 | 1 | 1 | 1.7 KB | 文档 1 | .md 1 | 0 | 0 | 1 |
