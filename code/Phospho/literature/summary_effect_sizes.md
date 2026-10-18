# 文献效应量汇总：磷酸化位点 / 蛋白的无序度差异

配套文件：`refs_effect_sizes.tsv`（27 行，16 列，已人工筛选，仅保留有具体数字的记录）；
`raw_search/topic1-4.tsv` 为四路检索的原始未筛选记录（含被剔除的"仅计数""无对照组"条目）。

## 300 字总结：我们的 δ ≈ 0.08（位点层）和 δ ≈ 0.2（蛋白层）在文献中的位置

文献几乎不报告标准化效应量，多为两组比例或 p 值。二分变量的 Cliff's δ 恰等于比例差，据此可统一换算。
位点层：Bludau 2022 用 AlphaFold 比较 PSP 调控位点与其余磷酸化位点，IDR 内比例 68.8% 对 62.6%，δ ≈ 0.06，OR 1.32；
结构区内暴露度两组几乎相同（72.6% 对 72.4%）。Tyanova 2013 只给出 p 值（KS p = 6.6e-13），且 90% 位点本就在无序区。
我们的 0.08 与这一"磷酸化位点对磷酸化位点"的量级一致：一旦对照组也是磷酸化位点，序列层无序度几乎不再提供区分力。
蛋白层：Valverde 2023 动态对非动态磷酸化蛋白的无序残基比例高 8.5 个百分点，由 W 统计量换算 δ ≈ 0.21；
信号蛋白对蛋白组（Iakoucheva 2002，66% 对 47%）δ ≈ 0.19；hub 对 end（Haynes 2006，78% 对 58%）δ ≈ 0.20。
我们的 0.2 正落在这一"蛋白类别"对比的典型区间（0.17–0.28）；激酶定义的类别（CDK1 底物，+21 个百分点）才明显更大。
空间蛋白质组学论文（Itzhak、Martinez-Val、Orre、Crook）均未报告重定位蛋白的无序度效应量，只报告计数（约 3% 蛋白移动）。
注意 Zhao 2021 的反向结果：多定位蛋白平均无序度更低（0.07 对 0.12）。

## 换算说明

- 二分变量（如"是否位于 IDR"）：Cliff's δ = p1 − p2，AUC = (δ+1)/2。我们的 δ=0.08 对应 AUC 0.54，δ=0.2 对应 AUC 0.60。
- Valverde 2023：R `wilcox.test` 的 W 即 Mann–Whitney U（含 tie 的 0.5 计数），AUC = W/(n1·n2)，δ = 2·AUC − 1。
  组大小 646（动态）/ 1,197（非动态）由作者仓库中的 `all_phosphosites_xenopus_nr.tab` 与 `ids/xenopuscluster[A-D]_mpi.txt` 复核。
- Bludau 2022 的"调控 vs 非调控"直接 OR（1.32）由其 Fisher 表中四格数自行计算（原文只给各组相对未修饰残基的 OR）。

## 来源级别（务必注意）

本环境的网络策略拒绝了全部文献全文站点（PMC、Europe PMC、Nature、Cell、Science、PLOS、eLife、bioRxiv、Crossref、OpenAlex 等），
只有搜索引擎摘要和 GitHub 可达。因此 `provenance` 列分三级：

- **A**：作者公开的 GitHub 分析仓库中的计算结果，已在本地克隆中逐行核对（Bludau 2022 `MannLabs/structuremap_analysis`；Valverde 2023 `gerodubra/CDK-mediated-phosphorylation-of-IDRs`）。数字可靠，但未与论文正文逐一对照图号。
- **B**：论文摘要原文（搜索片段）。
- **C**：正文/图注片段（搜索片段），数字可靠性次之，图表号部分未确认。

每行保留 `verbatim_quote` 便于复核。建议在开放网络后优先核对：Tyanova 2013 Fig 2（PMC3542066，取中位数/AUC）、
Ochoa 2020 Fig 2b 各特征 AUC、Xiao 2016 Sci Rep（PMC4835696）、Stenström 2020 Fig 5、Itzhak 2016 Fig 5/Supp 9、Martinez-Val 2021 Fig 4B。

## 检索到但无可用数字的文献（节选，详见 raw_search）

Olsen 2006/2010；Sharma 2014；Humphrey 2013/2015；Bekker-Jensen 2020；Francavilla 2016；Kanshin 2015；Leutert 2023（定性："regulated phosphosites enriched in disordered regions"）；
Studer 2016；Landry 2009；Strumillo 2019（仅样本量）；Itzhak 2016；Martinez-Val 2021；Geladaki 2019；Christopher 2021；Mendik 2019（仅计数）；
Kim/Gerstein 2008；Dosztányi 2006；Peng 2015；Xie 2007（仅关键词计数）。

## 我们这边的 n（site_meta.tsv）

11,826 条位点记录，10,873 个唯一位点，3,768 个蛋白（multiplicity 1/2/3 = 9,715 / 2,019 / 92）。
