# 05 四项定义核查（只读事实记录）

范围：仅 `code/` 下的脚本，未改任何文件。行号为仓库内文件行号。查不到的写"未找到"。
重复副本：`protein_contour_wiki/raw/code/*.py` 与 `Protein contour/Phospho/` 下同名文件 md5 一致（three_analyses、length_matched_temporal_analysis、task6_s5c_detection），下文只引后者。
`code/Phospho/0210trying.ipynb`（22115 行）与 `code/Protein contour/Phospho/0210trying.ipynb`（14549 行）内容不同，分别标注。

---

## 1. DynaMine 方向（biophys_json / b2bTools 读取方式）

### 1.1 总结

- **原始字段名**：所有脚本读的都是 b2bTools JSON 里 `data['residues'][i]` 下的键 `backbone`（同级键：`sidechain`、`disoMine`、`helix`、`sheet`、`coil`、`earlyFolding`、`ppII`、`aa`、`seqpos`）。
- **重命名**：test.py 家族把 `backbone` 重命名为 `backbone_dynamics`；`run_layerA_full.py` 同样重命名；`run_scaffold_check*.py` 重命名为 `Backbone`/`avg_Backbone`；其余脚本保留 `backbone`。
- **数值变换**：未找到任何取反、`1 - x`、乘 `-1` 或分箱。出现的变换只有：(a) 位点 ±5 窗口均值；(b) 蛋白级 S/T/Y 残基均值；(c) PCA 前 `StandardScaler` z-score（0331.py）；(d) stage6 中按匹配用的 z 标准化。
- **"刚性/柔性"解释原文**：仅在 test.py 家族的报告模板里出现一句 `范围: 0-1 (0=刚性, 1=高度灵活)`，以及 0331.py / 0228.py PCA 解读段落中的 `"dynamic/flexible vs structured" 轴`、`higher dynamics`。其余脚本无方向性注释。

### 1.2 逐脚本

**(A) test.py 家族**（`Test/test.py`、`test.py`、`0220.py`、`0228.py`、`0331.py`；提取函数体 md5 均为 fb536dfc，完全相同）

读取与重命名 — `code/Protein contour/Phospho/Test/test.py:166-178`（同一函数在 `test.py:260`、`0220.py:260`、`0228.py:260`、`0331.py:296` 起）：
```python
    residue_data = residues[position_idx]
    # 步骤7: 提取我们需要的生物物理特征
    features = {
        'backbone_dynamics': residue_data.get('backbone', None),
        'sidechain_dynamics': residue_data.get('sidechain', None),
        'disorder_propensity': residue_data.get('disoMine', None),
        'helix_propensity': residue_data.get('helix', None),
        'sheet_propensity': residue_data.get('sheet', None),
        'coil_propensity': residue_data.get('coil', None),
        'earlyFolding': residue_data.get('earlyFolding', None),
```
另两处窗口提取器同样用 `residue_data.get('backbone', None)` / `np.nan`：`Test/test.py:3456`、`Test/test.py:3586`；`0331.py:3597`、`0331.py:3727`。

解释原文 — `code/Protein contour/Phospho/Test/test.py:534-536`（同文：`test.py:673-675`、`0220.py:673-675`、`0228.py:673-675`、`0331.py:709-711`）：
```
1. Backbone Dynamics (骨架动力学)
   - 描述: 蛋白质主链的灵活性
   - 范围: 0-1 (0=刚性, 1=高度灵活)
```

PCA 前标准化 — `code/Phospho/0331.py:10065-10067`（另两处：`0331.py:10381`、`0331.py:10477`）：
```python
X = df[FEAT_SIG].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
```

PCA 解读原文 — `code/Phospho/0331.py:10288-10295`（同文 `0228.py:9821-9828`）：
```
1. PC1 ({ev[0]*100:.0f}%) 由 backbone/side-chain dynamics (-) 和 coil (+) 驱动
   → 这是 "dynamic/flexible vs structured" 轴
2. PC2 ({ev[1]*100:.0f}%) 几乎完全由 disorder (+0.9) 驱动
   → 这是 "ordered vs disordered" 轴
3. 三组的 centroid 偏移方向:
   → High-mob phospho: 偏向 higher disorder, higher dynamics
   → Low-mob phospho:  偏向 lower disorder, lower dynamics
```

原值直接绘图（无变换）：`Test/test.py:1553`、`:1614`、`:2158`、`:2519`，例 `ax.plot(residue_df['seqpos'], residue_df['backbone'], label='Backbone dynamics',`。

**(B) 0210trying.ipynb**（两份）
- `code/Phospho/0210trying.ipynb:612`、`:5637`、`:5830`：`'backbone_dynamics': residue_data.get('backbone', None/np.nan)`；解释原文 `:1187` `范围: 0-1 (0=刚性, 1=高度灵活)`；PCA 解读 `:14487`、`:14820`。
- `code/Protein contour/Phospho/0210trying.ipynb:531`、`:5611`、`:5804`；解释 `:1106`；PCA 解读 `:13847`、`:14180`。
- 无其他变换。

**(C) three_analyses.py**（`code/Protein contour/Phospho/three_analyses.py`）
- 字段：`:24` `FEATURES = ['backbone', 'sidechain', 'disoMine', 'helix', 'sheet', 'coil', 'earlyFolding']`；标签 `:26` `'backbone': 'Backbone dynamics'`。
- 变换：蛋白级非磷酸化 S/T/Y 残基均值，`:90-96`：
```python
    for r in data['residues']:
        if r['aa'] in ('S', 'T', 'Y') and r['seqpos'] not in phos_pos:
            sty_residues.append(r)
    ...
    means = {f: np.mean([r[f] for r in sty_residues if f in r]) for f in FEATURES}
```
- 刚性/柔性注释：未找到。

**(D) length_matched_temporal_analysis.py**（`code/Protein contour/Phospho/length_matched_temporal_analysis.py`）
- 字段 `:22`，标签 `:24` `'backbone': 'Backbone dynamics'`。均值 `:62-68`：
```python
    sty = [r for r in data['residues'] if r['aa'] in ('S', 'T', 'Y')]
    ...
    for f in FEATURES:
        vals = [r[f] for r in sty if f in r]
        rec[f] = np.mean(vals) if vals else np.nan
```
- 刚性/柔性注释：未找到。

**(E) Zhihan/Test/13_breadth_features.py**
- `:17` `CH = ["backbone","sidechain","ppII","coil","sheet","helix","earlyFolding","disoMine"]`
- `:153-156` 单点值 + ±5 窗口均值，无其他变换：
```python
            row = {f"single_{c}": float(rr[c]) for c in CH}
            for c in CH:
                vv = [float(pm[p][c]) for p in range(pos-5, pos+6) if p in pm]
                row[f"win5_{c}"] = sum(vv)/len(vv)
```
- 刚性/柔性注释：未找到。

**(F) Zhihan/Test/09_feature_comparison.py** — `:18` CHANNELS 同上；`:215-220` single_/win5_ 同上。注释：未找到。
**(G) Zhihan/Test/09b_paired_patch.py** — `:10-11`；`:96-99` 同上。注释：未找到。

**(H) scope3p data/observed_vs_unobserved/run_observed_vs_unobserved.py**
- `:19` `FEATURES = ["backbone", "sidechain", "disoMine", "helix", "sheet", "coil", "earlyFolding"]`
- `:171` `values = tuple(finite_float(record.get(feature)) for feature in FEATURES)`；无变换。
- `:501` 文字："The seven biophysical values are **predictions, not measurements**."；刚性/柔性：未找到。

**(I) Phospho/pilot_layerA/run_pilot.py** — `:14` FEATS 'backbone'…；`:113-115`：
```python
        ok=all(rr.get(f) is not None for f in FEATS)
        if not ok: continue
        for f in FEATS: row[f]=float(rr[f])
```
**(J) Phospho/layerA_full/run_layerA_full.py** — 重命名表 `:34-42`：
```python
FEATURES = [
    ("backbone", "backbone_dynamics"),
    ("sidechain", "sidechain_dynamics"),
    ("disoMine", "disorder_propensity"),
```
写入 `:482-483` `for src, label in FEATURES: row[label] = float(rr[src])`。无其他变换；注释未找到。

**(K) Phospho/pilot_layerA_codex/run_layerA_pilot_codex.py** — `:21-29` FEATURES 'backbone'…；`:159-160` `row[feature] = float(rr[feature])`。

**(L) Phospho/martinez_network_check/stage6_biophysics.py**
- `:10` FEATS；`:74` `arr={f:np.array([r.get(f,np.nan) for r in res],dtype=float) for f in FEATS}`
- 全长均值 `:79`、位点 ±5 局部均值 `:84-85`；匹配时 z 标准化 `:127-128`：
```python
    mu,sd=B[c].mean(),B[c].std()
    movb[c+'_z']=(movb[c]-mu)/sd; nonb[c+'_z']=(nonb[c]-mu)/sd
```
**(M) Phospho/martinez_network_check/step2345_presence_biophys.py** — 不从 JSON 读 backbone；用已提取 CSV 列 `:18-19` `FEATS_IN_DF = ['backbone_dynamics','sidechain_dynamics',...]`，JSON 只补 earlyFolding `:38`。

**(N) Phospho/scaffold_check/run_scaffold_check.py** — `:21` `JKEYS = {'backbone':'Backbone',...}`；蛋白级均值 `:35-37`：
```python
    for jk,lab in JKEYS.items():
        v=[r.get(jk) for r in res if r.get(jk) is not None]
        e['avg_'+lab]=np.mean(v) if v else np.nan
```
**(O) Phospho/scaffold_check/run_scaffold_check_v2.py** — `:40-48` `FEATURE_KEYS = {"backbone": "Backbone", ...}`；`:183-185` 同上均值。

**(P) Phospho/anomaly_probe/probe.py** — `:104` `vals=[r.get('backbone') for r in res if r.get('backbone') is not None]`，`:105` 均值；`:72` 硬编码数值 `feat_collapse={'Backbone':(-0.293,-0.160),...}`（来源未在脚本内说明）。

**(Q) Phospho/anomaly_first_codex/run_anomaly_scan.py** — 不读 JSON；`:114` `"backbone_dynamics": "backbone_dynamics"` 仅名字映射，输入为 `stage6_biophysics_features.csv`（`:70`）。

**(R) martinez_network_check/kl_core.py `:50-51`、test_window.py `:14-15`** — 只读 JSON 的 `aa`/`seqpos` 拼序列，不读 backbone。

**(S) protein_contour_wiki/wiki/figures/fig01_overview/fig01_overview.py** — 仅图注文字：`:49` `label='b2bTools\n7 per-residue\nfeatures'`；`:77` `label='F004  [V]\nPhospho proteins\nLOWER backbone-dyn,\n+disoMine\n(7/7 p<0.001)'`。

### 1.3 R 脚本
`--include=*.R` 内未找到 `backbone`/`biophys_json`/`b2bTools`。

---

## 2. DIA-NN 参数

### 2.1 仓库内唯一的 diann 命令：`code/Protein contour/Zhihan/Test/stage1_lib_check.sh`
用途（`:3-5`）："阶段 1: 小 FASTA 建库 + 验证 UniMod:21 是否进库"，即**建库测试**，不是正式搜库。二进制 `:14` `DIANN="$VSC_DATA/software/diann/diann-2.2.0/diann-linux"`。

完整命令 `:42-65`：
```bash
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
```
（续 `:52-65`）`--max-pr-charge 4`、`--min-pr-mz 300`、`--max-pr-mz 1800`、`--min-fr-mz 200`、`--max-fr-mz 1800`、`--fixed-mod "UniMod:4,57.021464,C"`、`--var-mod "UniMod:35,15.994915,M"`、`--var-mod "UniMod:1,42.010565,*n"`、`--var-mod "UniMod:21,79.966331,STY"`、`--var-mod-max 3`、`--threads "$THREADS"`、`--out-lib "$WORKDIR/output/stage1_lib.parquet"`、`--out "$WORKDIR/output/stage1_lib_gen.tsv"`。

第二条（库格式转换）`:93`：
```bash
    "$DIANN" --lib "$LIB_FILE" --out-lib "$WORKDIR/output/stage1_lib.parquet" --threads "$THREADS" \
```
注释 `:39`："--var-mod-max 保持 3, 与生产配置一致 —— 阶段 1 验证的就是这套配置本身"。

### 2.2 逐项
| 项 | 结果 | 位置 |
|---|---|---|
| `--reanalyse`（MBR） | 未找到 | 全仓库 0 次 |
| `--matrices` | 未找到（"matrices" 3 次均为 numpy/注释，无关） | — |
| 归一化相关 flag（`--no-norm` 等） | 未找到 | — |
| `--qvalue` | 未找到；下游 Python 自己过滤 `Q.Value<=0.01` | `07_main_analysis.py:21`、`:149`；`05_build_occupancy.py:7` |
| PTM 定位/打分 flag | 未找到；下游用 `PTM.Site.Confidence>=0.75` | `07_main_analysis.py:20`、`:155`；`validate_run.py:15` |
| `--fasta-search` | 有 | `stage1_lib_check.sh:44` |
| `--predictor` | 有 | `:45` |
| `--gen-spec-lib` / `--out-lib` | 有 | `:46`、`:63` |
| `--lib` | 仅用于 speclib→parquet 转换 | `:93` |
| 质量精度（`--mass-acc`、`--mass-acc-ms1`） | 未找到 | — |
| 修饰 | `--fixed-mod UniMod:4`；`--var-mod UniMod:35 / UniMod:1 / UniMod:21`；`--var-mod-max 3` | `:57-61` |
| 酶切/长度/电荷/mz | `--cut "K*,R*,!*P"`、missed 2、pep 7-30、charge 2-4、pr mz 300-1800、fr mz 200-1800 | `:47-56` |

### 2.3 正式搜库命令不在仓库
`list_pride_files.sh:6` 引用 `sbatch stage2_full_search.slurm <RAW_URL>`，但仓库内无任何 `.slurm` 文件，`stage2_full_search` 未找到。

### 2.4 引用正式搜库输出（`$VSC_DATA/rerun_2026-08/{phospho,proteome}/*/report.parquet`）的脚本
全部在 `code/Protein contour/Zhihan/Test/`：
- `04_qc_summary.py`（`report.stats.tsv` `:28-29`；`report.parquet` `:49`）
- `validate_run.py`（`:2`；`:8` "phosphosites_90/_99 矩阵存在性检查"）
- `validate_unimod.py`（`:2` 谱库 parquet/主报告 parquet/矩阵 tsv）
- `05_build_occupancy.py`（`:6`、`:90`）
- `files1/06_occupancy_variants.py`（`:76`）
- `07_main_analysis.py`（`:84`）
- `11_martinez_counts.py`（`:7`、`:58`）
- `11b_dump_site_membership.py`（`:7`、`:59`）
- `13_breadth_features.py`（`:102`）
- `15_fpr_null.py`（`:76`）

使用的列 — `07_main_analysis.py:92-97`：
```python
    c["mod"] = pick_col(cols, ["Modified.Sequence"], ...)
    c["strip"] = pick_col(cols, ["Stripped.Sequence"], ...)
    c["prot"] = pick_col(cols, ["Protein.Ids"], lambda x: x == "Protein.Group")
    c["qty"] = pick_col(cols, ["Precursor.Quantity"], lambda x: False)
    c["qval"] = pick_col(cols, ["Q.Value"], lambda x: x.endswith("Q.Value") and all(t not in x for t in ("PG", "PTM", "Lib", "Global")))
    c["conf"] = pick_col(cols, ["PTM.Site.Confidence"], lambda x: "Site.Confidence" in x) if need_conf else None
```
强度用 `Precursor.Quantity`（非 `Precursor.Normalised`），07 自行做 per-run log2 中位数中心化 `:238-247`。

---

## 3. 08 的 profile / profile_shift 定义

### 3.1 08 不计算轮廓，只读 07 的产物
`code/Protein contour/Zhihan/Test/08_q1q3_cross.py:4-5`："输入 (07 产物, 不重扫 parquet): analysis07/site_traj.tsv, analysis07/site_profiles.tsv"。

读取与阈值 `:96-100`、`:108`、`:124`：
```python
            if r["phos_shift_vs_ctrl"] != "":
                v = float(r["phos_shift_vs_ctrl"])
                all_shift.append(v)
                if d["max_shift"] is None or v > d["max_shift"]:
                    d["max_shift"], d["shift_tp"] = v, r["timepoint"]
...
    thr_shift = pctl(all_shift, PCT)
...
        mover = p.get("max_shift") is not None and thr_shift is not None and p["max_shift"] >= thr_shift
```
PCT 默认 90（`:45`）。输出列名 `max_profile_shift`（`:154`），汇总行 `profile_shift P{PCT}`（`:160`）。即 08 的 "profile_shift" = 07 的 `phos_shift_vs_ctrl`。

### 3.2 公式在 07：位点自身 6 组分组成的 CLR / Aitchison 距离（不是相对蛋白）
`code/Protein contour/Zhihan/Test/07_main_analysis.py:103-109`：
```python
def clr(props):
    logs = [math.log(p) for p in props]
    m = sum(logs) / len(logs)
    return [x - m for x in logs]

def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
```
`:317-321`（6 组分向量：每组分取重复的中位数，缺失填 0，加伪计数，**除以自身总和**，再 CLR）：
```python
    def clr_profile(frvals, pseudo):
        vec = [median(frvals[fr]) if frvals.get(fr) else 0.0 for fr in FRS]
        vec = [v + pseudo for v in vec]
        s = sum(vec)
        return clr([v / s for v in vec])
```
伪计数 `:323-324` `pseudo = 0.5 * min(pos_all)`。
位移 `:344`、`:351-352`：
```python
            clr_ctrl = clr_profile(ctrl_s, pseudo)
...
                clr_tp = clr_profile(s, pseudo)
                shift = dist(clr_tp, clr_ctrl)
```
→ `phos_shift_vs_ctrl` = 该位点 tp 时的 6 组分 CLR 轮廓 与 其 CTRL 轮廓 的欧氏距离（Aitchison 距离）。输入强度 `num` 是 per-run 中位数中心化后的磷酸肽前体强度和（`:245`）。

相对蛋白的量是另外两列，`:355-359`：
```python
                if clr_pc and (acc, tp) in prof_prot:
                    clr_pt = clr_profile(prof_prot[(acc, tp)], pseudo_q)
                    pshift = round(dist(clr_pt, clr_pc), 4)
                    div = round(dist(clr_tp, clr_pt), 4)
                    ddelta = round(dist(clr_tp, clr_pt) - div_ctrl, 4)
```
对应列 `prot_shift_vs_ctrl`、`site_vs_prot_divergence`、`divergence_delta_vs_ctrl`（`:333-335`）。08 把 `divergence_delta_vs_ctrl` 的绝对值单独作为 `divergence_mover`（`08:11`、`:101-105`、`:125`），与 profile_shift 分开。

### 3.3 07 的 intensity-responsive / occupancy-responsive 判定原文
常量 `:20-23`：
```python
LOC_CUTOFF = 0.75
QVAL_CUTOFF = 0.01
FC_THRESH = 1.0
MIN_REPS = 2
```
两个量的定义 `:266-272`（inten = log2 磷酸肽强度；occ = log2(磷酸肽强度 / proteome 分母)）：
```python
            for tp, rep in obs_use:
                key = (tp, fr, rep)
                p = num[(sid, key)]
                inten[tp].append(math.log2(p))
                d = den_res.get((sid, key)) if tier == "residue" else den_prot.get((acc, key))
                if d:
                    occ[tp].append(math.log2(p / d))
```
判定 `:274-290`：
```python
            def responder(vals):
                if len(vals.get("CTRL", [])) < MIN_REPS:
                    return False, ""
                base = median(vals["CTRL"])
                for tp in EGF_TPS:
                    v = vals.get(tp, [])
                    if len(v) < MIN_REPS:
                        continue
                    fc = median(v) - base
                    if abs(fc) >= FC_THRESH and all((x - base) * fc > 0 for x in v):
                        return True, tp
                return False, ""
            o_resp, o_tp = responder(occ)
            i_resp, i_tp = responder(inten)
            cls = ("both" if o_resp and i_resp else
                   "occupancy_only" if o_resp else
                   "intensity_only" if i_resp else "none")
```
docstring `:12`："筛选标准: |log2FC|>=1, 每侧 >=2 个重复, 且方向一致 —— 明确为候选筛选, 非最终统计检验。"
单位：位点 × 组分（fraction）轨迹（`:250-252`）。08 把它聚合到位点级（任一组分达标即算）`08:60-64`，类名重映射 `08:18-19`。

---

## 4. 线 B 的 mobility 定义

### 4.1 主定义（Movement Score，Cell 20/21），五代 .py 完全相同
Cell 21 代码块 md5 = a352a89d（`Test/test.py:4264-4377`、`test.py:4437-4550`、`0220.py:4437-4550`、`0228.py:4437-4550`、`0331.py:4439-4552`）；`assign_mobility_level` md5 = f3ae5fc8；Cell 20 读表块 md5 = e94dc34f，0331.py 仅路径不同（`D:\博士\Phospho\Full_data.xlsx` → `D:\博士\Protein contour\Phospho\Full_data.xlsx`）。两份 0210trying.ipynb 中函数文本相同（`code/Phospho/0210trying.ipynb:7116-7170`；`code/Protein contour/Phospho/0210trying.ipynb:7090`）。以下引 `code/Protein contour/Phospho/Test/test.py`。

**输入表**：S5D 位点表（log2），不是 S5A 蛋白表。`:4136-4140`：
```python
df_intensity = pd.read_excel(
    r"D:\博士\Phospho\Full_data.xlsx",
    sheet_name="S5D-Log2 proc HeLa+EGF PHOS",
    engine='openpyxl'
)
```
**组分聚合**：按 FR1+2=Cyt、FR3+4=Mem、FR5+6=Nuc（`:4150-4157`），对 2 组分 × 4 重复的 **log2 值直接取平均**，无 2^ 回线性，`:4177-4187`：
```python
        for fr in fractions:
            for rep in range(1, 5):  # Rep1-Rep4
                col_name = f"EGF_{time_point}_{fr}_Rep{rep}"
                if col_name in df_intensity.columns:
                    all_cols.append(col_name)
        if len(all_cols) > 0:
            intensity_data[f'{time_label}_{compartment}_intensity'] = df_intensity[all_cols].mean(axis=1)
```
**百分比**（位点自身三区室之和归一）`:4271-4287`：
```python
    success_df[total_col] = (
        success_df[f'{time_label}_Cyt_intensity'] +
        success_df[f'{time_label}_Mem_intensity'] +
        success_df[f'{time_label}_Nuc_intensity']
    )
        success_df[pct_col] = np.where(
            success_df[total_col] > 0,
            success_df[intensity_col] / success_df[total_col],
            0.0
        )
```
**公式**（相对 0min，每个后续时间点取三区室 |Δpct| 中最大两项之和 ×100，再取各时间点最大值）`:4330-4346`：
```python
    for time_label in ['2min', '8min', '20min', '90min']:
        if row[f'{time_label}_total_intensity'] == 0:
            continue
        deltas = {}
        for comp in compartments:
            delta = abs(row[f'{time_label}_{comp}_pct'] - baseline[comp])
            deltas[comp] = delta
        sorted_deltas = sorted(deltas.values(), reverse=True)
        mobility = (sorted_deltas[0] + sorted_deltas[1]) * 100  # 转成百分比
        if mobility > max_mobility:
```
0min 总强度为 0 → `movement_score = NaN`（`:4307-4313`）。单位：按 `PTM_collapse_key` 的**磷酸肽（位点）轮廓**，不是蛋白轮廓。

**阈值**：High = `movement_score >= 10`，Low = `< 5`（`:4675`、`:4678`；全局常量 `test.py:90` `MOBILITY_THRESHOLD = 10`，`0228.py:8446` / `0331.py:8448` `LOW_MOB_THRESHOLD  = 5`）。三档 `:4584-4592`：
```python
def assign_mobility_level(score):
    if pd.isna(score):
        return 'Undetectable'
    elif score >= 10:
        return 'High'
    elif score >= 5:
        return 'Medium'
    else:
        return 'Low'
```
细分档（Cell 26a）`Test/test.py:5886-5892`：≥20 very_high、15–20 high、10–15 moderate_high、2–5 low、<2 very_low（`0331.py:6275-6281` 同）。

### 4.2 第二份实现（"全自动数据构建与修复"，重建 success_df 时重算），公式一致
`code/Phospho/0331.py:6099`（同块 `test.py:6172`、`0220.py:6172`、`0228.py:6172`）读同一 S5D 表；区室均值仍为 log2 均值 `:6109-6115`；公式 `:6195-6200`：
```python
        curr_pct = {c: row[f'{tl}_{c}_intensity']/total_t for c in ['Cyt', 'Mem', 'Nuc']}
        # 计算差异 (sum of absolute differences / 2 就是移动比例，或者取最大的两个变动之和)
        # 这里沿用您之前的逻辑：取变动最大的两个部分之和
        deltas = sorted([abs(curr_pct[c] - base_pct[c]) for c in ['Cyt', 'Mem', 'Nuc']], reverse=True)
        mobility = (deltas[0] + deltas[1]) * 100
```
与 4.1 同公式，不输出方向。

### 4.3 第三种定义（基于 0/1 出现次数，"Cell 1: 数据准备与STY分组"），与前两者不同
`code/Protein contour/Phospho/Test/test.py:4958-4976`（同块 `test.py:5158`、`0220.py:5158`、`0228.py:5158`、`0331.py:5160`、`code/Phospho/0210trying.ipynb:8205`）：
```python
compartment_cols = [
    '0min_Cyt', '0min_Mem', '0min_Nuc',
    ...
    '90min_Cyt', '90min_Mem', '90min_Nuc'
]
# 计算每个位点的移动性得分 (出现的区室数量)
def calculate_mobility_score(row):
    """计算位点在不同区室出现的次数"""
    return row[compartment_cols].sum()
final_df['mobility_score'] = final_df.apply(calculate_mobility_score, axis=1)
mobility_threshold = final_df['mobility_score'].quantile(0.90)  # 可调整
```
输入是 15 个 0/1 出现列（来自位点表，见 4.5），阈值 = 90 百分位，`>=` 为 High（`:4983-4985`）。注释 `:4955-4956`："根据对话记录,高移动性定义为:在不同时间点和区室之间有转移的位点 / 这里使用ModScore或者根据0/1模式判断"。

### 4.4 前代 R 实现（不在要求的链内，但同名 "Movement Score"，定义不同）
`code/Protein contour/Phospho/Test/Phos_Functions.R`：输入 S5D（`Phos_Main_Analysis.R:13`）；**先 2^log2 回线性再按组分平均**，`:58-59`、`:71`、`:75-85`：
```r
  mat_log2 <- as.matrix(sub_df[ , -1, drop = FALSE])
  mat_lin  <- 2 ^ mat_log2
  fr_means <- t(apply(mat_lin, 1, function(x) tapply(x, frac_labels, mean, na.rm = TRUE)))
  cyto    <- fr_means$FR1 + fr_means$FR2
  memb    <- fr_means$FR3 + fr_means$FR4
  nucleus <- fr_means$FR5 + fr_means$FR6
  total   <- cyto + memb + nucleus
    Cyto  = cyto / total,
```
分数 = 三区室 |Δ| 之和（非最大两项），`:103-107`：
```r
  d_cyto <- df[[cyto_stim]] - df[[cyto_ctrl]]
  d_memb <- df[[memb_stim]] - df[[memb_ctrl]]
  d_nuc  <- df[[nuc_stim]]  - df[[nuc_ctrl]]
  mv_score <- abs(d_cyto) + abs(d_memb) + abs(d_nuc)
```
阈值 `Phos_Main_Analysis.R:70-73`：`filter(Movement_Score >= 0.1, q_value < 0.05)`（分数为 0–2 的比例尺，0.1 对应 Python 版的 10%）。与 Python 版的差别：线性 vs log2 平均；三项和 vs 最大两项和；每时间点一行 vs 取最大。

### 4.5 presence pattern（1/2/3 区室分组）
**输入**：`Phospho_Site_Compartment_Table.csv`（`Test/test.py:12` `SITES_TABLE = r"D:\博士\Phospho\Phospho_Site_Compartment_Table.csv"`；`:100` `sites_df = pd.read_csv(SITES_TABLE)`）。其列（notebook 输出 `code/Phospho/0210trying.ipynb:6676`）：`PTM_collapse_key, Gene, Site, AA, Position, ModNumber, 0min_Cyt, 0min_Mem, 0min_Nuc, 2min_Cyt, ..., 90min_Nuc, ...`。**生成该 CSV 的脚本：未找到**（仓库内只有消费者）。

相关但非来源：`task6_s5c_detection.py:27-28` 读 `S5C-Collapse PHOS Hela+EGF`，`:76-84` 用"8 列中非 NaN ≥2 → 1"重建二值矩阵，用于与 success_df 比对（`:2` "S5C raw data vs imputed success_df"）：
```python
# For each site, each (tp, comp): count non-NaN values. If >= 2 -> detected (1), else 0
            non_nan_count = s5c[cols].notna().sum(axis=1)
            binary_results[col_name] = (non_nan_count >= 2).astype(int)
```

**定义 A（时间折叠，1/2/3 区室）** — `Test/test.py:1219-1221`（`0331.py:102-104` 同）：
```python
    success_df['Cytany'] = success_df[['0min_Cyt', '2min_Cyt', '8min_Cyt', '20min_Cyt', '90min_Cyt']].max(axis=1)
    success_df['Memany'] = success_df[['0min_Mem', '2min_Mem', '8min_Mem', '20min_Mem', '90min_Mem']].max(axis=1)
    success_df['Nucany'] = success_df[['0min_Nuc', '2min_Nuc', '8min_Nuc', '20min_Nuc', '90min_Nuc']].max(axis=1)
```
分组 `Test/test.py:662-684`：
```python
def assign_group_A(row):
    c = row['Cyt_any']; m = row['Mem_any']; n = row['Nuc_any']
    if c == 1 and m == 0 and n == 0:   return "C-only"
    elif c == 0 and m == 1 and n == 0: return "M-only"
    elif c == 0 and m == 0 and n == 1: return "N-only"
    elif c == 1 and m == 1 and n == 0: return "C&M"
    elif c == 1 and m == 0 and n == 1: return "C&N"
    elif c == 0 and m == 1 and n == 1: return "M&N"
    elif c == 1 and m == 1 and n == 1: return "C&M&N"
    else: return "None"
```
同一逻辑的另两处写法：`:1207-1215`（标签 CM/CN/MN/CMN）、`:904-930` `assign_localization_group`（用 `(row[cols] > 0).any()`，标签 CMN/CM/CN/MN/C-only/M-only/N-only/Unknown）。0331.py 只分析 ≥2 区室组：`0331.py:125` `analysis_groups = ['C&M', 'C&N', 'M&N', 'C&M&N']`。

**定义 C（逐时间点状态与轨迹）** — `Test/test.py:688-713` `get_state_at_time` 返回 C/M/N/CM/CN/MN/CMN/None；`:716-752` `assign_group_C` 由 5 个状态串 `trajectory` 给出 sustained-single / sustained-multi / early-only / late-appearing / transient / translocation / other：
```python
    if len(set(states)) == 1 and states[0] in ["C", "M", "N"]:
        return "sustained-single"
    if len(set(states)) == 1 and states[0] in ["CN", "CM", "MN", "CMN"]:
        return "sustained-multi"
    if len(non_none_states) >= 2 and all(s == "None" for s in states[2:]):
        return "early-only"
    if all(s == "None" for s in states[:2]) and len(non_none_states) >= 2:
        return "late-appearing"
    if len(set(non_none_states)) > 1:
        return "translocation"
```

---

## 附：未找到清单
- 正式搜库脚本 `stage2_full_search.slurm` 及任何 `.slurm`。
- DIA-NN flag：`--reanalyse`、`--matrices`、`--qvalue`、`--mass-acc*`、任何归一化或 PTM 打分 flag。
- `Phospho_Site_Compartment_Table.csv` 的生成脚本。
- 对 `backbone` 做取反 / `1 - x` / 分箱的代码。
- R 脚本中对 b2bTools/DynaMine 的读取。
