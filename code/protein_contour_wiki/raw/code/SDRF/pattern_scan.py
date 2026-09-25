"""
Phase 2 v3: 扫描所有 raw filename,生成模式报告
Output: out/pattern_report.md
"""
import json
import re
from collections import Counter
from pathlib import Path

OUT_DIR = Path("./out")
files = json.loads((OUT_DIR / "raw_files.json").read_text())

# ============ Token 分类 regex (v3) ============

CELL_PATTERNS = ["HeLa", "U2OS"]
ORGAN_PART_PATTERNS = ["Kidney", "Liver", "Muscle", "Brain", "Heart", "Lung", "Spleen"]

STIM_RE = re.compile(r"^(EGF|OSMSTRS{0,2}|OSM|TGF|IL\d*)$", re.I)
TIME_RE = re.compile(r"^(\d+(min|h|hour)|CTRL|control)$", re.I)
COMBINED_TIME_RE = re.compile(r"^(EGF|OSMSTRS{0,2}|OSM|TGF)(\d+(min|h|hour)|CTRL)$", re.I)
OSMSTR_TIME_RE = re.compile(r"^\d+h(SORB|REL|RELEASE).*$", re.I)
RELEASE_TOKEN_RE = re.compile(r"^\d+min(REL|RELEASE)$", re.I)  # 30minRELEASE 单独 token

FRAC_RE = re.compile(r"^(?:fraction|frac|FR|F)(\d+)$", re.I)
REP_RE = re.compile(r"^Rep\d+$", re.I)
TRAIL_DIGIT_RE = re.compile(r"^\d{2,3}$")
DATE_RE = re.compile(r"^\d{8}$")
WT_ID_RE = re.compile(r"^WT-(Stim-)?\d+[LR]?$")
GELBAND_COORD_RE = re.compile(r"^[A-J]\d{1,2}$")

MACHINE_PARAM_RES = [
    re.compile(r"^FAIMS(CV)?(-\d+)?$", re.I),  # FAIMS / FAIMS-45 / FAIMSCV-45
    re.compile(r"^\d+spd$", re.I),
    re.compile(r"^\d+m$"),
    re.compile(r"^\d+k\d+s$", re.I),
    re.compile(r"^EXPL\d+$", re.I),
    re.compile(r"^Evo\d+$", re.I),
    re.compile(r"^AMV(-CLK)?$", re.I),
]
MACHINE_PARAM_LITERAL = {"SA", "DAU", "DIA", "DDA"}

EXPDESIGN_RE = re.compile(r"^(SUBCELL|SubCell|subcell|WCL|GelBand)$", re.I)
MODTYPE_RE = re.compile(r"^(prot|phos|phospho|wcl|total)$", re.I)
GELBAND_MARKER_RE = re.compile(r"^GelBand$", re.I)
PERC_RE = re.compile(r"^\d+perc$", re.I)
TREATMENT_LITERAL = {"stimulated", "serumstarv"}


def is_machine_param(tok):
    if tok.upper() in MACHINE_PARAM_LITERAL:
        return True
    return any(r.match(tok) for r in MACHINE_PARAM_RES)


# ============ 拆 token ============

file_tokens = []
all_tokens = []
for f in files:
    name = f["fileName"].rsplit(".raw", 1)[0]
    tokens = name.split("_")
    file_tokens.append((f["fileName"], tokens))
    all_tokens.extend(tokens)

token_counts = Counter(all_tokens)

# ============ 启发式分类 ============

mouse_files = [fn for fn, toks in file_tokens
               if any(t.lower() in ("mouse", "mus") for t in toks)]
human_files = [fn for fn, _ in file_tokens if fn not in set(mouse_files)]

organ_part_files = {op: [] for op in ORGAN_PART_PATTERNS}
for fn, toks in file_tokens:
    if fn in mouse_files:
        for op in ORGAN_PART_PATTERNS:
            if any(op.lower() == t.lower() for t in toks):
                organ_part_files[op].append(fn)

cell_files = {cp: [] for cp in CELL_PATTERNS}
for fn, toks in file_tokens:
    if fn in mouse_files:
        continue
    for cp in CELL_PATTERNS:
        if any(cp.lower() == t.lower() for t in toks):
            cell_files[cp].append(fn)

phospho_files = [fn for fn, toks in file_tokens
                 if any("phos" in t.lower() for t in toks)]
total_files = [fn for fn, toks in file_tokens
               if any(t.lower() in ("prot", "wcl", "total") for t in toks)
               and fn not in phospho_files]
gelband_files = [fn for fn, toks in file_tokens
                 if any(GELBAND_MARKER_RE.match(t) for t in toks)]
mod_unclassified = (set(fn for fn, _ in file_tokens)
                    - set(phospho_files) - set(total_files) - set(gelband_files))

subcell_files = [fn for fn, toks in file_tokens
                 if any(t.lower() == "subcell" for t in toks)
                 and fn not in gelband_files]
wcl_files = [fn for fn, toks in file_tokens
             if any(t.lower() == "wcl" for t in toks)]


def parse_time(tokens):
    # 1) combined (EGF8min)
    for tok in tokens:
        if COMBINED_TIME_RE.match(tok):
            return tok
    # 2) STIM + TIME 分离
    for i, tok in enumerate(tokens):
        if STIM_RE.match(tok) and i + 1 < len(tokens):
            nxt = tokens[i + 1]
            if TIME_RE.match(nxt) or OSMSTR_TIME_RE.match(nxt) or RELEASE_TOKEN_RE.match(nxt):
                return f"{tok}_{nxt}"
    # 3) (Patch 4) 独立 CTRL/Ctrl 在前 — 优先于 STIM-only
    for tok in tokens:
        if tok.lower() in ("ctrl", "control"):
            return "control"
    # 4) (Patch 6) 独立 treatment token
    for tok in tokens:
        if tok.lower() in TREATMENT_LITERAL:
            return tok.lower().replace("serumstarv", "serum_starved")
    # 5) (Patch 8) WT-Stim-* / WT-* (mouse muscle individual ID)
    for tok in tokens:
        if re.match(r"^WT-Stim-\d+[LR]?$", tok):
            return "stimulated"
        if re.match(r"^WT-\d+[LR]?$", tok):
            return "control"
    # 6) STIM 独立无 time → "_treated"
    for tok in tokens:
        if STIM_RE.match(tok):
            return f"{tok}_treated"
    return ""


time_values = Counter()
no_time_files = []
for fn, toks in file_tokens:
    tv = parse_time(toks)
    if tv:
        time_values[tv] += 1
    else:
        no_time_files.append(fn)

frac_tokens = Counter()
for tok in all_tokens:
    if FRAC_RE.match(tok):
        frac_tokens[tok] += 1

rep_explicit = Counter()
rep_trailing = Counter()
both_count = trail_only_count = rep_only_count = 0
for fn, toks in file_tokens:
    has_rep = any(REP_RE.match(t) for t in toks)
    has_trail = bool(toks and TRAIL_DIGIT_RE.match(toks[-1]))
    for t in toks:
        if REP_RE.match(t):
            rep_explicit[t.lower()] += 1
    if has_trail:
        rep_trailing[toks[-1]] += 1
    if has_rep and has_trail:
        both_count += 1
    elif has_trail:
        trail_only_count += 1
    elif has_rep:
        rep_only_count += 1

classified_lower = set()
classified_lower.update([cp.lower() for cp in CELL_PATTERNS])
classified_lower.update([op.lower() for op in ORGAN_PART_PATTERNS])
classified_lower.update(["mouse", "mus", "subcell", "wcl", "gelband"])
classified_lower.update(TREATMENT_LITERAL)

unknown = {}
for tok, cnt in token_counts.most_common():
    if cnt < 3:
        break
    tl = tok.lower()
    if (tl in classified_lower
        or is_machine_param(tok)
        or DATE_RE.match(tok)
        or STIM_RE.match(tok)
        or TIME_RE.match(tok)
        or COMBINED_TIME_RE.match(tok)
        or OSMSTR_TIME_RE.match(tok)
        or RELEASE_TOKEN_RE.match(tok)
        or FRAC_RE.match(tok)
        or REP_RE.match(tok)
        or TRAIL_DIGIT_RE.match(tok)
        or MODTYPE_RE.match(tok)
        or EXPDESIGN_RE.match(tok)
        or PERC_RE.match(tok)
        or GELBAND_COORD_RE.match(tok)
        or WT_ID_RE.match(tok)):
        continue
    unknown[tok] = cnt

# ============ 报告 ============

R = []
R.append("# Pattern Scan Report — PXD023690 (v3)\n")
R.append(f"**Total .raw files**: {len(files)}\n")
R.append(f"**Total unique tokens**: {len(token_counts)}\n\n---\n")

R.append("## 1. Organism\n")
R.append(f"- Mouse: **{len(mouse_files)}** | Human (incl. GelBand): **{len(human_files)}**\n")

R.append("## 2. Organism part (mouse only)\n")
for op, fs in organ_part_files.items():
    if fs:
        R.append(f"- {op}: **{len(fs)}**")

R.append("\n## 3. Cell type (HeLa/U2OS)\n")
for cp, fs in cell_files.items():
    R.append(f"- {cp}: **{len(fs)}**")

R.append("\n## 4. Modification type\n")
R.append(f"- Phospho: **{len(phospho_files)}**")
R.append(f"- Total proteome: **{len(total_files)}**")
R.append(f"- 🟡 GelBand (default total, **highlighted yellow,等博后确认**): **{len(gelband_files)}**")
R.append(f"- ⚠️  其他未分类: **{len(mod_unclassified)}**(理想 = 0)")

R.append("\n## 5. Experimental design\n")
R.append(f"- SUBCELL: **{len(subcell_files)}** | WCL: **{len(wcl_files)}** | GelBand: **{len(gelband_files)}**")
R.append(f"- 总和: {len(subcell_files) + len(wcl_files) + len(gelband_files)} / {len(files)}")

R.append("\n## 6. Sampling time tokens (前 30)\n")
for tv, cnt in time_values.most_common(30):
    R.append(f"- `{tv}`: {cnt}")
R.append(f"\n- 无 time 信息的文件: **{len(no_time_files)}**(理想 ≈ GelBand 的 36)")
if no_time_files:
    R.append("  示例:")
    for fn in no_time_files[:3]:
        R.append(f"    - `{fn}`")

R.append("\n## 7. Fraction tokens\n")
for tok, cnt in frac_tokens.most_common():
    R.append(f"- `{tok}`: {cnt}")
R.append("\n→ 映射: fraction 1,2 → **Cyt**;3,4 → **Mem**;5,6 → **Nuc**;WCL → **WholeCell**;GelBand → **GelBand_25perc_{coord}**")

R.append("\n## 8. Replicate 共现统计\n")
R.append(f"- Rep + trailing: **{both_count}** | trailing only: **{trail_only_count}** | Rep only: **{rep_only_count}**")

R.append("\n## 9. Instrument\n")
meta = json.loads((OUT_DIR / "project_metadata.json").read_text())
for i in meta.get("instruments", []):
    R.append(f"- {i}")

R.append("\n## 10. Unknown high-freq tokens (≥3 次)\n")
if not unknown:
    R.append("✅ 无遗漏 token。")
else:
    for tok, cnt in sorted(unknown.items(), key=lambda x: -x[1])[:50]:
        R.append(f"- `{tok}`: {cnt}")

(OUT_DIR / "pattern_report.md").write_text("\n".join(R))
print(f"✅ Phase 2 v3 done. Report: out/pattern_report.md")
print(f"   GelBand: {len(gelband_files)} files (will be highlighted yellow in Phase 4)")
print(f"   未分类 mod type: {len(mod_unclassified)} (理想 0)")
print(f"   无 time: {len(no_time_files)} (理想 ≈ {len(gelband_files)})")
