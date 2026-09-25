# -*- coding: utf-8 -*-
"""
博士课题工作目录盘点 + 脚本追溯。
用法（Windows，Python 3.8+，只依赖标准库）：
    python inventory.py
    python inventory.py --out "D:\\博士\\inventory_2026-09-25.md" --root "D:\\博士\\Phospho" ...
不传 --root 时扫描默认的四个路径。只写一个输出文件，表之外不写任何文字。
"""
import argparse
import hashlib
import json
import os
import re
import sys
from datetime import datetime

DEFAULT_ROOTS = [
    r"D:\博士\博后_大导开会会议记录",
    r"D:\博士\protein_contour_wiki",
    r"D:\博士\Protein contour",
    r"D:\博士\Phospho",
]
DEFAULT_OUT = r"D:\博士\inventory_2026-09-25.md"
BIG = 500 * 1024 * 1024

EXT_RAW = {".raw", ".mzml", ".mzxml", ".wiff", ".scan", ".baf", ".tdf", ".tdf_bin",
           ".mgf", ".dia", ".yep", ".fid", ".mzdata", ".pkl.raw", ".ms1", ".ms2", ".lcd", ".qgd"}
RAW_DIR_EXT = {".d"}  # Bruker .d 是目录，整体按一个原始数据文件记
EXT_SCRIPT = {".py", ".r", ".ipynb", ".sh", ".rmd", ".bat", ".ps1", ".pl", ".m", ".jl", ".sql", ".cmd"}
SCRIPT_READ = {".py", ".r", ".ipynb", ".sh", ".rmd"}
EXT_TABLE = {".csv", ".tsv", ".txt", ".xlsx", ".xls", ".xlsm", ".tab", ".rds", ".rdata", ".rda",
             ".parquet", ".feather", ".h5", ".hdf5", ".pkl", ".pickle", ".npy", ".npz", ".json",
             ".fasta", ".fa", ".faa", ".gz", ".zip", ".7z", ".rar", ".sqlite", ".db", ".xml",
             ".gmt", ".gff", ".gtf", ".bed", ".sav", ".mat"}
EXT_FIG = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".svg", ".ai", ".eps", ".gif", ".bmp", ".emf", ".wmf", ".pdf"}
EXT_DOC = {".docx", ".doc", ".md", ".pptx", ".ppt", ".tex", ".bib", ".rtf", ".html", ".htm", ".odt", ".pages", ".key", ".one", ".enex", ".eml", ".msg"}
MEETING_HINT = re.compile(r"(会议|开会|纪要|meeting|minutes|记录)", re.I)

# 路径字面量里可能出现的扩展名（用来在脚本里找输入/输出）
DATA_EXT_RE = r"(?:csv|tsv|txt|xlsx?|xlsm|tab|rds|rdata|rda|parquet|feather|h5|hdf5|pkl|pickle|npy|npz|json|fasta|fa|faa|gz|zip|sqlite|db|xml|gmt|gff|gtf|bed|mat|png|jpe?g|tiff?|svg|pdf|eps|ai|raw|mzml|mzxml|wiff|mgf|d|py|R|r|ipynb|sh|Rmd|rmd|docx?|pptx?|md|html?|sav)"
LITERAL_RE = re.compile(r"""(?P<q>["'])(?P<s>[^"'\n]{1,400}?\.""" + DATA_EXT_RE + r""")(?P=q)""")
FSTRING_RE = re.compile(r"""[fF]["'][^"'\n]*\{[^"'\n]*["']""")
SEP_RE = re.compile(r"[\\/]")

WRITE_HINT = re.compile(
    r"(to_csv|to_excel|to_parquet|to_pickle|to_hdf|to_json|to_feather|savefig|save\(|savetxt|np\.save|dump\(|"
    r"write\.csv|write\.table|write_csv|write_tsv|write_xlsx|writexl|saveRDS|save\.image|ggsave|"
    r"pdf\(|png\(|tiff\(|svg\(|jpeg\(|dev\.copy|fwrite|write_rds|write\.xlsx|openxlsx|"
    r"open\([^)]*['\"][wax]\+?b?['\"]|\bwrite\b|\bwriteLines\b|\bsink\(|"
    r">\s*['\"]|>>\s*['\"]|-o\s|--out|output|\bout\b|\bsave\b|export)",
    re.I)
READ_HINT = re.compile(
    r"(read_csv|read_excel|read_table|read_parquet|read_pickle|read_hdf|read_json|read_feather|"
    r"np\.load|loadtxt|genfromtxt|pd\.read|open\([^)]*['\"]r|open\(\s*['\"]|json\.load|"
    r"read\.csv|read\.table|read\.delim|readRDS|readxl|read_xlsx|read_excel|fread|load\(|source\(|"
    r"readLines|read_tsv|read_lines|import_|\bread\b|\bload\b|input|\bin\b|cat\s|<\s*['\"])",
    re.I)


def norm(p):
    return os.path.normcase(os.path.normpath(p))


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return ("%d %s" % (n, unit)) if unit == "B" else ("%.1f %s" % (n, unit))
        n /= 1024.0


def ext_of(name):
    low = name.lower()
    if low.endswith(".wiff.scan"):
        return ".scan"
    if low.endswith(".tar.gz"):
        return ".gz"
    return os.path.splitext(low)[1]


def classify(path, name, is_dir_unit=False):
    e = ext_of(name)
    if is_dir_unit or e in EXT_RAW:
        return "原始质谱数据"
    if e in EXT_SCRIPT:
        return "脚本"
    if MEETING_HINT.search(path) and (e in EXT_DOC or e in {".txt", ".pdf", ".m4a", ".mp3", ".wav", ".mp4"}):
        return "会议记录"
    if e in EXT_TABLE:
        return "数据表"
    if e in EXT_FIG:
        return "图"
    if e in EXT_DOC:
        return "文档"
    return "其他"


# ---------------------------------------------------------------- 第一步：盘点
def walk(roots):
    """返回 {norm_path: info}。.d 目录整体记为一个原始数据单元，不进入。"""
    files = {}
    for root in roots:
        if not os.path.exists(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            keep = []
            for d in dirnames:
                full = os.path.join(dirpath, d)
                if ext_of(d) in RAW_DIR_EXT:
                    size = 0
                    mtime = 0
                    for dp, _, fns in os.walk(full):
                        for fn in fns:
                            try:
                                st = os.stat(os.path.join(dp, fn))
                                size += st.st_size
                                mtime = max(mtime, st.st_mtime)
                            except OSError:
                                pass
                    files[norm(full)] = dict(path=full, name=d, size=size, mtime=mtime,
                                             type="原始质谱数据", root=root, ext=".d", noopen=True)
                elif d in {".git", "__pycache__", ".ipynb_checkpoints", ".Rproj.user", "node_modules", ".snakemake"}:
                    continue
                else:
                    keep.append(d)
            dirnames[:] = keep
            for fn in filenames:
                full = os.path.join(dirpath, fn)
                try:
                    st = os.stat(full)
                except OSError:
                    continue
                t = classify(full, fn)
                files[norm(full)] = dict(path=full, name=fn, size=st.st_size, mtime=st.st_mtime,
                                         type=t, root=root, ext=ext_of(fn),
                                         noopen=(t == "原始质谱数据" or st.st_size > BIG))
    return files


# ---------------------------------------------------------------- 第二步：读脚本
def script_text(info):
    p = info["path"]
    try:
        if info["ext"] == ".ipynb":
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                nb = json.load(f)
            cells = nb.get("cells", [])
            return "\n".join("".join(c.get("source", [])) for c in cells if c.get("cell_type") == "code")
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            return f.read()
    except Exception:
        return ""


def has_variable(line_before, literal):
    """字面量是否被变量拼接：f-string、+、%、format、paste/paste0、file.path、os.path.join(非全字面量)、sprintf。"""
    seg = line_before
    if re.search(r"""[fF]["']$""", seg):
        return True
    if "{" in literal and "}" in literal:
        return True
    if re.search(r"(paste0?|file\.path|sprintf|glue|os\.path\.join|Path\(|joinpath|format\(|\.format|%\s*\(|\+\s*$)", seg):
        return True
    return False


def extract_refs(text, script_path):
    """返回 [(literal, mode, variable_built)]，mode ∈ {'r','w','?'}。"""
    refs = []
    for m in LITERAL_RE.finditer(text):
        lit = m.group("s")
        start = m.start()
        ls = text.rfind("\n", 0, start) + 1
        le = text.find("\n", m.end())
        le = len(text) if le < 0 else le
        line = text[ls:le]
        before = text[ls:start]
        # 忽略 import 行、URL、注释
        if line.lstrip().startswith(("#", "//")) or "://" in lit:
            continue
        var = has_variable(before, lit)
        # 判断读写：先看同一行，再看前一行
        prev_ls = text.rfind("\n", 0, max(ls - 1, 0)) + 1
        ctx = text[prev_ls:le]
        w = bool(WRITE_HINT.search(ctx))
        r = bool(READ_HINT.search(ctx))
        if w and not r:
            mode = "w"
        elif r and not w:
            mode = "r"
        elif w and r:
            # 同时命中：以字面量前的最近函数名为准
            pw = [x.start() for x in WRITE_HINT.finditer(before)]
            pr = [x.start() for x in READ_HINT.finditer(before)]
            mode = "w" if (pw and (not pr or pw[-1] > pr[-1])) else ("r" if pr else "?")
        else:
            mode = "?"
        refs.append((lit, mode, var, line.strip()))
    # 脚本调用其他脚本
    calls = []
    for m in re.finditer(r"""(?:source|exec|execfile|%run|Rscript|python3?|bash|sh|subprocess\.\w+\([^)]*|system\()\s*\(?\s*["']?([^"'\s)]+\.(?:py|R|r|sh|ipynb|Rmd))""", text):
        calls.append(m.group(1))
    for m in re.finditer(r"""^\s*(?:from\s+\.?([\w\.]+)\s+import|import\s+([^\n#]+))""", text, re.M):
        if m.group(1):
            calls.append(m.group(1).split(".")[0] + ".py")
        else:
            for part in m.group(2).split(","):
                mod = part.strip().split(" as ")[0].strip().split(".")[0]
                if mod:
                    calls.append(mod + ".py")
    for m in re.finditer(r"""(?:library|require)\(["']?([\w\.]+)["']?\)""", text):
        calls.append(m.group(1) + ".R")
    return refs, calls


# ---------------------------------------------------------------- 路径解析
def resolve(literal, script_path, files, by_name, roots):
    """把脚本里的字面量对应到盘点表里的文件。
    返回 (matches, kind)：kind ∈ {'abs','rel','name'}；matches 为 norm_path 列表。"""
    lit = literal.strip().replace("\\\\", "\\")
    is_abs = bool(re.match(r"^[A-Za-z]:[\\/]", lit)) or lit.startswith(("\\\\", "/"))
    if is_abs:
        n = norm(lit)
        return ([n] if n in files else []), "abs"
    cands = set()
    sdir = os.path.dirname(script_path)
    bases = [sdir, os.path.dirname(sdir), os.path.dirname(os.path.dirname(sdir))] + list(roots)
    for b in bases:
        n = norm(os.path.join(b, lit))
        if n in files:
            cands.add(n)
    if cands:
        return sorted(cands), "rel"
    # 只按文件名匹配
    base = os.path.basename(SEP_RE.sub(os.sep, lit))
    return sorted(by_name.get(base.lower(), [])), "name"


def sha1(path, limit=BIG):
    h = hashlib.sha1()
    try:
        with open(path, "rb") as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                h.update(b)
    except OSError:
        return None
    return h.hexdigest()


# ---------------------------------------------------------------- 第二步：追溯
def trace(files, roots):
    by_name = {}
    for n, info in files.items():
        by_name.setdefault(info["name"].lower(), []).append(n)

    scripts = {n: i for n, i in files.items() if i["type"] == "脚本" and i["ext"] in SCRIPT_READ and not i["noopen"]}
    evidence = {n: [] for n in files}          # 确定引用：(脚本名, r/w)
    uncertain = {n: [] for n in files}         # 不确定引用：(脚本名, 原因)
    script_called = {n: set() for n in scripts}
    script_uncertain_refs = {n: [] for n in scripts}
    unresolved = {n: [] for n in scripts}      # 脚本里指向不存在文件的路径

    for sn, sinfo in scripts.items():
        text = script_text(sinfo)
        refs, calls = extract_refs(text, sinfo["path"])
        sname = sinfo["name"]
        for lit, mode, var, line in refs:
            matches, kind = resolve(lit, sinfo["path"], files, by_name, roots)
            tag = {"r": "读", "w": "写", "?": "引用"}[mode]
            if not matches:
                unresolved[sn].append(lit)
                continue
            if var:
                for m in matches:
                    uncertain[m].append((sname, "变量拼接路径：%s" % line[:120]))
            elif kind == "abs" or (kind == "rel" and len(matches) == 1):
                for m in matches:
                    evidence[m].append((sname, tag))
            elif kind == "rel":
                for m in matches:
                    uncertain[m].append((sname, "相对路径 %s 可对应多个文件" % lit))
            else:  # name
                if len(matches) == 1:
                    uncertain[matches[0]].append((sname, "相对路径 %s 未能按脚本所在目录解析，仅文件名一致" % lit))
                else:
                    for m in matches:
                        uncertain[m].append((sname, "相对路径 %s 仅文件名一致且有 %d 个同名文件" % (lit, len(matches))))
        for c in calls:
            matches, kind = resolve(c, sinfo["path"], files, by_name, roots)
            for m in matches:
                if m in scripts and m != sn:
                    script_called[m].add(sname)

    # 同名 / 同内容 的更新版本
    newer_by_name = {}
    for name, lst in by_name.items():
        if len(lst) < 2:
            continue
        newest = max(lst, key=lambda n: files[n]["mtime"])
        for n in lst:
            if n != newest and files[n]["mtime"] < files[newest]["mtime"]:
                newer_by_name[n] = newest
    by_hash = {}
    for n, info in files.items():
        if info["noopen"] or info["type"] == "原始质谱数据" or info["size"] == 0:
            continue
        key = (info["size"],)
        by_hash.setdefault(key, []).append(n)
    newer_by_content = {}
    for key, lst in by_hash.items():
        if len(lst) < 2:
            continue
        groups = {}
        for n in lst:
            h = sha1(files[n]["path"])
            if h:
                groups.setdefault(h, []).append(n)
        for h, g in groups.items():
            if len(g) < 2:
                continue
            newest = max(g, key=lambda n: files[n]["mtime"])
            for n in g:
                if n != newest and files[n]["mtime"] < files[newest]["mtime"]:
                    newer_by_content[n] = newest
    return scripts, evidence, uncertain, script_called, unresolved, newer_by_name, newer_by_content


ENTRY_HINT = re.compile(r"(main|run|pipeline|analysis|analyse|analyze|全流程|主|workflow|__main__|snakefile|makefile)", re.I)


def decide(files, scripts, evidence, uncertain, script_called, newer_by_name, newer_by_content):
    rows = []
    unsure_rows = []
    for n in sorted(files, key=lambda x: files[x]["path"].lower()):
        info = files[n]
        ev = evidence[n]
        un = uncertain[n]
        newer = newer_by_content.get(n) or newer_by_name.get(n)
        why = ""
        if info["type"] == "脚本":
            callers = script_called.get(n, set())
            if callers:
                status, basis = "活", "被调用：" + "、".join(sorted(callers))
            elif info["ext"] not in SCRIPT_READ or info["noopen"]:
                status, basis = "不确定", "未读取内容"
                why = "扩展名不在追溯范围（.py .R .ipynb .sh）或文件过大，未解析"
            elif ENTRY_HINT.search(info["name"]) or (info["ext"] == ".ipynb"):
                status, basis = "活", "入口脚本（未被其他脚本调用，按文件名/notebook 判定）"
            else:
                status, basis = "不确定", "未被其他脚本调用"
                why = "脚本未被任何脚本调用，也不能从文件名判断是否为入口脚本"
            if newer:
                basis += "；另有更新版本：" + files[newer]["path"]
        else:
            if ev:
                names = sorted({s for s, _ in ev})
                tags = sorted({t for _, t in ev})
                status = "活"
                basis = "、".join(names) + "（" + "/".join(tags) + "）"
                if newer:
                    basis += "；另有更新版本：" + files[newer]["path"]
            elif newer:
                status, basis = "死", "有更新版本：" + files[newer]["path"]
            elif un:
                status = "不确定"
                basis = "、".join(sorted({s for s, _ in un})) + "（引用无法确认）"
                why = "；".join(sorted({"%s：%s" % (s, r) for s, r in un}))
            elif info["type"] in ("文档", "会议记录"):
                status, basis = "不确定", "无引用"
                why = info["type"] + "本来就不被脚本引用"
            elif info["type"] == "其他":
                status, basis = "不确定", "无引用"
                why = "扩展名无法归类，不能判断是否应被脚本引用"
            else:
                status, basis = "死", "无引用"
        rows.append((info, status, basis))
        if status == "不确定":
            unsure_rows.append((info, why))
    return rows, unsure_rows


def md_cell(s):
    return str(s).replace("|", "\\|").replace("\n", " ")


def write_md(out, rows, unsure_rows):
    lines = ["| 路径 | 大小 | 修改日期 | 类型 | 状态 | 依据 |", "|---|---|---|---|---|---|"]
    for info, status, basis in rows:
        d = datetime.fromtimestamp(info["mtime"]).strftime("%Y-%m-%d") if info["mtime"] else ""
        lines.append("| %s | %s | %s | %s | %s | %s |" % (
            md_cell(info["path"]), human(info["size"]), d, info["type"], status, md_cell(basis)))
    lines += ["", "| 路径 | 为什么不确定 |", "|---|---|"]
    for info, why in unsure_rows:
        lines.append("| %s | %s |" % (md_cell(info["path"]), md_cell(why)))
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", help="要扫描的母文件夹，可多次传；默认四个路径")
    ap.add_argument("--out", default=DEFAULT_OUT)
    a = ap.parse_args()
    roots = a.root or DEFAULT_ROOTS
    files = walk(roots)
    scripts, evidence, uncertain, script_called, unresolved, nb, nc = trace(files, roots)
    rows, unsure = decide(files, scripts, evidence, uncertain, script_called, nb, nc)
    write_md(a.out, rows, unsure)


if __name__ == "__main__":
    main()
