"""
Phase 1: 拉 PRIDE 项目的 file list + project metadata
Output: raw_files.json, project_metadata.json
"""
import json
import requests
import time
from pathlib import Path

PROJECT_ID = "PXD023690"
BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
OUT_DIR = Path("./out")
OUT_DIR.mkdir(exist_ok=True)


def fetch_files():
    """分页拉 file list,只保留 .raw 文件"""
    all_files = []
    page = 0
    while True:
        url = f"{BASE}/projects/{PROJECT_ID}/files?pageSize=100&page={page}"
        print(f"  Fetching page {page}...")
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        data = r.json()

        # PRIDE v3 API 返回结构: 直接是 list,或 {"_embedded": {"files": [...]}}
        if isinstance(data, list):
            batch = data
        elif "_embedded" in data:
            batch = data["_embedded"].get("files", [])
        else:
            batch = data.get("files", data.get("content", []))

        # On the first page, dump schema so we can debug if needed
        if page == 0:
            print(f"  First page: type={type(data).__name__}, "
                  f"len(batch)={len(batch)}")
            if batch:
                print(f"  Sample keys in first item: {list(batch[0].keys())[:8]}")

        if not batch:
            break

        # 过滤 .raw 文件 + 提取 download URL
        for f in batch:
            fname = f.get("fileName", "")
            if not fname.lower().endswith(".raw"):
                continue
            # publicFileLocations 是 list of {"name": "FTP Protocol", "value": "ftp://..."} 之类
            uri = ""
            for loc in f.get("publicFileLocations", []):
                v = loc.get("value", "")
                if v.startswith("http") or v.startswith("ftp"):
                    uri = v
                    break
            all_files.append({
                "fileName": fname,
                "uri": uri,
                "fileCategory": f.get("fileCategory", {}),
            })

        page += 1
        time.sleep(0.3)
        if len(batch) < 100:
            break  # 最后一页

    return all_files


def fetch_project_metadata():
    """拉 project-level metadata,主要为了 instrument model"""
    url = f"{BASE}/projects/{PROJECT_ID}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.json()


if __name__ == "__main__":
    print(f"[1/2] Fetching file list for {PROJECT_ID}...")
    files = fetch_files()
    print(f"  -> {len(files)} .raw files")
    (OUT_DIR / "raw_files.json").write_text(json.dumps(files, indent=2))

    print(f"[2/2] Fetching project metadata...")
    meta = fetch_project_metadata()
    (OUT_DIR / "project_metadata.json").write_text(json.dumps(meta, indent=2))
    instruments = meta.get("instruments", [])
    print(f"  -> instruments: {[i.get('name', i) if isinstance(i, dict) else i for i in instruments]}")

    print(f"\n[OK] Phase 1 done. Files: {len(files)}, see out/raw_files.json")
