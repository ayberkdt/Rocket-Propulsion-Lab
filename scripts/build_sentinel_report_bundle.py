"""Build printable offline HTML and a source-bound archive of the Sentinel assignment.

Authoring dependencies: mistune, beautifulsoup4. Calculators/tests use only Python's stdlib.
Supply the source directory containing source_receipts.json; no network access is performed.
"""
import argparse
import base64
import hashlib
import html
import json
import mimetypes
import re
import shutil
import subprocess
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import mistune
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
CSS = """
*{box-sizing:border-box}body{margin:0;background:#eef3f7;color:#162d3c;
font:17px/1.65 'Segoe UI',Arial,sans-serif}main{max-width:1160px;margin:24px auto;
padding:38px 52px;background:white;border-top:6px solid #176b91}
h1{line-height:1.2;font-size:2rem}h2{margin-top:2.2em;border-bottom:2px solid #d7e3ea}
h3{margin-top:1.8em}a{color:#075d83;overflow-wrap:anywhere}img{max-width:100%;height:auto;
display:block;margin:24px auto}nav{background:#edf5f9;padding:18px 24px}
nav ul{columns:2;padding-left:20px}nav li{break-inside:avoid}.table-scroll{overflow-x:auto}
table{width:100%;border-collapse:collapse;font-size:.87rem;line-height:1.5;margin:18px 0}
td,th{border:1px solid #cbdce5;padding:9px;vertical-align:top;overflow-wrap:anywhere}
th{background:#e5f0f6;text-align:left}tr:nth-child(even){background:#f6f9fb}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f0f5f8;padding:18px}
code{font: .87em Consolas,monospace}footer{border-top:1px solid #cbdce5;padding-top:20px;
margin-top:36px;color:#4a6473;font-size:.85rem}
@media(max-width:750px){main{margin:0;padding:24px 16px}nav ul{columns:1}h1{font-size:1.6rem}}
@media print{body{background:white;font-size:10pt}main{margin:0;padding:0;border:0;
max-width:none}nav{display:none}table{font-size:8pt}.table-scroll{overflow:visible}
h2,h3{break-after:avoid}tr,img,pre{break-inside:avoid}a{color:inherit}}
"""


def render(path: Path, output_root: Path) -> Path:
    markdown = mistune.create_markdown(plugins=["table"])
    soup = BeautifulSoup(markdown(path.read_text(encoding="utf-8")), "html.parser")
    toc = []
    for i, heading in enumerate(soup.find_all(["h2", "h3"]), 1):
        heading["id"] = f"section-{i}"
        if heading.name == "h2":
            toc.append(f'<li><a href="#section-{i}">{html.escape(heading.get_text())}</a></li>')
    for table in list(soup.find_all("table")):
        table.wrap(soup.new_tag("div", attrs={"class": "table-scroll", "tabindex": "0"}))
    for picture in soup.find_all("img"):
        target = (path.parent / picture["src"]).resolve()
        if not target.is_relative_to(output_root.resolve()):
            raise ValueError(f"image escaped bundle: {target}")
        mime = mimetypes.guess_type(target.name)[0]
        picture["src"] = f"data:{mime};base64," + base64.b64encode(target.read_bytes()).decode()
    for link in soup.find_all("a", href=True):
        href = link["href"]
        if href.endswith(".md") and (path.parent / href).is_file():
            link["href"] = href[:-3] + ".html"
    title = soup.find("h1").get_text()
    home = "../" * len(path.relative_to(output_root).parent.parts) + "BASLANGIC.html"
    body = f'<p><a href="{home}">← Karşılaştırmalı ödev paketi</a></p>' + str(soup)
    doc = (f'<!doctype html><html lang="tr"><head><meta charset="utf-8">'
           f'<meta name="viewport" content="width=device-width,initial-scale=1">'
           f'<title>{html.escape(title)}</title><style>{CSS}</style></head><body><main>'
           f'<nav aria-label="İçindekiler"><ul>{"".join(toc)}</ul></nav>{body}'
           '<footer>V4 · 29 Ekim 2026 · Kaynak / miras / hesap / varsayım ayrıdır. '
           'Eksik uçuş bütçeleri açıkça işaretlidir. Görseller dosyaya gömülüdür.</footer>'
           '</main></body></html>')
    target = path.with_suffix(".html")
    target.write_text(doc, encoding="utf-8")
    return target


def check_html(path: Path, output_root: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if "\ufffd" in text:
        raise ValueError(f"encoding replacement character: {path}")
    soup = BeautifulSoup(text, "html.parser")
    ids = {node["id"] for node in soup.find_all(id=True)}
    for anchor in soup.find_all("a", href=True):
        href = anchor["href"]
        if href.startswith("#"):
            if href[1:] not in ids:
                raise ValueError(f"broken fragment: {path} {href}")
        elif not re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", href):
            target = (path.parent / href).resolve()
            if not target.is_relative_to(output_root.resolve()) or not target.is_file():
                raise ValueError(f"broken local link: {path} {href}")
    if any(not pic["src"].startswith("data:image/") for pic in soup.find_all("img")):
        raise ValueError("report pictures must be embedded for offline HTML")
    return {"path": path.relative_to(output_root).as_posix(),
            "tables": len(soup.find_all("table")), "images": len(soup.find_all("img"))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    if out == ROOT or ROOT.is_relative_to(out):
        raise ValueError("output must not contain the source repository")
    out.mkdir(parents=True, exist_ok=True)
    receipt_path = out / "source_receipts.json"
    records = json.loads(receipt_path.read_text(encoding="utf-8"))
    for rec in records:
        if rec["status"] == "downloaded":
            p = out / rec["local_path"]
            if not p.resolve().is_relative_to(out) or hashlib.sha256(p.read_bytes()).hexdigest() != rec["sha256"]:
                raise ValueError(f"source hash mismatch: {p}")
    relative_docs = Path("docs/missions/sentinel")
    shutil.copytree(ROOT / relative_docs, out / relative_docs, dirs_exist_ok=True)
    shutil.copy2(ROOT / "docs/SASS_PROJECT_PREPARATION.md", out / "docs")
    for directory, names in {"scripts": ["sentinel_series_budget.py", "plot_sentinel_series.py",
                                         "build_sentinel_report_bundle.py"],
                             "tests": ["test_sentinel_series_budget.py"]}.items():
        (out / directory).mkdir(exist_ok=True)
        for name in names:
            shutil.copy2(ROOT / directory / name, out / directory / name)
    base = "docs/missions/sentinel/"
    entries = [("Sentinel-1: A/B/C/D", "SENTINEL-1_SERISI"),
               ("Sentinel-2: A/B/C/D", "SENTINEL-2_SERISI"),
               ("Denklemler ve hesaplar", "HESAPLAR"),
               ("27 maddelik kontrol listesi", "KONTROL_LISTESI"),
               ("Kaynaklar ve erişim", "KAYNAKLAR"),
               ("Ayrıntılı Sentinel-1C", "SENTINEL-1C"),
               ("Ayrıntılı Sentinel-2C", "SENTINEL-2C")]
    links = "".join(f'<li><a href="{base}{name}.html">{label}</a></li>' for label, name in entries)
    index = (f'<!doctype html><html lang="tr"><meta charset="utf-8">'
             f'<meta name="viewport" content="width=device-width,initial-scale=1">'
             f'<title>SASS — Sentinel serileri</title><style>{CSS}</style><main>'
             '<h1>SASS: Sentinel-1 ve Sentinel-2 serileri</h1>'
             '<p>29 Ekim 2026 · A/B/C/D karşılaştırması · İtki ve besleme odağı</p>'
             f'<nav><ul>{links}</ul></nav><p>Gerçek görev fotoğrafları, işlevsel şemalar, '
             'kaynaklı donanım verileri, CONOPS, bütçeler ve yeniden üretilebilir hesaplar.</p>'
             '<p>Yakıt hesapları ölçüm değildir; eksik uçuş bütçeleri ve varsayımlar raporda belirtilir.</p>'
             '<p><a href="source_receipts.json">İndirme / URL / SHA-256 kayıtları</a> · '
             f'<a href="{base}SERI_SONUCLARI.json">Hesap sonuçları</a> · '
             f'<a href="{base}README.html">Hesapları çalıştırma</a></p>'
             '<p>HTML çevrimdışı okunabilir; görseller gömülüdür. Tam kaynaklar sources/ dizinindedir. '
             'Dosyaların bulunması yeniden dağıtım lisansı vermez.</p></main></html>')
    (out / "BASLANGIC.html").write_text(index, encoding="utf-8")
    pages = [render(path, out) for path in (out / "docs").rglob("*.md")]
    checks = [check_html(path, out) for path in [out / "BASLANGIC.html", *pages]]
    completed = subprocess.run([sys.executable, "-I", "-S",
                                str(out / "tests/test_sentinel_series_budget.py"), "-q"],
                               capture_output=True, encoding="utf-8", timeout=60, check=True)
    verification = {"source_hashes_checked": sum(r["status"] == "downloaded" for r in records),
                    "html": checks, "isolated_test_output": completed.stderr.strip(),
                    "scope": "Numerical, source-hash and structural checks; browser layout and flight hardware not validated"}
    (out / "VALIDATION.json").write_text(json.dumps(verification, ensure_ascii=False, indent=2)
                                          + "\n", encoding="utf-8")
    excluded = {"SASS_Sentinel_Serileri_V4.zip", "SASS_Sentinel_Serileri_V4.sha256", "MANIFEST.json"}
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name not in excluded
                   and "__pycache__" not in p.parts)
    manifest = {p.relative_to(out).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    manifest_path = out / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    archive = out / "SASS_Sentinel_Serileri_V4.zip"
    with ZipFile(archive, "w", ZIP_DEFLATED, compresslevel=6) as z:
        for path in [*files, manifest_path]:
            z.write(path, path.relative_to(out).as_posix())
    with ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ValueError("ZIP CRC verification failed")
        for name, digest in manifest.items():
            if hashlib.sha256(z.read(name)).hexdigest() != digest:
                raise ValueError(f"ZIP hash mismatch: {name}")
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (out / "SASS_Sentinel_Serileri_V4.sha256").write_text(f"{digest}  {archive.name}\n", encoding="ascii")
    print(json.dumps({"archive": str(archive), "bytes": archive.stat().st_size,
                      "sha256": digest, "files_verified": len(manifest),
                      "sources_verified": verification["source_hashes_checked"]}, indent=2))


if __name__ == "__main__":
    main()
