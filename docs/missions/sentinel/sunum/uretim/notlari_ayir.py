"""Okur metnini (konuşmacı notlarını) sunumdan ayırır: önce ozel/ altına yazar, sonra desteden siler.

Kullanım: python notlari_ayir.py <sunum.pptx> <okur_metni.md>   (sunum yerinde güncellenir)
Destede not yoksa .md dosyasına dokunmaz.
"""
import sys
from pathlib import Path

from pptx import Presentation
from pptx.opc.constants import RELATIONSHIP_TYPE as RT


def _title(slide):
    for sh in slide.shapes:   # İTÜ şablonunda başlık bandı metni (8.5", 0.54")
        if sh.has_text_frame and abs(sh.left / 914400 - 8.5) < 0.01 and abs(sh.top / 914400 - 0.54) < 0.01:
            return sh.text_frame.text
    return ""


def ayir(prs, md_path, surum):
    """prs içindeki notları md_path'e yazar ve siler. Yazılan not sayısını döndürür."""
    slides = list(prs.slides)
    texts = [(s.notes_slide.notes_text_frame.text.strip() if s.has_notes_slide else "") for s in slides]
    count = sum(1 for t in texts if t)
    if count:
        lines = [f"# Okur metni · {surum}", "",
                 "Konuşmacı notları. Sunum dosyasında yer almaz; bu klasör git'e girmez.", ""]
        for i, (s, t) in enumerate(zip(slides, texts), 1):
            lines += [f"## Slayt {i} · {_title(s)}", "", t or "(not yok)", ""]
        md_path = Path(md_path)
        md_path.parent.mkdir(exist_ok=True)
        md_path.write_text("\n".join(lines), encoding="utf-8")
    for s in slides:
        for rId, rel in list(s.part.rels.items()):
            if rel.reltype == RT.NOTES_SLIDE:
                s.part.drop_rel(rId)
    return count


if __name__ == "__main__":
    src, md = Path(sys.argv[1]), Path(sys.argv[2])
    prs = Presentation(src)
    n = ayir(prs, md, md.stem.split("_")[-1])
    prs.save(src)
    print(f"{src.name}: {n} slayt notu {md if n else '(yazılmadı)'} dosyasına ayrıldı ve desteden silindi")
