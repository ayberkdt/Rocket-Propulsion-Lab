"""v23: görseller büyütülür, alt boşluklar kapatılır, okur metni (konuşmacı notları) desteden çıkarılır.

Girdi : surumler/Sentinel_1_2_Propulsion_ITU_v22.pptx + grafik/s2c_yakit_v23.png
Çıktı : surumler/Sentinel_1_2_Propulsion_ITU_v23.pptx
        ozel/okur_metni_v23.md   (git'e girmez, hocaya giden dosyada yok)
Slayt 4, 5 ve 8'e dokunulmaz.
"""
import sys
from pathlib import Path

from pptx import Presentation
from pptx.enum.text import MSO_ANCHOR
from pptx.util import Inches, Pt

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v22.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v23.pptx"
NOTES = SUNUM / "ozel" / "okur_metni_v23.md"
CHART = HERE / "grafik" / "s2c_yakit_v23.png"

prs = Presentation(SRC)
slides = list(prs.slides)


def shapes(slide):
    return {sh.shape_id: sh for sh in slide.shapes}


def place(sh, x=None, y=None, w=None, h=None):
    if x is not None: sh.left = Inches(x)
    if y is not None: sh.top = Inches(y)
    if w is not None: sh.width = Inches(w)
    if h is not None: sh.height = Inches(h)


def inch(v):
    return v / 914400


# ---------- okur metni: önce ozel/ altına yazılır, sonra desteden silinir ----------
sys.path.insert(0, str(HERE))
from notlari_ayir import ayir  # noqa: E402

ayir(prs, NOTES, "v23")

# ---------- 1 · Mission Details: tablo boşluğu doldurur ----------
S = shapes(slides[0])
tbl_shape = S[62]
tbl = tbl_shape.table
scale = 4.6 / inch(tbl_shape.height)
for row in tbl.rows:
    row.height = int(row.height * scale)
tbl_shape.height = int(tbl_shape.height * scale)
for row in tbl.rows:
    for cell in row.cells:
        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        for p in cell.text_frame.paragraphs:
            for r in p.runs:
                if r.font.size:
                    r.font.size = Pt(16.5)

# ---------- 2 · Spacecraft and Launches: büyük fotoğraflar ----------
S = shapes(slides[1])
place(S[71], 0.59, 1.3, 5.5, 3.6)
place(S[72], 6.96, 1.3, 5.5, 3.6)
S[72].crop_bottom = 0.035                        # fotoğrafın altındaki küçük künye şeridi
place(S[73], 0.59, 4.95, 5.5, 0.3)
place(S[74], 6.96, 4.95, 5.5, 0.3)
for sid, y in ((970, 5.38), (972, 5.38), (971, 5.83), (973, 5.83)):
    place(S[sid], y=y)
place(S[974], y=6.33)

# ---------- 3 · Concept of Operations: üç görsel tam genişlik ----------
S = shapes(slides[2])
place(S[2], 0.59, 1.5, 3.85, 2.57)
place(S[14], 4.74, 1.5, 3.85, 2.56)
place(S[26], 8.89, 1.69, 3.85, 2.17)
for sid, x in ((500, 0.59), (501, 4.74), (502, 8.89)):
    place(S[sid], x=x, y=4.1, w=3.85)
row_map = {4.08: (4.38, None), 4.5: (4.8, 0.82), 5.52: (5.68, 0.8)}
for sid in range(503, 524):
    sh = S.get(sid)
    if sh is None:
        continue
    top = round(inch(sh.top), 2)
    if top in row_map:
        y, h = row_map[top]
        place(sh, y=y, h=h)
place(S[1001], y=6.53)

# ---------- 6 · 1 N Thruster Catalogue: çizimler büyür ----------
S = shapes(slides[5])
place(S[126], 0.59, 1.3, 5.78, 2.83)
place(S[127], 7.85, 1.3, 4.0, 2.83)
place(S[128], y=4.16)
place(S[129], y=4.16)
place(S[130], y=4.5)

# ---------- 7 · Sentinel-1 Manoeuvres: grafik sağ sütunu doldurur ----------
S = shapes(slides[6])
place(S[107], 7.75, 1.5, 4.98, 4.4)
place(S[1008], 7.8, 5.95, 4.9, 0.45)
for sid, y in ((1000, 1.5), (1001, 2.05), (1002, 2.05), (1003, 2.9), (1004, 2.9),
               (1005, 4.0), (1006, 4.5), (1007, 5.35)):
    place(S[sid], y=y)

# ---------- 9 · Propellant Budget: aynı düzen, daha büyük grafik ----------
s9 = slides[8]
card = next(sh for sh in s9.shapes if sh.shape_id == 1022 and not sh.text_frame.text)
caption = next(sh for sh in s9.shapes if sh.shape_id == 1022 and sh.text_frame.text)
caption._element.find(".//{*}cNvPr").set("id", str(max(sh.shape_id for sh in s9.shapes) + 1))  # v22'deki çift id
old_pic = next(sh for sh in s9.shapes if sh.shape_type == 13)
old_pic._element.getparent().remove(old_pic._element)
place(card, 7.75, 1.5, 4.98, 4.4)
pic = s9.shapes.add_picture(str(CHART), Inches(7.9), Inches(1.71), width=Inches(4.68))
pic.name = "S2C hydrazine waterfall"
place(caption, 7.8, 5.95, 4.9, 0.45)
S = shapes(s9)
for sid, y in ((1010, 1.5), (1014, 2.05), (1016, 2.05), (1015, 2.9), (1017, 2.9),
               (1018, 4.0), (1019, 4.5), (1020, 5.15), (1021, 5.7)):
    place(S[sid], y=y)

prs.save(OUT)
print("wrote", OUT)
print("wrote", NOTES)
