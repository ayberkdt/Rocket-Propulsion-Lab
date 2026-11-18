"""v24: kaynak makalelerden daha iyi görseller.

Slayt 3: Fig.16 (okunmayan tandem grafiği) ve ESA 1B infografiği yerine
         S1A LEOP açılış sırası (ISSFD 2015 #137 Fig.4) ve S1B bertaraf izi (SpaceOps 2025 #226 Fig.9).
Slayt 5: tablonun yanına S1B itici yerleşimi (SpaceOps 2025 #360 Fig.1).
Slayt 6: ArianeGroup çizimi yerine uyduya monte 1 N itici fotoğrafı (ArianeGroup).
Diğer slaytlara dokunulmaz. Girdi: surumler/..._v23.pptx  ·  Çıktı: surumler/..._v24.pptx
"""
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v23.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v24.pptx"
FIG = HERE / "gorsel_adaylari"

prs = Presentation(SRC)
slides = list(prs.slides)


def shapes(slide):
    return {sh.shape_id: sh for sh in slide.shapes}


def place(sh, x=None, y=None, w=None, h=None):
    if x is not None: sh.left = Inches(x)
    if y is not None: sh.top = Inches(y)
    if w is not None: sh.width = Inches(w)
    if h is not None: sh.height = Inches(h)


def set_text(shape, text):
    runs = shape.text_frame.paragraphs[0].runs
    runs[0].text = text
    for r in runs[1:]:
        r._r.getparent().remove(r._r)


def replace_picture(slide, old, path, x, y, w=None, h=None, name=None):
    """Yeni resmi eskisinin z-sırasına koyar, eskisini siler."""
    pic = slide.shapes.add_picture(str(path), Inches(x), Inches(y),
                                   Inches(w) if w else None, Inches(h) if h else None)
    old._element.addprevious(pic._element)
    old._element.getparent().remove(old._element)
    if name:
        pic.name = name
    return pic


# ---------- 3 · Concept of Operations ----------
s = slides[2]
S = shapes(s)
top, band = 1.5, 2.65                                   # görsel şeridi
place(S[2], 0.59, top + (band - 2.0) / 2, 3.0, 2.0)    # 1C fırlatma fotoğrafı
replace_picture(s, S[14], FIG / "D_s1a_leop_acilis_issfd137_fig4.png", 3.765, top, w=3.6,
                name="S1A deployment sequence")
replace_picture(s, S[26], FIG / "C_s1b_bertaraf_spaceops226_fig9.png", 7.54, top + (band - 2.39) / 2, w=5.2,
                name="S1B orbit lowering plan vs actual")
place(S[500], 0.59, 4.18, 3.0)
place(S[501], 3.765, 4.18, 3.6)
place(S[502], 7.54, 4.18, 5.2)
set_text(S[501], "Sentinel-1A deployment sequence · ISSFD 2015")
set_text(S[502], "Sentinel-1B orbit lowering, planned vs actual · SpaceOps 2025")
row_map = {4.38: (4.45, None), 4.8: (4.87, 0.78), 5.68: (5.71, 0.78)}
for sid in range(503, 524):
    sh = S.get(sid)
    if sh is None:
        continue
    t = round(sh.top / 914400, 2)
    if t in row_map:
        y, hh = row_map[t]
        place(sh, y=y, h=hh)
place(S[1001], y=6.53)
set_text(S[525], "Sources   Arianespace VV24/VV25, ESA Flight Operations, ISSFD 2015 paper 137, "
                 "SpaceOps 2025 papers 226 and 617, Airbus GNC 2014")

# ---------- 5 · Hydrazine Propulsion ----------
s = slides[4]
S = shapes(s)
tbl_shape = S[119]
old_w = tbl_shape.width
new_w = Inches(9.3)
for col in tbl_shape.table.columns:
    col.width = int(col.width * new_w / old_w)
tbl_shape.width = sum(c.width for c in tbl_shape.table.columns)
layout_h = 4.12
pic = s.shapes.add_picture(str(FIG / "A_s1b_itici_yerlesimi_spaceops360_fig1.png"),
                           Inches(0.59 + 9.3 + 0.25), Inches(2.25), height=Inches(layout_h))
pic.name = "S1B thruster layout"
cap = S[120]
set_text(cap, "2C tank charging pressure: 21 bar. ESA gives no absolute or gauge reference.  "
              "Thruster layout: Sentinel-1B, SpaceOps 2025 paper 360, Fig. 1.")

# ---------- 6 · 1 N Thruster Catalogue ----------
s = slides[5]
S = shapes(s)
photo = replace_picture(s, S[127], FIG / "G_arianegroup_1n_itici_monte.jpg", 8.6, 1.3, h=2.83,
                        name="ArianeGroup 1 N thruster photo")
photo.crop_top = 0.2                                   # üstteki kablo demetini kırp
photo.width = int(Inches(2.83) * (900 / (1254 * 0.8)))
photo.height = Inches(2.83)
photo.left = int(Inches(6.96) + (Inches(5.78) - photo.width) / 2)
set_text(S[129], "ArianeGroup 1 N thruster, integrated")
set_text(S[131], "Supplier product specifications. Rafael p. 46 and ArianeGroup brochure p. 6. "
                 "Photo: ArianeGroup.")

prs.save(OUT)
print("wrote", OUT)
