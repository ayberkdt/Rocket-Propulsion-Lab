"""v25: hoca değerlendirmesindeki zorunlu düzeltmeler.

Slayt 1: 2A uzatması end-2027 değil end-2026 (ESA, Mayıs 2026); filo durumuna kaynak.
Slayt 2: yörüngeye bırakılış irtifalarına kaynak.
Slayt 3: Troll ve 1C/1A fazlama satırlarına kaynak.
Slayt 7: doğrulanamayan tork grafiği yerine S1B itici yönleri (SpaceOps 2025 #226 Fig.3);
         metin kaynaktaki "yaklaşık 10 kat" ifadesine indirildi; virgüllü ondalık kalktı.
Slayt 8: 7.97 m/s satırı "planlanan, uçulmamış" olarak yazıldı.
Slayt 9: grafik ve büyük sayılar tek ondalık; 47 kg planın ideal 38 kg'a göre payı yazıldı.
Girdi: surumler/..._v24.pptx  ·  Çıktı: surumler/..._v25.pptx
"""
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v24.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v25.pptx"
FIG = HERE / "gorsel_adaylari"
CHART = HERE / "grafik" / "s2c_yakit_v25.png"

prs = Presentation(SRC)
slides = list(prs.slides)


def shapes(slide):
    return {sh.shape_id: sh for sh in slide.shapes}


def set_text(shape, text):
    p = shape.text_frame.paragraphs[0]
    p.runs[0].text = text
    for r in p.runs[1:]:
        r._r.getparent().remove(r._r)


def find_text(slide, start):
    hits = [sh for sh in slide.shapes if sh.has_text_frame and sh.text_frame.text.startswith(start)]
    assert len(hits) == 1, (start, len(hits))
    return hits[0]


def replace_cell(table, old, new):
    n = 0
    for row in table.rows:
        for cell in row.cells:
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    if old in r.text:
                        r.text = r.text.replace(old, new)
                        n += 1
    assert n == 1, (old, n)


def replace_picture(slide, old, path, x, y, w=None, h=None, name=None):
    pic = slide.shapes.add_picture(str(path), Inches(x), Inches(y),
                                   Inches(w) if w else None, Inches(h) if h else None)
    old._element.addprevious(pic._element)
    old._element.getparent().remove(old._element)
    if name:
        pic.name = name
    return pic


# ---------- 1 · Mission Details ----------
s = slides[0]
replace_cell(shapes(s)[62].table, "2A extended to end-2027", "2A extended to end-2026")
set_text(find_text(s, "Fleet status on"),
         "Fleet status on 12 Nov 2026. ESA: 1A ended 29 Jun 2026; 2A extended to end-2026; "
         "2D in storage (CEOS WGISS-61).")

# ---------- 2 · Spacecraft and Launches ----------
s = slides[1]
set_text(find_text(s, "Kourou, French Guiana"),
         "Kourou. Photos: Thales Alenia Space/ImagIN, ESA–CNES–Arianespace. "
         "Insertion: VV25 p.2 (1C); SpaceOps 2025 paper 617 p.15 (2C).")

# ---------- 3 · Concept of Operations ----------
s = slides[2]
set_text(find_text(s, "Sources"),
         "Sources   Arianespace VV24/VV25, ESA Operations (2C LEOP via Troll), DLR (1C/1A phasing), "
         "ISSFD 2015 paper 137, SpaceOps 2025 papers 226 and 617, Airbus GNC 2014")

# ---------- 7 · Sentinel-1 Manoeuvres ----------
s = slides[6]
S = shapes(s)
fig_h = 4.4
fig_w = fig_h * 755 / 1180
replace_picture(s, S[107], FIG / "B_s1b_itici_yonleri_spaceops226_fig3.png",
                7.75 + (4.98 - fig_w) / 2, 1.5, h=fig_h, name="S1B thruster directions")
set_text(S[1006], "1A: Torques about 10× higher than predicted in early tests led to a 30 s limit "
                  "on in-plane firings.")
set_text(S[1008], "Sentinel-1B thruster directions in pitch-slew disposal attitude · SpaceOps 2025, Fig. 3")
set_text(S[113], "ISSFD 2015 papers 134 (p.8: torque, −21 % thrust) and 137 (p.14: 30 s limit). "
                 "SpaceOps 2025, papers 226 and 360.")

# ---------- 8 · Sentinel-2C Δv ----------
s = slides[7]
set_text(find_text(s, "Planning case"),
         "Planned, not flown: two plane corrections (0.57 m/s) → ≈7.97 m/s.")

# ---------- 9 · Propellant Budget ----------
s = slides[8]
S = shapes(s)
old_pic = next(sh for sh in s.shapes if sh.shape_type == 13)
replace_picture(s, old_pic, CHART, 7.9, 1.71, w=4.68, name="S2C hydrazine waterfall")
set_text(S[1016], "82.1 kg")
set_text(find_text(s, "Estimate, not telemetry"),
         "Estimate, not telemetry  ·  plan 47 kg vs ideal 38 kg (73 m/s)")

prs.save(OUT)
print("wrote", OUT)
