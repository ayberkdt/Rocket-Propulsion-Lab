"""v26: ödevin eksik maddeleri, yeni sayfa açmadan birer satır.

D9 güvenlik  : slayt 2, dolum fotoğrafı altyazısı (toksik hidrazin, SCAPE giysisi).
D9 düzenleme : slayt 3, 1B bertaraf hücresi (25 yıl kuralı).
D10 TRL      : slayt 6, kaynak satırının başı.
C2/D6 (S1)   : slayt 7, S1A için roket denklemi: 6.3 m/s -> ideal 6.7 kg, uçuşta 8.4 kg.
D5 ve D8     : slayt 9, kütle debisi / yanma süresi ve anlık yanma kabulü.
Girdi: surumler/..._v25.pptx  ·  Çıktı: surumler/..._v26.pptx
Hesaplar [H]: g0 = 9.80665 m/s², Isp = 220 s varsayımı.
"""
import copy
import math
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v25.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v26.pptx"

G0, ISP = 9.80665, 220.0
VE = ISP * G0
# S1A: ISSFD 2015 #134 Table 2 (4.1 + 2.2 m/s), SP-1322/1 p.24 tasarım kütlesi 2300 kg
S1A_DV, S1A_M0 = 4.1 + 2.2, 2300.0
S1A_IDEAL = S1A_M0 * (1 - math.exp(-S1A_DV / VE))
# S2C: 7.40 m/s, 1143 kg (VV24); 1 N için debi ve eşdeğer yanma süresi
S2C_MP = 1143.0 * (1 - math.exp(-(6.33 + 0.012 + 1.056442) / VE))
MDOT = 1.0 / VE
BURN_H = S2C_MP / MDOT / 3600
PERIOD_MIN = 2 * math.pi * math.sqrt((6378.137e3 + 786e3) ** 3 / 3.986004418e14) / 60
assert abs(S1A_IDEAL - 6.7) < 0.05 and abs(MDOT * 1000 - 0.46) < 0.01 and abs(BURN_H - 2.3) < 0.05
assert abs(314.1 / (PERIOD_MIN * 60) - 0.052) < 0.002

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


def place(sh, x=None, y=None, w=None, h=None):
    if x is not None: sh.left = Inches(x)
    if y is not None: sh.top = Inches(y)
    if w is not None: sh.width = Inches(w)
    if h is not None: sh.height = Inches(h)


def clone(slide, src, x, y, w, h, text):
    el = copy.deepcopy(src._element)
    el.find(".//{*}cNvPr").set("id", str(max(sh.shape_id for sh in slide.shapes) + 1))
    slide.shapes._spTree.append(el)
    shp = slide.shapes[-1]
    place(shp, x, y, w, h)
    set_text(shp, text)
    return shp


def set_runs(shape, parts):
    """parts: [(metin, None|'sub')]; ilk koşunun biçimini kopyalar."""
    p = shape.text_frame.paragraphs[0]
    first = p.runs[0]._r
    for r in p.runs[1:]:
        r._r.getparent().remove(r._r)
    for k, (txt, kind) in enumerate(parts):
        r = first if k == 0 else copy.deepcopy(first)
        if k:
            p._p.append(r)
        r.find("{*}t").text = txt
        rpr = r.find("{*}rPr")
        if kind == "sub":
            rpr.set("baseline", "-25000")
        elif "baseline" in rpr.attrib:
            del rpr.attrib["baseline"]
    end = p._p.find("{*}endParaRPr")
    if end is not None:
        p._p.remove(end)
        p._p.append(end)


# ---------- 2 · D9 güvenlik ----------
s = slides[1]
set_text(find_text(s, "Sentinel-2C during fuelling"), "Sentinel-2C fuelling · toxic hydrazine, SCAPE suits")

# ---------- 3 · D9 düzenleme ----------
s = slides[2]
set_text(find_text(s, "ESA lowered 1B"), "1B lowered to ~580 × 610 km for 25-year re-entry, then passivated.")

# ---------- 6 · D10 TRL ----------
s = slides[5]
set_text(find_text(s, "Supplier product specifications"),
         "TRL 9: both families flown on Sentinel-1/2 since 2014 (ISO 16290). "
         "Specs: Rafael p. 46, ArianeGroup brochure p. 6. Photo: ArianeGroup.")

# ---------- 7 · S1A roket denklemi ----------
s = slides[6]
S = shapes(s)
clone(s, S[1006], 0.59, 3.45, 6.93, 0.42,
      f"Needed {S1A_DV:.1f} m/s → ideal {S1A_IDEAL:.1f} kg [H]; flown 8.4 kg (plume losses)")
set_text(S[113], "ISSFD 2015 papers 134 (p.8: Table 2 Δv, torque, −21 % thrust) and 137 (p.14: 30 s limit). "
                 "SpaceOps 2025, papers 226 and 360. [H] m₀ 2,300 kg (SP-1322/1), Isp 220 s.")

# ---------- 9 · D5 ve D8 ----------
s = slides[8]
S = shapes(s)
for sid, y in ((1018, 3.75), (1019, 4.2), (1020, 4.78)):
    place(S[sid], y=y)
mdot = clone(s, S[1020], 0.59, 5.24, 6.93, 0.4, "x")
set_runs(mdot, [("ṁ = F / (I", None), ("sp", "sub"), (" g", None), ("0", "sub"),
                (f") = {MDOT * 1000:.2f} g/s for 1 N → ≈ {BURN_H:.1f} h of firing", None)])
clone(s, S[1020], 0.59, 5.7, 6.93, 0.4,
      f"Longest burn 314 s ≈ {314.1 / (PERIOD_MIN * 60) * 100:.0f} % of a {PERIOD_MIN:.0f}-min orbit → impulsive")
place(S[1021], y=6.16)

prs.save(OUT)
print("wrote", OUT)
print(f"S1A ideal {S1A_IDEAL:.2f} kg, mdot {MDOT*1000:.3f} g/s, burn {BURN_H:.2f} h, period {PERIOD_MIN:.1f} min")
