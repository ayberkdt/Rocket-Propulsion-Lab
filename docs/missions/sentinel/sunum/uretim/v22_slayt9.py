"""v22: V21 aynen korunur; yalnız slayt 9 (Propellant Budget) V21'in kendi şekil stilleriyle yeniden kurulur.

Girdi: surumler/..._v21.pptx + grafik/s2c_yakit_v22.png  ·  Çıktı: surumler/..._v22.pptx
"""
import copy
import sys
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Inches, Pt

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "surumler" / "Sentinel_1_2_Propulsion_ITU_v21.pptx"
CHART = HERE / "grafik" / "s2c_yakit_v22.png"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "surumler" / "Sentinel_1_2_Propulsion_ITU_v22.pptx"

prs = Presentation(SRC)
s4, s7, s9 = prs.slides[3], prs.slides[6], prs.slides[8]
by_id = lambda slide: {sh.shape_id: sh for sh in slide.shapes}
S4, S7, S9 = by_id(s4), by_id(s7), by_id(s9)
assert S9[145].text_frame.text == "Propellant Budget"

next_id = max(sh.shape_id for sh in s9.shapes) + 1


def set_text(shape, text):
    runs = shape.text_frame.paragraphs[0].runs
    runs[0].text = text
    for r in runs[1:]:
        r._r.getparent().remove(r._r)
    for p in shape.text_frame.paragraphs[1:]:
        p._p.getparent().remove(p._p)


def set_runs(shape, parts):
    """parts: list of (text, None|'sub'|'sup'); copies the first run's formatting."""
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
        if kind:
            rpr.set("baseline", "-25000" if kind == "sub" else "30000")
        elif "baseline" in rpr.attrib:
            del rpr.attrib["baseline"]
    # keep endParaRPr last
    end = p._p.find("{*}endParaRPr")
    if end is not None:
        p._p.remove(end); p._p.append(end)


def clone(src, x, y, w, h, text):
    global next_id
    el = copy.deepcopy(src._element)
    el.find(".//{*}cNvPr").set("id", str(next_id))
    next_id += 1
    s9.shapes._spTree.append(el)
    shp = s9.shapes[-1]
    shp.left, shp.top, shp.width, shp.height = Inches(x), Inches(y), Inches(w), Inches(h)
    set_text(shp, text)
    return shp


# remove the old table and the two explanatory lines
for sid in (155, 1011, 1012):
    el = S9[sid]._element
    el.getparent().remove(el)

set_text(S9[1010], "Sentinel-2C hydrazine estimate")
S9[1010].top, S9[1010].width = Inches(1.57), Inches(7.3)

clone(S7[1001], 0.59, 2.1, 3.04, 0.83, "3.9 kg")
clone(S7[1003], 0.64, 2.95, 2.92, 0.33, "Flown so far (7.40 m/s)")
clone(S7[1002], 3.57, 2.1, 3.91, 0.83, "82 kg")
clone(S7[1004], 3.61, 2.95, 4.01, 0.38, "Left after 47 kg disposal")

clone(S7[1005], 0.59, 3.88, 7.14, 0.38, "Rocket equation")
eq = clone(S4[1022], 0.59, 4.36, 7.14, 0.47, "x")
set_runs(eq, [("Δv = I", None), ("sp", "sub"), (" g", None), ("0", "sub"), (" ln(m", None), ("0", "sub"),
              (" / m", None), ("f", "sub"), (")", None)])
inp = clone(S7[1006], 0.59, 4.95, 6.93, 0.4, "x")
set_runs(inp, [("m", None), ("0", "sub"), (" = 1,143 kg,  I", None), ("sp", "sub"), (" = 220 s (assumed)", None)])
clone(S7[1007], 0.59, 5.45, 6.93, 0.4, "Propellant sized for 12 years; design life 7 years")

card = s9.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(7.83), Inches(1.78), Inches(4.89), Inches(3.17))
card.adjustments[0] = 0.04
card.fill.solid(); card.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
card.line.color.rgb = RGBColor(0xE0, 0xE7, 0xEF); card.line.width = Pt(1)
card.shadow.inherit = False
card.name = "Chart card"
pic = s9.shapes.add_picture(str(CHART), Inches(7.98), Inches(1.93), width=Inches(4.59))
pic.name = "S2C hydrazine waterfall"
clone(S7[1008], 7.88, 5.07, 4.74, 0.52, "Estimate, not telemetry  ·  disposal: A/B plan, 73 m/s")

set_text(S9[1013], "Sources   Arianespace VV24, ESA fuelling report, SpaceOps 2025 paper 617, "
                   "Airbus GNC 2014, ESA SP-1322/2")

# Slayt 9 okur metni ozel/okur_metni_v22.md dosyasında; desteye not yazılmaz.

prs.save(OUT)
print("wrote", OUT, "slides:", len(prs.slides))
