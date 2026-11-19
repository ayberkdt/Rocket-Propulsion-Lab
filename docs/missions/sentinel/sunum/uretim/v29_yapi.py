"""v29: deste itki merkezli ve tek şablonla baştan kurulur (10 slayt).

Kullanıcı geri bildirimi (2026-11-18): ayrıntı ana bölümlerin önüne geçiyordu; propulsion için düzgün bir
tablo yoktu; yapı "AI sunumu" gibiydi (KPI kartları, iddia/kanıt satırları, dağınık dipnotlar).
Şablon: bant başlık · (varsa) alt başlık · tablo ya da tek blok · "Table n." altyazısı · tek kaynak satırı.

 1 Mission Overview              Table 1
 2 Spacecraft and Launch         iki fotoğraf
 3 Propulsion Architecture       Table 2  (ana slayt)
 4 Feed System and Layout        akış şeridi (v28'den) + iki S1B figürü + 4 satır
 5 Thruster Performance          Table 3
 6 Concept of Operations         evre şeridi + 3 satır
 7 Flight Experience             Table 4
 8 ΔV Budget                     Table 5
 9 Propellant Budget             grafik + roket denklemi
10 Conclusions

Mimari tablosunun hücreleri: SENTINEL-1C.md D1 (14 x 1 N, OCT +X/-X/-Y, PEPT-590GB 204 L,
27.4 -> 5.5 bar, tek tank çıkarımı), SENTINEL-2_SERISI.md D1-D3 ve SENTINEL-2C.md D1 (2 x 4 eğik
CHT1N-6, He blow-down, modül arayüz halkası içinde, dört itici birlikte), ESA dolum haberleri.
Girdi: surumler/..._v28.pptx (şablon, görseller, akış şeridi) · Çıktı: surumler/..._v29.pptx
"""
import copy
import io
import math
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt
from lxml import etree

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v28.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v29.pptx"
CHART = HERE / "grafik" / "s2c_yakit_v29.png"

FONT = "Helvetica Neue"
NAVY, INK, SRC_INK, TEAL = "092E52", "131A22", "082D55", "066764"
HEAD_FILL, BAND, LINE = "1C5284", "EDF4FA", "AABCCD"
X0, W = 0.59, 12.15

# ---------------------------------------------------------------- hesaplar (v27/v28 ile aynı kaynaklar)
G0, ISP = 9.80665, 220.0


def mp(m0, dv, isp=ISP):
    return m0 * (1 - math.exp(-dv / (isp * G0)))


YIL = 2024.21 - 2015.50                                       # 2A: 2 Tem 2015 - Mar 2024 (#617 s.5)
DV_ACQ = 6.33 + 0.012 + 0.0146 + math.hypot(0.9344, 0.0537) + 0.0025 + 0.0207 + 0.0206 + 0.0309 + 0.0029 + 0.0283
DV_OPS = (64 * 0.0156 + 24 * 2 * 0.04) / YIL * 12
m = 1143.0
KG_ACQ = mp(m, DV_ACQ); m -= KG_ACQ
KG_OPS = mp(m, 64 * 0.0156 / YIL * 12); m -= KG_OPS
k2 = mp(m, 24 * 2 * 0.04 / YIL * 12); m -= k2; KG_OPS += k2
EOL_IDEAL = mp(m, 73.0)
DV_TOT = DV_ACQ + DV_OPS + 73.0
KG_TOT = KG_ACQ + KG_OPS + 47.0
SHARE = 73.0 / DV_TOT
OPS_PER_YEAR = DV_OPS / 12
S1A_IDEAL = mp(2157.0, 4.1, 214.0)
S1A_IDEAL += mp(2157.0 - S1A_IDEAL, 2.2, 214.0)
PERIOD = 2 * math.pi * math.sqrt((6378.137e3 + 786e3) ** 3 / 3.986004418e14)
ARC = 2 * math.pi * 314.1 / PERIOD
LOSS = 1 - math.sin(ARC / 2) / (ARC / 2)
assert abs(DV_ACQ - 7.398) < 0.001 and abs(DV_OPS - 4.02) < 0.01 and abs(KG_ACQ - 3.91) < 0.01
assert abs(KG_OPS - 2.12) < 0.01 and abs(DV_TOT - 84.42) < 0.02 and abs(SHARE - 0.865) < 0.001
assert abs(KG_TOT - 53.03) < 0.02 and abs(EOL_IDEAL - 37.9) < 0.1 and abs(S1A_IDEAL - 6.47) < 0.01
assert abs(LOSS - 0.00445) < 0.0001 and abs(OPS_PER_YEAR - 0.335) < 0.001

# ---------------------------------------------------------------- şablondan alınacaklar
prs = Presentation(SRC)
old = list(prs.slides)
assert len(old) == 11


def by_id(slide):
    return {sh.shape_id: sh for sh in slide.shapes}


LAYOUT = old[4].slide_layout                                   # "Başlıklı İçerik": bant + İTÜ logosu
TITLE_XML = copy.deepcopy(by_id(old[4])[115]._element)        # beyaz başlık metni (bant üstünde)
AUTHOR_XML = copy.deepcopy(by_id(old[4])[116]._element)
FLOW_XML = [copy.deepcopy(by_id(old[4])[i]._element) for i in range(300, 309)]   # akış şeridi
IMG = {
    "s1c": by_id(old[1])[71].image.blob,
    "s2c": by_id(old[1])[72].image.blob,
    "layout": by_id(old[4])[309].image.blob,
    "rafael": by_id(old[5])[126].image.blob,
    "ag": by_id(old[5])[132].image.blob,
    "dirs": by_id(old[6])[1009].image.blob,
}

sld_ids = prs.slides._sldIdLst
for sid in list(sld_ids):                                      # eski slaytları çıkar
    prs.part.drop_rel(sid.rId)
    sld_ids.remove(sid)


# ---------------------------------------------------------------- şablon yardımcıları
def new_slide(title):
    s = prs.slides.add_slide(LAYOUT)
    for ph in list(s.placeholders):
        ph._element.getparent().remove(ph._element)
    for xml in (TITLE_XML, AUTHOR_XML):
        s.shapes._spTree.append(copy.deepcopy(xml))
    t = s.shapes[0]
    t.text_frame.paragraphs[0].runs[0].text = title
    return s


def text(s, x, y, w, h, lines, size=16.5, bold=False, color=NAVY, align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.TOP, space=0):
    """lines: str ya da [(metin, punto, kalın, renk)] / [str]; her öğe bir paragraf."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.05)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    items = [lines] if isinstance(lines, str) else lines
    for k, item in enumerate(items):
        if isinstance(item, str):
            item = (item, size, bold, color)
        txt, sz, bd, col = item
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        p.alignment = align
        if k and space:
            p.space_before = Pt(space)
        for j, part in enumerate(txt.split("|")):              # "|" = satır kırılımı
            if j:
                p.add_line_break()
            r = p.add_run()
            r.text = part
            r.font.size, r.font.bold, r.font.name = Pt(sz), bd, FONT
            r.font.color.rgb = RGBColor.from_string(col)
    return tb


def rect(s, x, y, w, h, fill, kind=MSO_SHAPE.RECTANGLE):
    shp = s.shapes.add_shape(kind, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = RGBColor.from_string(fill)
    shp.line.fill.background()
    shp.shadow.inherit = False
    return shp


def _border(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    for tag in ("a:lnL", "a:lnR", "a:lnT", "a:lnB"):
        for old_ln in tcPr.findall(qn(tag)):
            tcPr.remove(old_ln)
        ln = etree.SubElement(tcPr, qn(tag), w="9525")
        fill = etree.SubElement(ln, qn("a:solidFill"))
        etree.SubElement(fill, qn("a:srgbClr"), val=LINE)
    # OOXML sırası: kenarlar dolgudan önce gelir
    fills = tcPr.findall(qn("a:solidFill"))
    for f in fills:
        tcPr.remove(f)
        tcPr.append(f)


def table(s, rows, widths, y, row_h=0.42, size=13.5, head_size=14.25, num_cols=(), head_h=None):
    nr, nc = len(rows), len(rows[0])
    assert abs(sum(widths) - W) < 1e-6, sum(widths)
    gf = s.shapes.add_table(nr, nc, Inches(X0), Inches(y), Inches(W), Inches(row_h * nr))
    tb = gf.table
    tb.first_row = False
    tb.horz_banding = False
    for c, w in enumerate(widths):
        tb.columns[c].width = Inches(w)
    for r, row in enumerate(rows):
        tb.rows[r].height = Inches(head_h if (r == 0 and head_h) else row_h)
        for c, val in enumerate(row):
            cell = tb.cell(r, c)
            cell.margin_left = cell.margin_right = Inches(0.09)
            cell.margin_top = cell.margin_bottom = Inches(0.04)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            head = r == 0
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string(HEAD_FILL if head else ("FFFFFF" if r % 2 else BAND))
            tf = cell.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.RIGHT if (c in num_cols and not head) else PP_ALIGN.LEFT
            if head and c in num_cols:
                p.alignment = PP_ALIGN.RIGHT
            rr = p.add_run()
            rr.text = val
            rr.font.name = FONT
            rr.font.size = Pt(head_size if head else size)
            rr.font.bold = head or c == 0
            rr.font.color.rgb = RGBColor.from_string("FFFFFF" if head else INK)
            _border(cell)
    return gf


def caption(s, y, n, body):
    tb = s.shapes.add_textbox(Inches(X0), Inches(y), Inches(W), Inches(0.32))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    for txt, bold in ((f"Table {n}. ", True), (body, False)):
        r = p.add_run()
        r.text = txt
        r.font.size, r.font.bold, r.font.name = Pt(12), bold, FONT
        r.font.color.rgb = RGBColor.from_string(SRC_INK)


def sources(s, body, y=6.78):
    text(s, X0, y, W, 0.30, [(f"Sources: {body}", 10, False, SRC_INK)])


def picture(s, key, x, y, w=None, h=None):
    return s.shapes.add_picture(io.BytesIO(IMG[key]), Inches(x), Inches(y),
                                Inches(w) if w else None, Inches(h) if h else None)


def subtitle(s, txt, y=1.22):
    text(s, X0, y, W, 0.45, [(txt, 20, True, INK)])


# ================================================================= 1 · Mission Overview
s = new_slide("Mission Overview")
subtitle(s, "Copernicus Sentinel-1 and Sentinel-2")
table(s, [
    ["Parameter", "Sentinel-1", "Sentinel-2"],
    ["Purpose", "All-weather radar mapping of land, ice and oceans", "Optical imaging of land and coastal waters"],
    ["Payload", "C-band synthetic aperture radar", "13-band multispectral imager"],
    ["Prime contractor", "Thales Alenia Space", "Airbus Defence and Space"],
    ["Orbit", "693 km, Sun-synchronous", "786 km, Sun-synchronous"],
    ["Revisit with two satellites", "6 days", "5 days"],
    ["Design life (C unit)", "7.25 years", "7 years"],
    ["In service, Nov 2026", "1C and 1D", "2B and 2C; 2A extended to end-2027"],
], [2.95, 4.60, 4.60], y=1.85, row_h=0.52)
caption(s, 6.10, 1, "Mission parameters of the two Copernicus series, operated by ESA.")
sources(s, "ESA; Copernicus (9 Oct 2026); Arianespace launch kits VV24 and VV25.")

# ================================================================= 2 · Spacecraft and Launch
s = new_slide("Spacecraft and Launch")
for k, (key, cap, name, l1, l2) in enumerate([
    ("s1c", "Sentinel-1C during radio-frequency testing", "Sentinel-1C",
     "Vega-C, 5 December 2024", "Launch mass 2,186 kg, released at ~700 km"),
    ("s2c", "Sentinel-2C fuelled with toxic hydrazine in protective suits", "Sentinel-2C",
     "Vega, 5 September 2024", "Launch mass 1,143 kg, released at ~775 km"),
]):
    x = X0 if k == 0 else 6.96
    picture(s, key, x, 1.30, 5.50, 3.60)
    text(s, x, 4.95, 5.50, 0.30, [(cap, 11.25, False, SRC_INK)], align=PP_ALIGN.CENTER)
    text(s, x, 5.40, 5.78, 1.10, [(name, 21, True, NAVY), (l1, 17.25, False, NAVY), (l2, 17.25, False, NAVY)])
sources(s, "Arianespace VV24, VV25; SpaceOps 2025 paper 617. Photos: Thales Alenia Space/ImagIN; "
           "ESA–CNES–Arianespace.")

# ================================================================= 3 · Propulsion Architecture
s = new_slide("Propulsion Architecture")
table(s, [
    ["Parameter", "Sentinel-1", "Sentinel-2"],
    ["Propulsion type", "Chemical monopropellant, catalytic", "Chemical monopropellant, catalytic"],
    ["Propellant", "Hydrazine (N₂H₄)", "Hydrazine (N₂H₄)"],
    ["Thruster configuration", "2 × 7 thrusters of 1 N (3 orbit, 4 attitude)", "2 × 4 canted thrusters of 1 N, cold-redundant"],
    ["Thruster model", "Rafael 1 N", "CHT1N-6 (ArianeGroup)"],
    ["Propellant storage", "1 titanium tank, 204 L*", "1 blowdown tank, model not published"],
    ["Pressurization", "Blowdown, 27.4 → 5.5 bar", "Blowdown, helium, 21 bar at fuelling"],
    ["Propellant load", "154 kg (1C)", "133 kg (2C)"],
    ["Thrust control", "Pulsed and steady; in-plane burns ≤ 30 s", "Four thrusters fire together in orbit burns"],
    ["Physical location", "Propulsion module; orbit thrusters on ±X, −Y", "Module inside the launcher interface ring"],
    ["Operational role", "Orbit control, attitude, avoidance, disposal", "Orbit control, avoidance, safe mode, disposal"],
], [2.95, 4.60, 4.60], y=1.30, row_h=0.44)
caption(s, 6.22, 2, "Propulsion architecture of the C units. * One tank inferred: the 154 kg fill equals its "
                    "nominal capacity.")
sources(s, "SpaceOps 2025 papers 360, 617; ISSFD 2015 paper 134; Rafael catalogue pp. 9, 41; "
           "Airbus GNC 2014; ESA SP-1322; ESA fuelling reports.")

# ================================================================= 4 · Feed System and Layout
s = new_slide("Feed System and Layout")
for xml in FLOW_XML:                                            # v28'deki akış şeridi, değişmeden
    s.shapes._spTree.append(copy.deepcopy(xml))
picture(s, "layout", X0, 2.30, 2.45, 3.90)
picture(s, "dirs", 3.25, 2.30, 2.50, 3.90)
text(s, X0, 6.24, 5.30, 0.45, [("Sentinel-1B thruster layout and firing directions (SpaceOps 2025 "
                                "papers 360 and 226)", 11.25, False, SRC_INK)])
text(s, 6.35, 2.45, 6.39, 3.80, [
    ("Blowdown: no regulator. Tank pressure, thrust and Isp fall as hydrazine is used.", 17.25, False, NAVY),
    ("Isolation (latch) valves separate the tank from two thruster branches; one branch is active.",
     17.25, False, NAVY),
    ("The catalyst bed decomposes hydrazine into hot gas; heaters warm it before firing.", 17.25, False, NAVY),
], space=14)
sources(s, "SpaceOps 2025 papers 226, 360; ESA (2C fuelling, latch valves); Airbus GNC 2014.")

# ================================================================= 5 · Thruster Performance
s = new_slide("Thruster Performance")
picture(s, "rafael", X0, 1.25, 5.20, 2.55)
picture(s, "ag", 8.85, 1.25, 2.29, 2.55)
text(s, X0, 3.82, 5.20, 0.28, [("Rafael 1 N thruster (drawing)", 11.25, False, SRC_INK)], align=PP_ALIGN.CENTER)
text(s, 7.50, 3.82, 5.00, 0.28, [("ArianeGroup 1 N thruster (integrated)", 11.25, False, SRC_INK)],
     align=PP_ALIGN.CENTER)
table(s, [
    ["Parameter", "Rafael 1 N (Sentinel-1)", "ArianeGroup 1 N (Sentinel-2)"],
    ["Inlet pressure", "5.5–24.5 bar", "5.5–22 bar"],
    ["Thrust", "1.0–1.4 N", "0.32–1.1 N"],
    ["Specific impulse", "205–214 s", "200–223 s"],
    ["Minimum impulse bit", "0.008 N·s", "0.010–0.043 N·s"],
    ["Flight heritage", "Sentinel-1 since 2014 (TRL 9)", "Sentinel-2 since 2015 (TRL 9)"],
], [2.95, 4.60, 4.60], y=4.22, row_h=0.36)
caption(s, 6.42, 3, "Catalogue performance of the two thruster families.")
sources(s, "Rafael catalogue p. 46; ArianeGroup 1 N brochure p. 6; TRL per ISO 16290.", y=6.80)

# ================================================================= 6 · Concept of Operations
s = new_slide("Concept of Operations")
LW = 1.55
PW = (W - LW) / 5
for k, (name, span) in enumerate(zip(
        ["Insertion", "LEOP", "Orbit acquisition", "Routine operations", "Disposal"],
        ["separation", "3 days", "2–4 months", "7–12 years", "end of life"])):
    kind = MSO_SHAPE.PENTAGON if k == 0 else MSO_SHAPE.CHEVRON
    ch = rect(s, X0 + LW + k * PW, 1.30, PW + 0.10 if k < 4 else PW, 0.58, NAVY, kind)
    ch.adjustments[0] = 0.28
    tf = ch.text_frame
    tf.margin_left = tf.margin_right = Inches(0.22)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = name
    r.font.size, r.font.bold, r.font.name = Pt(15), True, FONT
    r.font.color.rgb = RGBColor.from_string("FFFFFF")
    text(s, X0 + LW + k * PW, 1.92, PW, 0.30, [(span, 13.5, False, SRC_INK)], align=PP_ALIGN.CENTER)
GRID = [
    ("Propulsion", BAND, 1.05, [
        "None: launcher|upper stage", "No burns planned;|avoidance if needed",
        "Correct injection|errors", "Ground track,|avoidance, inclination", "Lower the orbit,|then passivate"]),
    ("Example", None, 1.05, [
        "2C released 10.6 km|below its orbit", "1A: one avoidance|burn on day 2",
        "1A: ~500 burns|in four months", "2A: 64 trims and 24|avoidance pairs", "1B lowered to|~580 × 610 km"]),
]
y = 2.42
for label, band, h, cells in GRID:
    if band:
        rect(s, X0, y, W, h, band)
    text(s, X0, y + 0.08, LW, h, [(label, 15, True, NAVY)])
    for k, c in enumerate(cells):
        text(s, X0 + LW + k * PW, y + 0.08, PW, h, [(c, 15, False, INK)])
    y += h + 0.08
rect(s, X0, y + 0.05, W, 0.025, NAVY)
y += 0.15
text(s, X0, y, LW, 0.75, [("ΔV, 2C", 15, True, NAVY)])
for k, (num, qual) in enumerate([("—", ""), ("—", ""), ("7.4 m/s", "flown"), (f"{DV_OPS:.1f} m/s", "projected, 12 yr"),
                                  ("73 m/s", "A/B design")]):
    items = [(num, 17.25, True, NAVY)] + ([(qual, 13.5, False, SRC_INK)] if qual else [])
    text(s, X0 + LW + k * PW, y, PW, 0.75, items)
sources(s, "ISSFD 2015 papers 134, 137; SpaceOps 2025 papers 226, 617; Airbus GNC 2014.")

# ================================================================= 7 · Flight Experience
s = new_slide("Flight Experience")
subtitle(s, "Sentinel-1 propulsion in flight")
table(s, [
    ["Event", "What happened", "Consequence"],
    ["1A orbit acquisition, 2014", "Plume impingement: in-plane thrust 21–30 % low, torque ~10× predicted",
     f"Burns ≤ 30 s: ~500 burns, 8.4 kg ({S1A_IDEAL:.1f} ideal)"],
    ["1B disposal, 2023–24", "Reaction-wheel overspeed and several aborted burns",
     "Lowering stopped at ~580 × 610 km"],
], [2.95, 4.70, 4.50], y=1.85, row_h=0.95, size=15, head_size=15, head_h=0.52)
caption(s, 4.42, 4, "Propulsion problems reported in flight and their effect on operations.")
sources(s, "ISSFD 2015 papers 134, 137; SpaceOps 2025 papers 226, 360. Ideal value: rocket equation, "
           "6.3 m/s, m₀ 2,157 kg, Isp 214 s.")

# ================================================================= 8 · ΔV Budget
s = new_slide("ΔV Budget")
subtitle(s, "Sentinel-2C, illustrative 12-year budget")
table(s, [
    ["Phase", "Δv (m/s)", "Propellant (kg)", "Status"],
    ["Acquisition, 2024", f"{DV_ACQ:.1f}", f"{KG_ACQ:.1f}", "Flown"],
    ["Operations, 12 years", f"{DV_OPS:.1f}", f"{KG_OPS:.1f}", "Projected from Sentinel-2A, 2015–2024"],
    ["Inclination control", "—", "—", "Not published"],
    ["Disposal, 786 → 653 km", "73.0", "47.0", f"A/B design value ({EOL_IDEAL:.1f} kg ideal)"],
    ["Total", f"≈ {DV_TOT:.1f}", f"≈ {KG_TOT:.1f}", "Illustrative, not a verified plan"],
], [3.35, 1.70, 2.20, 4.90], y=1.85, row_h=0.62, size=15, head_size=15, num_cols=(1, 2))
caption(s, 5.70, 5, f"Disposal is {SHARE * 100:.0f} % of the illustrative Δv. Propellant: rocket equation, "
                    "m₀ = 1,143 kg, Isp = 220 s.")
sources(s, "SpaceOps 2025 paper 617 pp. 5, 15–17; Airbus GNC 2014 p. 24.")

# ================================================================= 9 · Propellant Budget
s = new_slide("Propellant Budget")
subtitle(s, "Sentinel-2C hydrazine, scenario", y=1.30)
eq = text(s, X0, 2.00, 6.90, 0.55, "")
p = eq.text_frame.paragraphs[0]
for part, sub in (("Δv = I", False), ("sp", True), (" g", False), ("0", True), (" ln(m", False), ("0", True),
                  (" / m", False), ("f", True), (")", False)):
    r = p.add_run()
    r.text = part
    r.font.size, r.font.bold, r.font.name = Pt(24), True, FONT
    r.font.color.rgb = RGBColor.from_string(TEAL)
    if sub:
        r.font._rPr.set("baseline", "-25000")
text(s, X0, 2.75, 6.90, 3.60, [
    "m₀ = 1,143 kg, Isp = 220 s (assumed)",
    f"Mass flow of one 1 N thruster: {1000 / (ISP * G0):.2f} g/s",
    f"Longest burn 314 s, {314.1 / PERIOD * 100:.1f} % of an orbit",
    f"Finite-burn penalty ≈ {LOSS * 100:.2f} % (node-centred estimate)",
    "Propellant sized for 12 years; design life 7 years",
], size=17.25, space=12)
card = rect(s, 7.75, 1.40, 4.98, 4.55, "FFFFFF", MSO_SHAPE.ROUNDED_RECTANGLE)
card.adjustments[0] = 0.04
card.line.color.rgb = RGBColor.from_string("D9E1EA")
card.line.width = Pt(0.75)
s.shapes.add_picture(str(CHART), Inches(7.90), Inches(1.55), Inches(4.68), Inches(3.97))
text(s, 7.75, 6.00, 4.98, 0.32, [("Scenario, not telemetry.", 12, False, SRC_INK)], align=PP_ALIGN.CENTER)
sources(s, "ESA fuelling report; Arianespace VV24; SpaceOps 2025 paper 617; Airbus GNC 2014; ESA SP-1322/2.")

# ================================================================= 10 · Conclusions
s = new_slide("Conclusions")
text(s, X0, 1.55, W, 4.90, [
    "Both series use 1 N hydrazine thrusters in a blowdown system, flown since 2014.",
    f"Routine Δv is small: {DV_ACQ:.1f} m/s for 2C's acquisition, ~{OPS_PER_YEAR:.1f} m/s per year after.",
    f"Disposal is the largest item, {SHARE * 100:.0f} % of the illustrative Δv budget.",
    "Plume impingement limited Sentinel-1A's in-plane burns to 30 s.",
    "Propellant for 12 years exceeds the 7-year design life; 1A flew 12 years.",
], size=20, space=22)

assert len(prs.slides) == 10
prs.save(OUT)
print("wrote", OUT)
print(f"ΔV acq {DV_ACQ:.3f} ops {DV_OPS:.3f} total {DV_TOT:.2f} share {SHARE:.3f}; "
      f"kg acq {KG_ACQ:.3f} ops {KG_OPS:.3f} total {KG_TOT:.2f}; EOL ideal {EOL_IDEAL:.2f}")
