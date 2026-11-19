"""v28: ikinci geri bildirim ve CONOPS slaytının yeniden tasarımı.

Slayt 3  : iki küçük makale figürü (sınıfta okunmuyordu) ve 5 x 2 yazı ızgarası yerine evre şeridi
           (5 chevron, süre satırı) + 4 satır: itkinin işi, Sentinel-1 örneği, Sentinel-2 örneği, 2C ΔV.
           Kaynaklar: ISSFD 2015 #137 (üç günlük LEOP, ikinci gün kaçınma yanması), #134 (Nis-Ağu 2014,
           ~500 manevra, 8.4 kg), SpaceOps 2025 #617 (2C 10.6 km alçak, 29 Eki ince edinim, 2A 64 OCM ve
           24 CAM çifti), #226 (1B ~580 x 610 km), Airbus GNC 2014 (786 -> 653 km, 73 m/s).
Slayt 7  : 7.6-8.3 kg aralığı "simplified estimate".
Slayt 9  : "Illustrative 12-year budget"; Basis sütunu durumla başlar (flown / planned / projected /
           A/B design); eğim satırı "Additional inclination control (12 yr)", 6 Kasım yanması satır 1'de.
Slayt 10 : yatay şelale, son çubuk "Unallocated (scenario)", EOL "A/B design"; 3.4 N "average
           effective force"; "impulsive model holds" yerine sonlu yanma cezası ve varsayımı.
Slayt 11 : %86 "illustrative 12-year Δv"ye göre.
Girdi: surumler/..._v27.pptx  ·  Çıktı: surumler/..._v28.pptx  ·  Grafik: grafik_s2c_yakit.py v28
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v27.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v28.pptx"
CHART = HERE / "grafik" / "s2c_yakit_v28.png"

NAVY, TEXT, TEAL = "092E52", "092E52", "066764"
BAND_1, BAND_2 = "EDF4FA", "F4F9F7"             # v27 satır zeminleri
FONT = "Helvetica Neue"

prs = Presentation(SRC)
slides = list(prs.slides)
assert len(slides) == 11


def shapes(slide):
    return {sh.shape_id: sh for sh in slide.shapes}


def set_text(shape, text):
    p = shape.text_frame.paragraphs[0]
    p.runs[0].text = text
    for r in p.runs[1:]:
        r._r.getparent().remove(r._r)
    for extra in shape.text_frame.paragraphs[1:]:
        extra._p.getparent().remove(extra._p)


def find_text(slide, start):
    hits = [sh for sh in slide.shapes if sh.has_text_frame and sh.text_frame.text.startswith(start)]
    assert len(hits) == 1, (start, len(hits))
    return hits[0]


def place(sh, x=None, y=None, w=None, h=None):
    if x is not None: sh.left = Inches(x)
    if y is not None: sh.top = Inches(y)
    if w is not None: sh.width = Inches(w)
    if h is not None: sh.height = Inches(h)


def delete(sh):
    sh._element.getparent().remove(sh._element)


def box(slide, x, y, w, h, text="", size=13.5, bold=False, color=TEXT, fill=None,
        align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, kind=MSO_SHAPE.RECTANGLE, inset=0.08):
    shp = slide.shapes.add_shape(kind, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = RGBColor.from_string(fill)
    else:
        shp.fill.background()
    shp.line.fill.background()
    shp.shadow.inherit = False
    tf = shp.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(inset)
    tf.margin_top = tf.margin_bottom = Inches(0.04)
    p = tf.paragraphs[0]
    p.alignment = align
    for k, part in enumerate(text.split("|") if text else []):   # "|" = satır kırılımı
        if k:
            p.add_line_break()
        r = p.add_run()
        r.text = part
        f = r.font
        f.size, f.bold, f.name = Pt(size), bold, FONT
        f.color.rgb = RGBColor.from_string(color)
    return shp


def add_line(shp, text, size, bold, color=TEXT):
    """Kutuya ikinci bir paragraf ekler (ΔV satırında sayı + niteleyici)."""
    p = shp.text_frame.add_paragraph()
    p.alignment = shp.text_frame.paragraphs[0].alignment
    r = p.add_run()
    r.text = text
    r.font.size, r.font.bold, r.font.name = Pt(size), bold, FONT
    r.font.color.rgb = RGBColor.from_string(color)


def replace_cell_text(table, row, col, text):
    cell = table.cell(row, col)
    p = cell.text_frame.paragraphs[0]
    p.runs[0].text = text
    for r in p.runs[1:]:
        r._r.getparent().remove(r._r)


# ---------- 3 · Concept of Operations: yeniden tasarım ----------
s = slides[2]
S = shapes(s)
keep = {97, 99, 525}                              # yazar, başlık, kaynak satırı
for sid, sh in S.items():
    if sid not in keep:
        delete(sh)

X0, LW = 0.59, 1.35
PW = (12.15 - LW) / 5
PHASES = ["Insertion", "LEOP", "Orbit acquisition", "Routine operations", "Disposal"]
SPANS = ["separation", "3 days (1A)", "2–4 months", "7–12 years", "end of life"]
ROWS = [
    ("Propulsion", BAND_1, TEXT, 0.85, [
        "None: launcher|upper stage",
        "No burns planned;|avoidance if required",
        "Correct injection errors,|reach reference orbit",
        "Ground track, collision|avoidance, inclination",
        "Lower orbit, passivate;|re-entry ≤ 25 yr",
    ]),
    ("Sentinel-1", None, NAVY, 0.95, [
        "1C: Vega-C, ~700 km|5 Dec 2024",
        "1A avoidance burn:|torque ~10× predicted",
        "1A: ~500 burns, 8.4 kg|Apr–Aug 2014",
        "1C + 1D: 6-day repeat|since June 2026",
        "1B lowered 2023–24:|693 → ~580 × 610 km",
    ]),
    ("Sentinel-2", BAND_2, TEAL, 1.10, [
        "2C: Vega, ~775 km|5 Sep 2024",
        "2C: solar array|deployment confirmed|via Troll",
        "2C: 10.6 km low|tandem with 2A|from end Oct 2024",
        "2A, 2015–24:|64 trims,|24 avoidance pairs",
        "A/B design:|786 → 653 km",
    ]),
]
DV_ROW = [("—", ""), ("—", ""), ("7.4 m/s", "flown"), ("≈ 4 m/s", "projected, 12 yr"), ("73 m/s", "A/B design")]

y = 1.30
for k, name in enumerate(PHASES):                # evre şeridi; son chevron ızgaranın sağ kenarında biter
    kind = MSO_SHAPE.PENTAGON if k == 0 else MSO_SHAPE.CHEVRON
    w = PW + 0.10 if k < 4 else PW
    chev = box(s, X0 + LW + k * PW, y, w, 0.55, name, 14.25, True, "FFFFFF", NAVY,
               PP_ALIGN.CENTER, MSO_ANCHOR.MIDDLE, kind, inset=0.22)
    chev.adjustments[0] = 0.28
    box(s, X0 + LW + k * PW, y + 0.58, PW, 0.30, SPANS[k], 12, False, "082D55",
        align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
y = 2.40
for label, band, label_color, h, cells in ROWS:
    if band:
        box(s, X0, y, 12.15, h, fill=band)
    box(s, X0, y, LW, h, label, 13.5, True, label_color, anchor=MSO_ANCHOR.TOP)
    for k, text in enumerate(cells):
        box(s, X0 + LW + k * PW, y, PW, h, text, 13.5)
    y += h + 0.08
box(s, X0, y, 12.15, 0.02, fill=NAVY)             # ΔV satırının üst çizgisi
box(s, X0, y + 0.06, LW, 0.66, "ΔV, 2C", 14.25, True, NAVY, anchor=MSO_ANCHOR.TOP)
for k, (num, qual) in enumerate(DV_ROW):
    cell = box(s, X0 + LW + k * PW, y + 0.06, PW, 0.66, num, 14.25, True, NAVY, anchor=MSO_ANCHOR.TOP)
    if qual:
        add_line(cell, qual, 12, False, NAVY)
y += 0.12
assert y + 0.66 < 6.42, y
src = shapes(s)[525]
place(src, X0, 6.42, 12.15, 0.42)
set_text(src, "Sources   Arianespace VV24/VV25; ESA Operations; ISSFD 2015 papers 134, 137; "
              "SpaceOps 2025 papers 226, 617; Airbus GNC 2014")

# ---------- 7 · Sentinel-1 Manoeuvres ----------
s = slides[6]
set_text(find_text(s, "In-plane thrust"), "Simplified estimate, in-plane −21 to −30 %: 7.6–8.3 kg")

# ---------- 9 · Sentinel-2C ΔV Budget ----------
s = slides[8]
set_text(find_text(s, "Known items over"), "Illustrative 12-year budget")
table = [sh for sh in s.shapes if sh.has_table][0].table
BASIS = {
    1: "Flown · SpaceOps 617 pp. 15–16, Table 3",
    2: "Planned (0.57) · execution not reported",
    3: "Projected · 2A: 64 burns × 1.56 cm/s",
    4: "Projected · 2A: 24 pairs × ~4 cm/s",
    5: "Not published · 6 Nov burn (0.93) is in row 1",
    6: "A/B design value · Airbus 2014",
    7: "Mixes flown, projected and design values",
}
for r, text in BASIS.items():
    replace_cell_text(table, r, 2, text)
replace_cell_text(table, 5, 0, "Additional inclination control, 12 yr")
for col, w in zip(table.columns, (3.80, 1.20, 5.40, 1.75)):   # ilk sütun tek satıra sığsın
    col.width = Inches(w)
replace_cell_text(table, 7, 0, "Illustrative total")
foot = find_text(s, "SpaceOps 2025 paper 617 pp. 5")
set_text(foot, "SpaceOps 2025 paper 617 pp. 5, 15–17; Airbus GNC 2014 p. 24. [H] m₀ 1,143 kg, Isp 220 s, "
               "ideal rocket equation; 2A rates (Jul 2015–Mar 2024) scaled to 12 years. 10 of 2A's avoidance burns also did "
               "maintenance, so rows 3 and 4 may overlap.")
place(foot, y=6.50, h=0.48)

# ---------- 10 · Propellant Budget ----------
s = slides[9]
set_text(find_text(s, "Disposal plan"), "A/B design (ideal 37.8 kg)")
set_text(find_text(s, "6 Nov burn"), "6 Nov: 0.936 m/s in 314 s → ≈ 3.4 N average effective force")
set_text(find_text(s, "Loss vs. impulse"), "Finite-burn penalty (node-centred, constant thrust): ≈ 0.45 %")
set_text(find_text(s, "Allocation estimate"),
         "Scenario, not telemetry. 47 kg is the A/B design value; unallocated excludes inclination "
         "control and residuals.")
old = [sh for sh in s.shapes if sh.shape_type == 13][0]
assert old.name == "S2C hydrazine allocation", old.name
pic = s.shapes.add_picture(str(CHART), old.left, old.top, old.width, old.height)
old._element.addprevious(pic._element)
pic.name = "S2C hydrazine allocation (scenario)"
delete(old)

# ---------- 11 · Design Drivers ----------
s = slides[10]
set_text(find_text(s, "73 m/s and 47 kg"), "73 m/s, 47 kg (A/B design): 86 % of the illustrative 12-year Δv.")
set_text(find_text(s, "0.46 mm/s²"), "0.46 mm/s² on 1C; ≈ 3.4 N average effective force on 2C's largest burn.")
CLAIMS = {                                         # slogan yerine düz teknik ifade
    "Disposal dominates": "Disposal is the largest item",
    "Plumes limit burn length": "Plume impingement limits burn length",
    "1 N is enough": "1 N thrust is sufficient",
    "Margin buys extensions": "Propellant margin allows extensions",
}
for old_claim, new_claim in CLAIMS.items():
    set_text(find_text(s, old_claim), new_claim)

prs.save(OUT)
print("wrote", OUT)
