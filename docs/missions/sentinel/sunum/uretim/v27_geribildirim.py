"""v27: dış geri bildirimin doğrulanan maddeleri, iki yeni slayt (11 slayt), dolgu cümlelerin silinmesi.

Doğrulama (2026-11-18, birincil kaynaktan okundu):
  CDSE 09.10.2026      : 2A uzatması 31 Aralık 2027 (v25'teki end-2026 yanlıştı).
  SpaceOps 2025 #617   : s.16 "a small OOP OCM was planned ... (total delta-v of 0.57 m/s)" -> uçtuğu/uçmadığı
                         yazmıyor; Tablo 3-4 (6 Kasım S2C 0.936, S2A 0.92 m/s); s.5 2A: 2 Tem 2015 - Mar 2024,
                         64 OCM, ort. 1.56 cm/s, 24 CAM çifti, örnek CAM 4 cm/s.
  ISSFD 2015 #134      : s.8 Tablo 2 enjeksiyon düzeltmesi 4.1 + 2.2 m/s; itici 1: -%21 performans, ~10x tork,
                         plume; IP iticilerde ~-%30; IP yanma <= 30 s (~1 cm/s); s.14 ~500 manevra, 8.4 kg.
  ISSFD 2015 #137      : s.11-14 40 s LEOP CAM'i, tork bir mertebe büyük, IP yanma 300 s -> 30 s.
  SENTINEL-2_SERISI.md : S2A fırlatma 23.06.2015 (slayt 6'daki "since 2014" S2 için yanlıştı).
  SENTINEL-1_SERISI.md : S1A kit kütlesi 2157 kg; 2300 kg eski tasarım bütçesi, A için kullanılmaz.

Slayt düzeni: 1-8 düzeltilir, 9 yeni (S2C ΔV bütçesi), 10 eski 9 (yakıt tahsisi), 11 yeni (tasarım etkenleri).
Girdi: surumler/..._v26.pptx  ·  Çıktı: surumler/..._v27.pptx  ·  Grafik: grafik_s2c_yakit.py v27
Hesaplar [H]: g0 = 9.80665 m/s², ideal roket denklemi, sıralı kütle güncellemesi.
"""
import copy
import math
from pathlib import Path

from pptx import Presentation
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from pptx.dml.color import RGBColor
from pptx.util import Inches

HERE = Path(__file__).resolve().parent
SUNUM = HERE.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v26.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v27.pptx"
CHART = HERE / "grafik" / "s2c_yakit_v27.png"

# ---------------------------------------------------------------- hesaplar
G0 = 9.80665


def mp(m0, dv, isp):
    return m0 * (1 - math.exp(-dv / (isp * G0)))


# S1A: ISSFD 2015 #134 Tablo 2; m0 VS07 kiti 2157 kg; Isp 214 s (Rafael, 22 bar kararlı) [H]
S1A_M0, S1A_ISP = 2157.0, 214.0
S1A_IP = mp(S1A_M0, 4.1, S1A_ISP)
S1A_IDEAL = S1A_IP + mp(S1A_M0 - S1A_IP, 2.2, S1A_ISP)


def s1a_degraded(eff):
    ip = mp(S1A_M0, 4.1, S1A_ISP * eff)            # IP itkisi plume ile düşer; OOP iticisi etkilenmedi (#134 s.8)
    return ip + mp(S1A_M0 - ip, 2.2, S1A_ISP)


S1A_LO, S1A_HI = s1a_degraded(0.79), s1a_degraded(0.70)
assert abs(S1A_IDEAL - 6.47) < 0.01 and abs(S1A_LO - 7.58) < 0.01 and abs(S1A_HI - 8.27) < 0.01

# S2C ΔV bütçesi, 220 s, 1143 kg (VV24)
ISP = 220.0
YIL = 2024.21 - 2015.50                              # 2 Tem 2015 - Mar 2024
DV = {
    "acq": 6.33 + 0.012 + 0.0146 + math.hypot(0.9344, 0.0537) + 0.0025 + 0.0207 + 0.0206 + 0.0309 + 0.0029 + 0.0283,
    "oop": 0.57,
    "maint": 64 * 0.0156 / YIL * 12,
    "cam": 24 * 2 * 0.04 / YIL * 12,
    "eol": 73.0,
}
m = 1143.0
KG = {}
for k in ("acq", "oop", "maint", "cam"):
    KG[k] = mp(m, DV[k], ISP)
    m -= KG[k]
EOL_IDEAL = mp(m, DV["eol"], ISP)
DV_TOTAL = sum(DV.values())
KG_TOTAL = sum(KG.values()) + 47.0
SHARE = DV["eol"] / DV_TOTAL
assert abs(DV["acq"] - 7.398) < 0.001 and abs(KG["acq"] - 3.91) < 0.01
assert abs(DV["maint"] - 1.375) < 0.01 and abs(DV["cam"] - 2.645) < 0.01
assert abs(DV_TOTAL - 84.99) < 0.02 and abs(SHARE - 0.859) < 0.001 and abs(EOL_IDEAL - 37.8) < 0.05
assert abs(KG_TOTAL - 53.33) < 0.02

MDOT = 1.0 / (ISP * G0)                              # kg/s, 1 N
BURN_H = KG["acq"] / MDOT / 3600
PERIOD = 2 * math.pi * math.sqrt((6378.137e3 + 786e3) ** 3 / 3.986004418e14)
FRAC = 314.1 / PERIOD
ARC = 2 * math.pi * FRAC
LOSS = 1 - math.sin(ARC / 2) / (ARC / 2)             # düğüm merkezli düzlem dışı yanmanın etkinlik kaybı
F_NET = (1143.0 - KG["acq"]) * math.hypot(0.9344, 0.0537) / 314.1
assert abs(MDOT * 1000 - 0.4635) < 0.0005 and abs(BURN_H - 2.34) < 0.01
assert abs(PERIOD / 60 - 100.58) < 0.01 and abs(FRAC - 0.052) < 0.001
assert abs(math.degrees(ARC) - 18.7) < 0.1 and abs(LOSS - 0.00445) < 0.0001 and abs(F_NET - 3.39) < 0.01
A_S1C, A_S2C = 1 / 2186, 1 / 1143

# ---------------------------------------------------------------- yardımcılar
prs = Presentation(SRC)
slides = list(prs.slides)
assert len(slides) == 9


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


def replace_in_table(table, old, new):
    n = 0
    for row in table.rows:
        for cell in row.cells:
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    if old in r.text:
                        r.text = r.text.replace(old, new)
                        n += 1
    assert n == 1, (old, n)


def set_cell(cell, text, align=None):
    p = cell.text_frame.paragraphs[0]
    if not p.runs:
        p.add_run()
    set_text(cell, text)
    if align is not None:
        p.alignment = align


def slide_from(src_slide, keep_ids):
    """src_slide'ın düzeniyle yeni slayt; yalnız keep_ids şekilleri (resim değil) kopyalanır."""
    s = prs.slides.add_slide(src_slide.slide_layout)
    for ph in list(s.placeholders):
        delete(ph)
    src = shapes(src_slide)
    for sid in keep_ids:
        s.shapes._spTree.append(copy.deepcopy(src[sid]._element))
    return s


# ---------------------------------------------------------------- yeni slaytlar (kaynaklar değişmeden önce)
s8, s9 = slides[7], slides[8]
budget = slide_from(s8, [139, 140, 141, 142, 950, 151])
drivers = slide_from(s9, [145, 146, 1013, 1018, 1021])

# ---------- 1 · Mission Details ----------
s = slides[0]
delete(find_text(s, "Spacecraft Architectures and Subsystems Design"))
replace_in_table([sh for sh in s.shapes if sh.has_table][0].table, "end-2026", "end-2027")
set_text(find_text(s, "Fleet status on"),
         "Fleet status on 12 Nov 2026. 1A ended 29 Jun 2026 (ESA); 2A extended to 31 Dec 2027 "
         "(Copernicus, 9 Oct 2026); 2D in storage (ESA FutureEO 2026).")

# ---------- 2 · Spacecraft and Launches ----------
s = slides[1]
S = shapes(s)
delete(S[974])                                       # dolgu cümle; yakıt miktarı slayt 5'te

# ---------- 3 · CONOPS: dolgu satırları ve fırlatma fotoğrafı çıkar, iki teknik görsel büyür ----------
s = slides[2]
S = shapes(s)
for sid in (1000, 1001, 2, 500):
    delete(S[sid])
H = 3.0
w_dep, w_plot = H * 3.60 / 2.65, H * 5.20 / 2.39
x0 = 0.59 + (12.15 - (w_dep + 0.35 + w_plot)) / 2
place(S[1002], x0, 1.25, w_dep, H)
place(S[1003], x0 + w_dep + 0.35, 1.25, w_plot, H)
place(S[501], x0, 4.30, w_dep, 0.25)
place(S[502], x0 + w_dep + 0.35, 4.30, w_plot, 0.25)
for sid in range(503, 524):
    if sid in S:
        S[sid].top = S[sid].top + Inches(0.27)

# ---------- 4 · Spacecraft Systems ----------
s = slides[3]
S = shapes(s)
replace_in_table(S[90].table, "<0.004° per axis", "< 70 µrad (0.004°) per axis")
set_text(S[1020], "One 1 N thruster")
delete(S[1021])
place(S[1022], y=5.10)
set_text(S[1022], f"a = 1 N / 2,186 kg ≈ {A_S1C * 1e4:.2f} × 10⁻⁴ m/s²")
place(S[1023], y=5.68)
set_text(S[1023], f"30 s firing (1A in-plane limit) → Δv ≈ {30 * A_S1C * 100:.1f} cm/s.  "
                  f"Sentinel-2C, 1,143 kg: {A_S2C * 1e4:.1f} × 10⁻⁴ m/s².")

# ---------- 5 · Hydrazine Propulsion ----------
s = slides[4]
set_text(find_text(s, "2C tank charging pressure"),
         "2C pressurised with helium to 21 bar (ESA; abs/gauge not stated). "
         "Thruster layout: Sentinel-1B, SpaceOps 2025 paper 360, Fig. 1.")

# ---------- 6 · Catalogue / TRL ----------
s = slides[5]
replace_in_table([sh for sh in s.shapes if sh.has_table][0].table,
                 "5.5–22 bar. Pressure reference unspecified.", "5.5–22 bar (abs/gauge not stated)")
trl = find_text(s, "TRL 9")
set_text(trl, "TRL 9 for the families: Rafael 1 N on Sentinel-1 since 2014, CHT1N on Sentinel-2 since 2015. "
              "C/D variants unpublished. Specs: Rafael p. 46, ArianeGroup brochure p. 6. Photo: ArianeGroup.")
place(trl, y=6.72, h=0.34)

# ---------- 7 · Sentinel-1 Manoeuvres ----------
s = slides[6]
S = shapes(s)
place(S[1010], y=3.40)
set_text(S[1010], f"Injection errors: 4.1 + 2.2 = 6.3 m/s → {S1A_IDEAL:.1f} kg ideal [H]")
clone(s, S[1010], 0.59, 3.82, 6.93, 0.42,
      f"In-plane thrust −21 to −30 % (plume) → {S1A_LO:.1f}–{S1A_HI:.1f} kg; flown 8.4 kg")
place(S[1005], y=4.40)
place(S[1006], y=4.85, h=0.70)
set_text(S[1006], "1A: torques ~10× predicted on a 40 s LEOP burn; burns cut from 300 s to 30 s (≈ 1 cm/s each).")
place(S[1007], y=5.70)
set_text(S[113], "ISSFD 2015 #134 pp. 8, 14 and #137 pp. 11–14; SpaceOps 2025 #226, #360. "
                 f"[H] m₀ {S1A_M0:,.0f} kg (VS07 kit), Isp {S1A_ISP:.0f} s (Rafael, 22 bar).")

# ---------- 8 · Sentinel-2C Δv ----------
s = s8
S = shapes(s)
set_text(S[144], "6 Nov 2024 inclination burn")
set_text(S[149], "Same day as Sentinel-2A's 0.92 m/s burn. Counted in Table 3.")
set_text(S[970], "Plane corrections after LEOP: 0.57 m/s planned, execution not reported → up to ≈7.97 m/s.")
set_text(S[950], "23 Nov avoidance pair: 0.0413 m/s, inside Table 3; net along-track ≈ 0.")
set_text(S[151], "SpaceOps 2025, paper 617, pp. 15–17, Tables 3–4")

# ---------- 9 (yeni) · Sentinel-2C ΔV Budget ----------
s = budget
S = shapes(s)
set_text(S[139], "Sentinel-2C ΔV Budget")
set_text(S[141], "Known items over the 12-year propellant horizon")
tbl = S[142].table._tbl
grid = tbl.find(qn("a:tblGrid"))
last_col = grid.findall(qn("a:gridCol"))[-1]
new_col = copy.deepcopy(last_col)
for col_id in new_col.iter("{*}colId"):                # a16:colId tekil olmalı
    col_id.set("val", str(int(col_id.get("val")) + 1000))
last_col.addnext(new_col)
for tr in tbl.findall(qn("a:tr")):
    last_tc = tr.findall(qn("a:tc"))[-1]
    last_tc.addnext(copy.deepcopy(last_tc))            # extLst'ten önce: tc'ler bitişik kalmalı
trs = tbl.findall(qn("a:tr"))
row_a, row_b, row_total = trs[1], trs[2], trs[4]
for tr in trs[1:]:
    tbl.remove(tr)
ROWS = [
    ("Acquisition and tandem, 2024", f"{DV['acq']:.1f}", "Flown; SpaceOps 617 pp. 15–16, Table 3", f"{KG['acq']:.1f}"),
    ("Plane corrections after LEOP", f"{DV['oop']:.1f}", "Planned (0.57); execution not reported", f"{KG['oop']:.1f}"),
    ("In-plane maintenance, 12 yr", f"{DV['maint']:.1f}", "2A rate: 64 burns × 1.56 cm/s, 2015–2024", f"{KG['maint']:.1f}"),
    ("Collision avoidance, 12 yr", f"{DV['cam']:.1f}", "2A: 24 pairs × ~4 cm/s, 2015–2024 [H]", f"{KG['cam']:.1f}"),
    ("Inclination control", "—", "Not published; one burn was 0.93 m/s (6 Nov 2024)", "—"),
    ("Disposal, 786 → 653 km", f"{DV['eol']:.1f}", "A/B design study (Airbus 2014)", f"47.0 ({EOL_IDEAL:.1f} ideal)"),
]
# datawiz tablo.denetle: sütun başına tek ondalık (önce 1-2 karışık → şimdi 1)
for k, _ in enumerate(ROWS):
    tbl.append(copy.deepcopy(row_a if k % 2 == 0 else row_b))
tbl.append(copy.deepcopy(row_total))
table = S[142].table
WIDTHS = [3.30, 1.30, 5.80, 1.75]
assert abs(sum(WIDTHS) - 12.15) < 1e-9
for col, w in zip(table.columns, WIDTHS):
    col.width = Inches(w)
HEADER = ("Item", "Δv (m/s)", "Basis", "Propellant (kg)")
TOTAL = ("Known total", f"≈ {DV_TOTAL:.1f}", f"Disposal is {SHARE * 100:.0f} % of the known Δv", f"≈ {KG_TOTAL:.1f}")
for r, values in enumerate([HEADER, *ROWS, TOTAL]):
    for c, text in enumerate(values):
        set_cell(table.cell(r, c), text, PP_ALIGN.RIGHT if c in (1, 3) else PP_ALIGN.LEFT)
for row in table.rows:
    row.height = Inches(0.44)
place(S[142], 0.59, 1.97, 12.15, 0.44 * 8)
delete(S[950])                                       # toplam satırı mesajı zaten söylüyor (Occam)
set_text(S[151], "SpaceOps 2025 paper 617 pp. 5, 15–17; Airbus GNC 2014 p. 24. "
                 "[H] m₀ 1,143 kg, Isp 220 s, ideal rocket equation, 12-year horizon.")

# ---------- 10 (eski 9) · Propellant Budget ----------
s = s9
S = shapes(s)
set_text(S[1010], "Sentinel-2C hydrazine allocation (estimate)")
set_text(S[1015], f"Flown 2024 ({DV['acq']:.2f} m/s)")
set_text(S[1016], "47 kg")
set_text(S[1017], f"Disposal plan (ideal {EOL_IDEAL:.1f} kg)")
delete(S[1018])
place(S[1019], y=3.45)
place(S[1020], y=4.00)
place(S[1027], y=4.43)
set_runs(S[1027], [("ṁ = F / (I", None), ("sp", "sub"), (" g", None), ("0", "sub"),
                   (f") = {MDOT * 1000:.2f} g/s at 1 N; {KG['acq']:.1f} kg ≈ {BURN_H:.1f} h at 1 N", None)])
clone(s, S[1020], 0.59, 4.86, 6.93, 0.40,
      f"6 Nov burn: 0.936 m/s in 314 s → {F_NET:.1f} N net from the thrusters")
place(S[1028], y=5.29)
set_text(S[1028], f"Longest burn 314 s = {FRAC * 100:.1f} % of orbit → {math.degrees(ARC):.0f}° of arc")
clone(s, S[1020], 0.59, 5.72, 6.93, 0.40,
      f"Loss vs. impulse at the node: {LOSS * 100:.2f} % → impulsive model holds")
place(S[1021], y=6.15)
set_text(S[1024], "Allocation estimate, not telemetry. Rest excludes inclination control and unusable residuals.")
place(S[1024], y=5.92, h=0.62)
old = S[1026]
pic = s.shapes.add_picture(str(CHART), old.left, old.top, old.width, old.height)
old._element.addprevious(pic._element)
pic.name = "S2C hydrazine allocation"
delete(old)

# ---------- 11 (yeni) · Design Drivers ----------
s = drivers
S = shapes(s)
set_text(S[145], "Design Drivers")
claim_src, text_src = S[1018], S[1021]
DRIVERS = [
    ("Disposal dominates", f"73 m/s and 47 kg of 133 kg on Sentinel-2C: {SHARE * 100:.0f} % of the known Δv."),
    ("Routine Δv is small", f"2C acquisition {DV['acq']:.1f} m/s; 2A maintenance {DV['maint'] / 12:.2f} m/s per year "
                            "(2015–2024)."),
    ("Plumes limit burn length", "1A in-plane burns cut from 300 s to 30 s → nearly 500 manoeuvres in four months."),
    ("1 N is enough", f"{A_S1C * 1000:.2f} mm/s² on 1C; 2C's thrusters gave {F_NET:.1f} N net on its largest burn."),
    ("Margin buys extensions", "Propellant for 12 years vs 7-year design life: 1A flew 12.2 years; "
                               "2A runs to end-2027."),
]
for k, (claim, evidence) in enumerate(DRIVERS):
    y = 1.55 + k * 1.0
    c = clone(s, claim_src, 0.59, y, 4.0, 0.80, claim)
    clone(s, text_src, 4.75, y, 7.99, 0.80, evidence)
    if k == 0:                                       # tek vurgu: destedeki denklem rengi (teal 066764)
        c.text_frame.paragraphs[0].runs[0].font.color.rgb = RGBColor(0x06, 0x67, 0x64)
delete(claim_src)
delete(text_src)
set_text(S[1013], "Sources   as on slides 3–10: SpaceOps 2025 papers 226, 360, 617; ISSFD 2015 papers 134, 137; "
                  "Airbus GNC 2014; ESA; Copernicus")

# ---------------------------------------------------------------- sıra: bütçe slaytı 9. sıraya
sld = prs.slides._sldIdLst
ids = list(sld)
assert len(ids) == 11
sld.remove(ids[9])
sld.insert(8, ids[9])

prs.save(OUT)
print("wrote", OUT)
print(f"S1A ideal {S1A_IDEAL:.2f} kg, degraded {S1A_LO:.2f}-{S1A_HI:.2f} kg")
print(f"S2C ΔV {', '.join(f'{k} {v:.3f}' for k, v in DV.items())} total {DV_TOTAL:.2f} share {SHARE:.3f}")
print(f"S2C kg {', '.join(f'{k} {v:.3f}' for k, v in KG.items())} EOL ideal {EOL_IDEAL:.2f} total {KG_TOTAL:.2f}")
print(f"mdot {MDOT*1000:.4f} g/s, {BURN_H:.2f} h, T {PERIOD/60:.2f} min, arc {math.degrees(ARC):.1f} deg, "
      f"loss {LOSS*100:.3f} %, F_net {F_NET:.2f} N")
