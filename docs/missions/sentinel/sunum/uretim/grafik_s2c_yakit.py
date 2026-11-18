"""S2C hidrazin şelalesi, V21 slayt 7 grafiğinin renkleriyle (mavi 2967B4, gri 8299AF, beyaz zemin).

Kullanım: python grafik_s2c_yakit.py v25   ->  grafik/s2c_yakit_v25.*
v22 varsayılan oranla üretildi; v23 slayttaki karta göre daha uzun (oran 0.85).
"""
import math
import sys
from pathlib import Path

import matplotlib as mpl
import datawiz as dw

SURUM = sys.argv[1] if len(sys.argv) > 1 else "v23"
ORAN = {"v22": None, "v23": 0.85, "v25": 0.85, "v27": 0.85}[SURUM]
DEC = {"v22": 0, "v23": 0, "v25": 1, "v27": 1}[SURUM]          # v25: slayttaki 3.9 kg ile aynı yazım
OUT = Path(__file__).parent / "grafik"
OUT.mkdir(exist_ok=True)

G0 = 9.80665
LOADED = 133.0                          # ESA S2C fuelling blog (Aug 2024)
DV_FLOWN = 6.33 + 0.012 + 1.056442      # SpaceOps 2025 #617 pp.15-16
FLOWN = 1143.0 * (1 - math.exp(-DV_FLOWN / (220.0 * G0)))   # m0 VV24 p.3, Isp 220 s assumed
DISPOSAL = 47.0                         # Airbus GNC 2014, A/B disposal plan (73 m/s)

if SURUM == "v27":
    # v27: "kalan yakıt" değil, tahsis tahmini. Planlanan kalemler açık ton + mavi kenar,
    # son çubuk taranmış: ölçülmüş rezerv olmadığını görsel olarak da söyler.
    # 12 yıllık ufuk [H]: 2A oranları SpaceOps 2025 #617 s.5 (2 Tem 2015 - Mar 2024, 8.71 yıl):
    # 64 OCM x 1.56 cm/s ; 24 CAM çifti x 2 x ~4 cm/s (makalenin örnek CAM boyu).
    YIL = 2024.21 - 2015.50
    PLAN_DV = [0.57, 0.998 / YIL * 12, 24 * 2 * 0.04 / YIL * 12]
    m = 1143.0 - FLOWN
    PLAN = 0.0
    for dv in PLAN_DV:
        p = m * (1 - math.exp(-dv / (220.0 * G0)))
        m -= p
        PLAN += p
    LEFT = LOADED - FLOWN - PLAN - DISPOSAL
    assert abs(FLOWN - 3.91) < 0.01 and abs(PLAN - 2.42) < 0.01 and abs(LEFT - 79.67) < 0.02, (FLOWN, PLAN, LEFT)
    fig, ax = dw.figur("sunum", hedef="sunum-yarim", palet="okyanus", oran=0.85)
    mpl.rcParams["savefig.facecolor"] = "white"
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    dw.teknik.selale(ax, ["Loaded", "Flown\n2024", "Ops\n12 yr", "EOL\nplan"],
                     [LOADED, -FLOWN, -PLAN, -DISPOSAL], baslangic="Loaded", bitis="Rest",
                     renkler=["#2967B4", "#2967B4"], bicim="{:+.1f}")
    bars = sorted(ax.patches, key=lambda b: b.get_x())
    stil = [("#8299AF", "#8299AF", None),          # Loaded: ölçülmüş dolum
            ("#2967B4", "#2967B4", None),          # Flown 2024: uçmuş manevraların ideal hesabı
            ("#C9D8EE", "#2967B4", None),          # Planned: açık ton, mavi kenar
            ("#C9D8EE", "#2967B4", None),          # Disposal plan
            ("white", "#5E7389", "///")]           # Not allocated: rezerv değil
    for bar, (fc, ec, h) in zip(bars, stil):
        bar.set_facecolor(fc)
        bar.set_edgecolor(ec)
        bar.set_linewidth(1.2)
        if h:
            bar.set_hatch(h)
    mpl.rcParams["hatch.color"] = "#5E7389"
    for text in ax.texts:
        raw = text.get_text().replace("−", "-")
        try:
            value = float(raw)
        except ValueError:
            continue
        text.set_text(f"{value:+.1f}".replace("-", "−") if raw.startswith(("+", "-")) else f"{value:.1f}")
        text.xyann = (0, 8)               # bağlayıcı çizgiden bir harf boşluğu (cizgi-nefes)
    ax.tick_params(axis="x", labelrotation=0)
    ax.set_ylim(0, 175)
    ax.set_yticks([0, 40, 80, 120, 160])     # 150 ızgarası 133 / 129 etiketlerine değiyordu (cizgi-nefes)
    ax.set_ylabel("Hydrazine (kg)")
    dw.kaydet(fig, OUT / f"s2c_yakit_{SURUM}", hedef="sunum-yarim", cvd=True)
    print(f"flown {FLOWN:.3f} plan {PLAN:.3f} disposal {DISPOSAL} not allocated {LEFT:.3f}")
    sys.exit(0)

fig, ax = dw.figur("sunum", hedef="sunum-yarim", palet="okyanus", oran=ORAN)
mpl.rcParams["savefig.facecolor"] = "white"
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
dw.teknik.selale(ax, ["Loaded", "Flown", "Disposal"], [LOADED, -FLOWN, -DISPOSAL],
                 baslangic="Loaded", bitis="Left", renkler=["#2967B4", "#2967B4"], bicim="{:+.1f}")
bars = sorted(ax.patches, key=lambda b: b.get_x())   # V21 slide 7 colours: totals grey, changes blue
for k, bar in enumerate(bars):
    bar.set_facecolor("#8299AF" if k in (0, len(bars) - 1) else "#2967B4")
for text in ax.texts:
    raw = text.get_text().replace("\u2212", "-")
    try:
        value = float(raw)
    except ValueError:
        continue
    text.set_text(f"{value:+.{DEC}f}".replace("-", "\u2212") if raw.startswith(("+", "-")) else f"{value:.{DEC}f}")
ax.tick_params(axis="x", labelrotation=0)
ax.set_ylim(0, 170)
ax.set_yticks([0, 50, 100, 150])
ax.set_ylabel("Hydrazine (kg)")
dw.kaydet(fig, OUT / f"s2c_yakit_{SURUM}", hedef="sunum-yarim", cvd=True)
