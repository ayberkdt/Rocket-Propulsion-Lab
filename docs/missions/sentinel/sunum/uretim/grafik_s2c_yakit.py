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
ORAN = {"v22": None, "v23": 0.85, "v25": 0.85, "v27": 0.85, "v28": 0.85, "v29": 0.85}[SURUM]
DEC = {"v22": 0, "v23": 0, "v25": 1, "v27": 1, "v28": 1, "v29": 1}[SURUM]          # v25: slayttaki 3.9 kg ile aynı yazım
OUT = Path(__file__).parent / "grafik"
OUT.mkdir(exist_ok=True)

G0 = 9.80665
LOADED = 133.0                          # ESA S2C fuelling blog (Aug 2024)
DV_FLOWN = 6.33 + 0.012 + 1.056442      # SpaceOps 2025 #617 pp.15-16
FLOWN = 1143.0 * (1 - math.exp(-DV_FLOWN / (220.0 * G0)))   # m0 VV24 p.3, Isp 220 s assumed
DISPOSAL = 47.0                         # Airbus GNC 2014, A/B disposal plan (73 m/s)

if SURUM in ("v27", "v28", "v29"):
    # v28: son çubuk "Unallocated (scenario)" (ikinci geri bildirim); başka değişiklik yok.
    # v29: planlanan 0.57 m/s (uçtuğu bildirilmedi) bütçeden çıktı; etiketler tam yazıldı.
    SON = "Rest" if SURUM == "v27" else "Unallocated\n(scenario)"
    # v27: "kalan yakıt" değil, tahsis tahmini. Planlanan kalemler açık ton + mavi kenar,
    # son çubuk taranmış: ölçülmüş rezerv olmadığını görsel olarak da söyler.
    # 12 yıllık ufuk [H]: 2A oranları SpaceOps 2025 #617 s.5 (2 Tem 2015 - Mar 2024, 8.71 yıl):
    # 64 OCM x 1.56 cm/s ; 24 CAM çifti x 2 x ~4 cm/s (makalenin örnek CAM boyu).
    YIL = 2024.21 - 2015.50
    PLAN_DV = [0.57, 0.998 / YIL * 12, 24 * 2 * 0.04 / YIL * 12]
    if SURUM == "v29":
        PLAN_DV = PLAN_DV[1:]
    m = 1143.0 - FLOWN
    PLAN = 0.0
    for dv in PLAN_DV:
        p = m * (1 - math.exp(-dv / (220.0 * G0)))
        m -= p
        PLAN += p
    LEFT = LOADED - FLOWN - PLAN - DISPOSAL
    if SURUM == "v29":
        assert abs(FLOWN - 3.91) < 0.01 and abs(PLAN - 2.12) < 0.01 and abs(LEFT - 79.97) < 0.02, (FLOWN, PLAN, LEFT)
    else:
        assert abs(FLOWN - 3.91) < 0.01 and abs(PLAN - 2.42) < 0.01 and abs(LEFT - 79.67) < 0.02, (FLOWN, PLAN, LEFT)
    if SURUM in ("v28", "v29"):
        # Yatay şelale: dikeyde "Unallocated (scenario)" x eksenine sığmadı (tik-cakisma HATA);
        # yatayda kategori adları tam yazılır. Renk/dolgu kodu v27 ile aynı.
        mpl.rcParams["hatch.color"] = "#5E7389"
        fig, ax = dw.figur("sunum", hedef="sunum-yarim", palet="okyanus", oran=0.85)
        mpl.rcParams["savefig.facecolor"] = "white"
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        adlar = ["Loaded", "Flown 2024", "Ops, 12 yr", "EOL (A/B design)", "Unallocated\n(scenario)"]
        if SURUM == "v29":
            adlar = ["Loaded", "Flown, 2024", "Operations,\n12 years", "Disposal\n(A/B design)",
                     "Unallocated\n(scenario)"]
        sol = [0, LOADED - FLOWN, LOADED - FLOWN - PLAN, LEFT, 0]
        boy = [LOADED, FLOWN, PLAN, DISPOSAL, LEFT]
        stil = [("#8299AF", "#8299AF", None), ("#2967B4", "#2967B4", None), ("#C9D8EE", "#2967B4", None),
                ("#C9D8EE", "#2967B4", None), ("white", "#5E7389", "///")]
        yk = list(range(len(adlar)))
        for k in yk:
            fc, ec, h = stil[k]
            ax.barh(k, boy[k], left=sol[k], height=0.62, color=fc, edgecolor=ec, linewidth=1.2, hatch=h)
        for k in yk[:-1]:                       # bağlayıcı: bir çubuğun bittiği seviye bir sonrakinin başı
            x = LOADED if k == 0 else sol[k]
            ax.plot([x, x], [k + 0.31, k + 0.69], color="#9AA7B4", lw=0.8, zorder=1)
        for k in yk:
            if k in (0, len(yk) - 1):
                ax.annotate(f"{sol[k] + boy[k]:.1f}", (sol[k] + boy[k], k), xytext=(11, 0),
                            textcoords="offset points", va="center", ha="left")
            else:
                ax.annotate(f"−{boy[k]:.1f}", (sol[k], k), xytext=(-11, 0),
                            textcoords="offset points", va="center", ha="right")
        ax.set_yticks(yk, adlar)
        ax.invert_yaxis()
        ax.set_xlim(0, 165)
        ax.set_xticks([0, 40, 80, 120, 160])
        ax.set_xlabel("Hydrazine (kg)")
        dw.izgara.uygula(ax, "silik", "x")      # değer ekseni x: ızgara yalnız x'te
        dw.kaydet(fig, OUT / f"s2c_yakit_{SURUM}", hedef="sunum-yarim", cvd=True)
        print(f"flown {FLOWN:.3f} plan {PLAN:.3f} disposal {DISPOSAL} unallocated {LEFT:.3f}")
        sys.exit(0)
    fig, ax = dw.figur("sunum", hedef="sunum-yarim", palet="okyanus", oran=0.85)
    mpl.rcParams["savefig.facecolor"] = "white"
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    dw.teknik.selale(ax, ["Loaded", "Flown\n2024", "Ops\n12 yr", "EOL\nplan"],
                     [LOADED, -FLOWN, -PLAN, -DISPOSAL], baslangic="Loaded", bitis=SON,
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
