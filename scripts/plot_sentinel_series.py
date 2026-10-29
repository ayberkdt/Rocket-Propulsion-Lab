"""Author the Sentinel report's scientific figures; matplotlib is an authoring dependency."""
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/missions/sentinel"
FIGURES = DOCS / "figures"


def functional(series, count, notes):
    fig, ax = plt.subplots(figsize=(12, 5.3))
    fig.patch.set_facecolor("#f4f8fb")
    ax.set(xlim=(0, 12), ylim=(0, 6))
    ax.axis("off")
    boxes = [(0.3, 3, "Basınçlı gaz / hidrazin\nBlow-down tank"),
             (3.3, 3, "İzolasyon / hat\nFiltre ve valf işlevleri"),
             (6.3, 3, "Ana / yedek dal\n" + count),
             (9.3, 3, "Katalitik iticiler\nKuvvet + moment"),
             (3.3, .7, "Basınç / sıcaklık\nAOCS / teker verileri"),
             (6.3, .7, "Kontrol / FDIR\nValf komutu / durdurma")]
    for x, y, label in boxes:
        ax.add_patch(FancyBboxPatch((x, y), 2.4, 1.25, boxstyle="round,pad=0.08",
                                   facecolor="white", edgecolor="#246681", linewidth=1.5))
        ax.text(x + 1.2, y + .625, label, ha="center", va="center", fontsize=10)
    for start, end in [((2.8, 3.625), (3.15, 3.625)), ((5.8, 3.625), (6.15, 3.625)),
                       ((8.8, 3.625), (9.15, 3.625)), ((5.8, 1.325), (6.15, 1.325)),
                       ((7.5, 2.05), (7.5, 2.85))]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                    color="#246681"))
    ax.text(6, 5.45, series + " — işlevsel itki/besleme şeması", ha="center", fontsize=16,
            weight="bold")
    ax.text(6, 4.85, "Kendi çizimimiz · uçuş P&ID'si veya kesin valf/hat sayısı değildir",
            ha="center", fontsize=10, color="#4d6470")
    ax.text(6, .1, notes, ha="center", fontsize=10, color="#4d6470")
    fig.savefig(FIGURES / f"S{series[-1]}_functional.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    FIGURES.mkdir(exist_ok=True)
    functional("Sentinel-1", "2 × 7 adet 1 N\nYönler/işlevler farklı",
               "Aile/B mimarisi · kaynaklar: SpaceOps 2025 #360 s.2; Rafael kataloğu")
    functional("Sentinel-2", "2 × 4 eğik 1 N\nSoğuk yedekleme",
               "A/B mimari mirası · C'ye özgü kabul haritası açık değil · kaynak: Airbus AOCS 2014")
    spec = importlib.util.spec_from_file_location("sentinel_calc", ROOT / "scripts/sentinel_series_budget.py")
    assert spec is not None and spec.loader is not None
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    data = json.loads((DOCS / "SERI_VERILERI.json").read_text(encoding="utf-8"))
    result = model.calculate(data)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), constrained_layout=True)
    ax = axes[0]
    labels = ["1A", "1B", "1C", "1D"]
    bars = ax.bar(labels, [2157, 2164, 2186, 2184], color=["#9db4bf", "#9db4bf", "#176b91", "#2e8a80"])
    ax.set(ylim=(0, 2500), ylabel="Yayımlanan araç kütlesi (kg)", title="Sentinel-1: ayrı kit sürümleri")
    ax.bar_label(bars, padding=3)
    ax.text(.5, -.18, "A/B/C: birincil yazarlı kit aynası; D: resmî kit\nYakıtlar bütün araçlara eşit atanmadı",
            transform=ax.transAxes, ha="center", fontsize=9)
    ax = axes[1]
    cases = result["s2c_fuel_sensitivity"]
    bars = ax.bar([f"{r['isp_s_assumed']} s\nη={r['axial_factor_assumed']}" for r in cases],
                  [r["fuel_used_ideal_kg"] for r in cases], color="#176b91")
    ax.set(ylabel="İdeal model yakıtı (kg)", ylim=(0, 7), title="Sentinel-2C: ≈7,97 m/s için duyarlılık")
    ax.bar_label(bars, fmt="%.2f", padding=3)
    ax.text(.5, -.18, "Isp/η varsayımları; ölçülmüş tüketim değil\nKaynak: SpaceOps #617; roket denklemi",
            transform=ax.transAxes, ha="center", fontsize=9)
    fig.savefig(FIGURES / "series_comparison.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
