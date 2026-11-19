"""hesap.py'ye bilinçli hatalar sokar, tests/test_sentinel_sunum.py'nin yakaladığını doğrular, dosyayı her seferinde geri yükler.

Kullanım (depo kökünden): python docs/missions/sentinel/sunum/uretim/mutasyon_denetimi.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]   # uretim → sunum → sentinel → missions → docs → depo kökü
HESAP = ROOT / "docs/missions/sentinel/sunum/uretim/hesap.py"
YEDEK = HESAP.with_suffix(".py.yedek")

MUTASYONLAR = [
    ("Isp varsayımı 220 → 222 s", '"isp_2c": (220.0,', '"isp_2c": (222.0,'),
    ("Sonlu yanma cezası yanlış formül", "return 1 - math.sin(x) / x if x else 0.0", "return 1 - math.cos(x) if x else 0.0"),
    ("2A bakım ortalaması 1,56 → 1,65 cm/s", '"occ_ort": (0.0156,', '"occ_ort": (0.0165,'),
    ("expm1 yerine saf 1 − exp", "return m0 * -math.expm1(-dv / (verim * isp * G0))",
     "return m0 * (1 - math.exp(-dv / (verim * isp * G0)))"),
    ("Sıralı defterde kütle güncellenmiyor", "        m -= mp\n        satirlar.append", "        satirlar.append"),
    ("Kaçınma çifti tek manevra sayılıyor", 'v("cam_cift") * 2 * v("cam_ornek")', 'v("cam_cift") * v("cam_ornek")'),
    ("1A düzlem dışı itici de kayıplı sayılıyor", "return ip_k + yakit(m0 - ip_k, v(\"dv_1a_inc\"), isp)",
     "return ip_k + yakit(m0 - ip_k, v(\"dv_1a_inc\"), isp, 1 - k)"),
    ("Hohmann ikinci darbe unutuluyor", "+ abs(math.sqrt(MU / r2) * (1 - math.sqrt(r1 / a)))", ""),
    ("Bertaraf tasarım yakıtı 47 → 45 kg", '"yakit_eol": (47.0,', '"yakit_eol": (45.0,'),
    ("Kaynak kimliği defterde yok", '"S2-SP1322"),', '"S2-SP1322-YOK"),'),
]

shutil.copy(HESAP, YEDEK)
sonuc = []
try:
    for ad, eski, yeni in MUTASYONLAR:
        metin = YEDEK.read_text(encoding="utf-8")
        assert metin.count(eski) == 1, (ad, metin.count(eski))
        HESAP.write_text(metin.replace(eski, yeni), encoding="utf-8")
        p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                            "tests/test_sentinel_sunum.py"], cwd=ROOT, capture_output=True, check=False, text=True, encoding="utf-8", errors="replace",
                           env={**__import__("os").environ, "PYTHONIOENCODING": "utf-8"})
        son = [l for l in p.stdout.splitlines() if "passed" in l or "failed" in l][-1:]
        yakalandi = p.returncode != 0
        sonuc.append(yakalandi)
        print(f"{'YAKALANDI' if yakalandi else 'KAÇTI    '}  {ad:45s} {son[0] if son else p.stdout[-200:]}")
finally:
    shutil.copy(YEDEK, HESAP)
    YEDEK.unlink()
print(f"\n{sum(sonuc)}/{len(sonuc)} mutasyon yakalandı; hesap.py geri yüklendi")
