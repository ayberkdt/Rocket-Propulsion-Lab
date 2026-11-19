"""kanitlar.json'daki alıntıları indirilmiş birincil belgelerin metnine karşı denetler.

Kullanım: python kanit_dogrula.py <metin_klasoru>
Klasörde her belgenin 'yerel' adlı düz metni bulunur (PDF'ten '==== PAGE n' işaretli). Alıntı, boşluklar
sadeleştirilerek önce belirtilen sayfada, bulunamazsa bütün belgede aranır. Çıktı: her alıntı için
BULUNDU (sayfada) / BULUNDU (başka sayfada) / YOK; YOK varsa çıkış kodu 1.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sade(s: str) -> str:
    s = s.replace("−", "-").replace("–", "-").replace("’", "'")
    return re.sub(r"\s+", " ", s).strip()


def sayfalar(metin: str) -> dict[int, str]:
    parcalar = re.split(r"\n==== PAGE (\d+)\n", metin)
    if len(parcalar) == 1:
        return {0: sade(metin)}
    return {int(parcalar[i]): sade(parcalar[i + 1]) for i in range(1, len(parcalar) - 1, 2)}


def main(klasor: Path) -> int:
    defter = json.loads((HERE / "kanitlar.json").read_text(encoding="utf-8"))
    belgeler = defter["belgeler"]
    eksik = 0
    for k in defter["kanitlar"]:
        if k["durum"] != "birincil":
            continue
        yol = klasor / belgeler[k["belge"]]["yerel"]
        if not yol.exists():
            print(f"{k['id']:14s} BELGE YOK: {yol.name}")
            eksik += 1
            continue
        sf = sayfalar(yol.read_text(encoding="utf-8"))
        m = re.search(r"\d+", k["sayfa"])
        hedef = int(m.group()) if m and 0 not in sf else None
        for a in k["alintilar"]:
            a2 = sade(a)
            if hedef is not None and a2 in sf.get(hedef, ""):
                durum = f"BULUNDU (s. {hedef})"
            else:
                yerler = [n for n, t in sf.items() if a2 in t]
                if yerler:
                    durum = f"BULUNDU (başka sayfa: {yerler[:3]})" if hedef is not None else "BULUNDU"
                else:
                    durum = "YOK"
                    eksik += 1
            print(f"{k['id']:14s} {durum:28s} «{a[:70]}»")
    print(f"\nsonuç: {'TAMAM' if not eksik else f'{eksik} alıntı bulunamadı'}")
    return 1 if eksik else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
