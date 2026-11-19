"""v32: slayt 7'de Sentinel-1B pasivasyonu "elektriksel" olarak düzeltilir.

v31 "Predicted within 25 years; passivated" diyordu ve tam pasivasyon gibi okunuyordu. Depodaki SENTINEL-1C.md
(SpaceOps 2025 #226 s.20'ye dayanarak): elektriksel pasivasyon Eylül 2024'te yapıldı; yakıt kaldığı için tam
pasivasyon sağlanamadı.
Girdi: surumler/..._v31.pptx · Çıktı: surumler/..._v32.pptx · Başka değişiklik yok.
"""
from pathlib import Path

from pptx import Presentation

SUNUM = Path(__file__).resolve().parent.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v31.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v32.pptx"

OLD = "Predicted within 25 years; passivated"
NEW = "Under 25 years; electrically passivated only"

prs = Presentation(SRC)
s = list(prs.slides)[6]
hits = 0
for sh in s.shapes:
    if sh.has_text_frame:
        for p in sh.text_frame.paragraphs:
            for r in p.runs:
                if r.text == OLD:
                    r.text = NEW
                    hits += 1
assert hits == 1, hits
prs.save(OUT)
print("wrote", OUT)
