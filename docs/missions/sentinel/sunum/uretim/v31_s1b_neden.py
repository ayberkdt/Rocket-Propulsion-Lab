"""v31: slayt 7'de Sentinel-1B bertarafının duruş nedeni düzeltilir.

v30 "Problem: wheel overspeed and several aborted burns → Result: lowering stopped" diyordu; bu, duruşu teker
aşırı hızına bağlıyordu. SpaceOps 2025 #226 Fig. 9 notu: "Apr 2024 Disposal suspended after RCT-3B
underperformance". Teker aşırı hızı (Temmuz 2023) yarıda kalan yanmalardan biriydi.
Girdi: surumler/..._v30.pptx · Çıktı: surumler/..._v31.pptx · Başka değişiklik yok.
"""
from pathlib import Path

from pptx import Presentation

SUNUM = Path(__file__).resolve().parent.parent
SRC = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v30.pptx"
OUT = SUNUM / "surumler" / "Sentinel_1_2_Propulsion_ITU_v31.pptx"

prs = Presentation(SRC)
s = list(prs.slides)[6]
REPS = {
    "Wheel overspeed and several aborted burns": "Aborted burns; RCT-3B thruster underperformed",
    "Lowering stopped at ~580 × 610 km": "Campaign stopped in April 2024 at ~580 × 610 km",
}
done = set()
for sh in s.shapes:
    if sh.has_text_frame:
        for p in sh.text_frame.paragraphs:
            for r in p.runs:
                if r.text in REPS:
                    done.add(r.text)
                    r.text = REPS[r.text]
assert done == set(REPS), done
prs.save(OUT)
print("wrote", OUT)
