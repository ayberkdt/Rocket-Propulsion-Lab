# Rapor kontrol listesi — kapsam matrisi (Revizyon V3)

**Kodlar:** **D** = görevin kendisine özgü birincil kaynakla doğrulandı · **M** = aile/miras
verisi (A/B, tasarım belgesi) · **H** = bu çalışmanın kontrol edilmiş hesabı · **AÇIK** =
kamuya açık kaynakta bulunamadı. "Başlık var" tamamlandı anlamına gelmez.

| Madde | Sentinel-1C | Sentinel-2C | Nerede |
|---|---|---|---|
| **A** Görev tanımı, operatör, üretici | D | D/M | 1C A1 · 2C A1 |
| A Fotoğraf | D (ESA sayfası; ESA Standard Licence, ek izin gerekebilir) | D (Airbus kredili, ESA lisansı yok; alternatif dolum fotoğrafı ESA Standard Licence) | KAYNAKLAR |
| A Kütle, boyut, güç | D (2186 kg, 154 kg) + M (boyut, güç) | D (1143 kg, 133 kg) + M | A2 |
| A Yörünge, tasarım ömrü | D (hedef 700 km; ömür ifadesi) + M | D (hedef 775 km; ömür ifadesi) + M | A2 |
| A Fırlatıcı, tarih, rampa, enjeksiyon | D | D | A3 |
| A Operasyonel durum ve tarih | D (29.10.2026) | D (29.10.2026) | A4 |
| A Tüm görev sistemi (fırlatıcı + üst kademe) | D (Vega-C P120C/Z40/Z9/AVUM+, Vampire) | D (Vega P80/Z23/Z9/AVUM) | A3 |
| A Genel tartışma | var | var | A4 |
| **B** Kütle bütçesi | M + AÇIK (alt kalem) | M (70 kg marj) + AÇIK | B1 |
| B Güç bütçesi | M (5900 W, 3870 W SAR) | M (1700 W, 1,4 kW, 87 Ah) | B2 |
| B Pointing bütçesi | M (gereksinimler) + AÇIK (1C hata bütçesi) | M (Airbus gereksinimleri) + AÇIK | B3 |
| B Link bütçesi | M (hızlar) + AÇIK (marj) | M (hızlar, modülasyon) + AÇIK (marj) | B4 |
| B ΔV bütçesi | M (S1B LEOP tabloları) + H (EOL) + AÇIK (1C gerçek değerleri) | **D (gerçekleşen 7,97 m/s)** + M (EOL 73 m/s) + H | B7 · HESAPLAR §4–5 |
| B Veri bütçesi | M + H | M + H | B5 |
| B Isıl bütçe | D (2026 model korelasyonu) + AÇIK (limitler) | AÇIK (limitler) | B6 |
| B Faydalı yük ve alt sistemler | M (C-SAR) + D (AIS) | M (MSI) | B8 |
| **C** CONOPS zaman çizelgesi | D (2024–2026 olayları) + M (1B bertarafı) | D (2024–2025 olayları, manevralar) + M (EOL planı) | C |
| C CONOPS'a göre ΔV | H (senaryo) + AÇIK | **D + H** | HESAPLAR §5 |
| **D** İtki tipi | D (hidrazin) | D (hidrazin + He blow-down) | D1 |
| D Kütle, güç, boyut, yerleşim | M (Rafael tank ve itici, itki modülü) | M (ArianeGroup itici, halka içi modül) | D1–D2 |
| D İşlev | M | M | D3 |
| D Tsiolkovsky (ΔV, kütle, itki) | H (kaynaklı Isp aralığı) | H (kaynaklı Isp) | D4 · HESAPLAR §2 |
| D Yanma süresi, kütle debisi | H + M (35 s / 600 s sınırları) | H + D (314 s eğim manevrası, 3,4 N etkin itki) | D4 · HESAPLAR §6 |
| D İtici kütlesi hesabı | H | H + D (4,2 kg kullanılan) | HESAPLAR §5 |
| D Rezerv ve ek yıl | H (senaryo) + AÇIK (gerçek kalan) | H (kalan ~75–81 kg) — yıl sayısı verilmez | D6 · HESAPLAR §12 |
| D Anlık / sürekli | H + M | H + D | D5 |
| D Güvenlik ve regülasyon | D (dolum) + M (1B bertarafı) | D (dolum) + M (EOL planı) | D7 |
| D TRL (önce/sonra) | M (uçuş mirası) + AÇIK (resmî kayıt) | M (uçuş mirası) + AÇIK | D8 |

**V2'ye göre seviye atlayan maddeler:** S2C ΔV bütçesi ve CONOPS ΔV'si artık gerçek uçuş
manevralarına dayanıyor; S2 ömür sonu (73 m/s, 47 kg) ve Isp (220 s, 200–223 s) kaynaklı;
S1 tank (PEPT-590GB) ve itici sayısı (2 × 7) kaynaklı; S1 ömür sonu hedefi 1B'nin gerçek
bertaraf hedefinden alındı.

**Hâlâ açık:** as-built kütle/güç/ısıl/link bütçeleri, 1C'nin gerçek manevra büyüklükleri ve
kalan yakıtı, resmî TRL değerlendirmeleri. Bunlar ESA veya üreticiden belge talebiyle
kapanabilir; kamuya açık kaynakta yoktur.

## Kaynaklar

Bütün kaynaklar ve doğrulama durumları: [KAYNAKLAR.md](KAYNAKLAR.md). Örnek birincil kaynak:
<https://publications.spaceops.org/2025/download_by_id.php?id=0617>.
