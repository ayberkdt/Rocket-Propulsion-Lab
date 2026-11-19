# Sentinel-1/2 itki sunumu

Güncel sürüm: **v32** → `surumler/Sentinel_1_2_Propulsion_ITU_v32.pptx`

Eski sürümler silinmez. Her değişiklik yeni bir numarayla `surumler/` altına yazılır.

| Sürüm | Tarih | Değişiklik | Üreten betik |
|---|---|---|---|
| v21 | 2026-11-12 | ChatGPT'nin teslim ettiği 9 slaytlık İTÜ destesi (temel) | — |
| v22 | 2026-11-18 | Slayt 9 (Propellant Budget): Tsiolkovsky tahmini ve S2C yakıt grafiği. Slayt 1–8 v21 ile aynı. | `uretim/v22_slayt9.py` |
| v23 | 2026-11-18 | Slayt 1, 2, 3, 6, 7 ve 9'da görseller büyütüldü, alt boşluklar kapatıldı. Slayt 4, 5, 8 v22 ile aynı. | `uretim/v23_yerlesim.py` |
| v24 | 2026-11-18 | Slayt 3: S1A açılış sırası (ISSFD 2015 #137) ve S1B bertaraf izi (SpaceOps 2025 #226). Slayt 5: S1B itici yerleşimi (SpaceOps 2025 #360). Slayt 6: ArianeGroup 1 N fotoğrafı. Diğerleri v23 ile aynı. | `uretim/v24_gorseller.py` |
| v25 | 2026-11-18 | Hoca değerlendirmesinin zorunlu düzeltmeleri: 2A uzatması end-2026 (end-2027 değil); filo, bırakılış ve CONOPS satırlarına kaynak; slayt 7'de doğrulanamayan tork grafiği yerine S1B itici yönleri ve "yaklaşık 10×"; slayt 8'de 7.97 m/s "planlanan"; slayt 9 tek ondalık ve 47 kg / 38 kg payı. Okur metni 420 kelime (`ozel/okur_metni_v25.md`). | `uretim/v25_duzeltmeler.py` |
| v26 | 2026-11-18 | Ödevin eksik maddeleri birer satır: D9 güvenlik (slayt 2) ve 25 yıl kuralı (slayt 3), D10 TRL 9 (slayt 6), S1A roket denklemi 6,3 m/s → ideal 6,7 kg, uçuşta 8,4 kg (slayt 7), D5 debi ve yanma süresi, D8 anlık yanma (slayt 9). Okur metni `ozel/okur_metni_v26.md`. | `uretim/v26_odev_maddeleri.py` |
| v27 | 2026-11-18 | Dış geri bildirim, birincil kaynakla doğrulanarak (11 slayt). 2A uzatması end-2027 (CDSE, 9 Eki 2026; v25'teki end-2026 yanlıştı). S2 uçuş mirası 2015'ten (S2A 23 Haz 2015). Slayt 7: 6,3 m/s enjeksiyon düzeltmesi, m₀ 2157 kg (VS07), Isp 214 s → ideal 6,5 kg; IP itkisi −%21…−%30 → 7,6–8,3 kg; uçuşta 8,4 kg (ISSFD 2015 #134/#137). Slayt 8: 0,57 m/s "planlanan, icrası bildirilmemiş" (SpaceOps #617 s.16). Yeni slayt 9: S2C 12 yıllık ΔV bütçesi (≈85 m/s, bertaraf %86). Slayt 10: 82,1 kg "kalan" yerine tahsis şelalesi, son çubuk "Rest" (rezerv değil); 3,4 N net itki ve 19° yay / %0,45 kayıp. Yeni slayt 11: tasarım etkenleri. SASS alt başlığı ve dolgu cümleler silindi. Okur metni `ozel/okur_metni_v27.md`. | `uretim/v27_geribildirim.py` |
| v28 | 2026-11-18 | İkinci geri bildirim. Slayt 3 yeniden tasarlandı: iki küçük makale figürü yerine 5 evreli şerit ve 4 satır (itkinin işi, Sentinel-1 ve Sentinel-2 örneği, 2C ΔV). Slayt 7: 7,6–8,3 kg "simplified estimate". Slayt 9: "Illustrative 12-year budget", Basis sütunu durumla başlıyor (flown / planned / projected / A/B design), "Additional inclination control, 12 yr". Slayt 10: yatay şelale, "Unallocated (scenario)", 3,4 N "average effective force", sonlu yanma cezası varsayımıyla. Slayt 11: iddialar düz teknik ifade. Okur metni `ozel/okur_metni_v28.md` (618 kelime). | `uretim/v28_conops.py` |
| v29 | 2026-11-18 | Deste itki merkezli ve tek şablonla baştan kuruldu (10 slayt): Mission Overview (Table 1), Spacecraft and Launch, Propulsion Architecture (Table 2), Feed System and Layout, Thruster Performance (Table 3), Concept of Operations, Flight Experience (Table 4), ΔV Budget (Table 5), Propellant Budget, Conclusions. KPI kartları, iddia/kanıt satırları, 6 Kasım ve 23 Kasım manevra ayrıntıları, Spacecraft Systems ve Sentinel-2C Δv slaytları çıktı. Planlanan 0,57 m/s bütçeden çıktı. Okur metni `ozel/okur_metni_v29.md` (364 kelime). | `uretim/v29_yapi.py` |
| v30 | 2026-11-18 | Görsel kalite: CONOPS dikeyde dengelendi ve 16 pt'ye çıktı; Flight Experience iki vaka sütunu (başlık bandı, kaynak figürü: ISSFD 2015 #134 Fig. 5 ve SpaceOps 2025 #226 Fig. 9, Problem/Fix/Cost satırları); Conclusions numaralı dört çıkarım ve dolum fotoğrafı. ΔV tablosu Table 4. İçerik v29 ile aynı. Okur metni `ozel/okur_metni_v30.md`. | `uretim/v30_kalite.py` |
| v31 | 2026-11-19 | Slayt 7: Sentinel-1B'nin duruş nedeni düzeltildi. v30 duruşu teker aşırı hızına bağlıyordu; SpaceOps 2025 #226 Fig. 9'a göre kampanya Nisan 2024'te RCT-3B iticisinin düşük performansından sonra askıya alındı. Okur metni `ozel/okur_metni_v31.md`; çalışma kağıtları `ozel/calisma/calisma_kagidi.pdf`. | `uretim/v31_s1b_neden.py` |
| v32 | 2026-11-19 | Slayt 7: 1B için "passivated" yerine "Under 25 years; electrically passivated only". Elektriksel pasivasyon Eylül 2024'te yapıldı, yakıt kaldığı için tam pasivasyon sağlanamadı (SENTINEL-1C.md, SpaceOps 2025 #226 s.20). | `uretim/v32_pasivasyon.py` |

Hiçbir sürümde okur metni (konuşmacı notu) yok; hepsi `uretim/notlari_ayir.py` ile `ozel/okur_metni_vNN.md` dosyalarına taşındı.

## Klasörler

- `surumler/`: sunum dosyaları, sürüm numaralı.
- `uretim/`: sürümleri üreten betikler. `uretim/grafik/` altında grafik ve datawiz denetim raporları. `notlari_ayir.py` bir sunumun notlarını `ozel/`e taşır.
- `uretim/gorsel_adaylari/`: kaynak makalelerden çıkarılan figürler (A–F) ve ArianeGroup fotoğrafı (G).
- `ozel/`: okur metni (konuşmacı notları) ve `calisma/` altında çalışma kağıtları (XeLaTeX, iki kez derlenir). Git'e girmez, hocaya giden dosyada yoktur.

## Yeniden üretmek

```bash
python uretim/grafik_s2c_yakit.py v23
python uretim/v23_yerlesim.py
python uretim/v24_gorseller.py
python uretim/grafik_s2c_yakit.py v25
python uretim/v25_duzeltmeler.py
python uretim/v26_odev_maddeleri.py
python uretim/grafik_s2c_yakit.py v27
python uretim/v27_geribildirim.py
python uretim/grafik_s2c_yakit.py v28
python uretim/v28_conops.py
python uretim/grafik_s2c_yakit.py v29
python uretim/v29_yapi.py
python uretim/v30_kalite.py
python uretim/v31_s1b_neden.py
python uretim/v32_pasivasyon.py
```

Yeni sürüm için `v23_yerlesim.py` dosyasını kopyalayın, girdi olarak bir önceki sürümü, çıktı olarak bir sonraki numarayı verin ve bu tabloya bir satır ekleyin.
