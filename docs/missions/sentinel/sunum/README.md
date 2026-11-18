# Sentinel-1/2 itki sunumu

Güncel sürüm: **v27** → `surumler/Sentinel_1_2_Propulsion_ITU_v27.pptx`

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

Hiçbir sürümde okur metni (konuşmacı notu) yok; hepsi `uretim/notlari_ayir.py` ile `ozel/okur_metni_vNN.md` dosyalarına taşındı.

## Klasörler

- `surumler/`: sunum dosyaları, sürüm numaralı.
- `uretim/`: sürümleri üreten betikler. `uretim/grafik/` altında grafik ve datawiz denetim raporları. `notlari_ayir.py` bir sunumun notlarını `ozel/`e taşır.
- `uretim/gorsel_adaylari/`: kaynak makalelerden çıkarılan figürler (A–F) ve ArianeGroup fotoğrafı (G).
- `ozel/`: okur metni (konuşmacı notları). Git'e girmez, hocaya giden dosyada yoktur.

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
```

Yeni sürüm için `v23_yerlesim.py` dosyasını kopyalayın, girdi olarak bir önceki sürümü, çıktı olarak bir sonraki numarayı verin ve bu tabloya bir satır ekleyin.
