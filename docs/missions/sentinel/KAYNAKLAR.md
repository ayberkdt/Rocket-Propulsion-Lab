# Kaynaklar ve doğrulama durumu (Revizyon V3)

**Kontrol tarihi:** 29 Ekim 2026. Her kaynak için **Doğrulama** sütunu şunu gösterir:

- **Ü (üçlü):** V2 paketi + bağımsız doğrulama okuması + birincil belgenin bu çalışmada
  yeniden indirilip ilgili cümle/tablonun okunması.
- **Ç (çift):** V2 paketi veya bir doğrulama okuması + bu çalışmada kaynağın yeniden
  açılması (alıntı karşılaştırıldı).
- **T (tek):** yalnızca bir okuma; raporda bu kaynağa dayanan cümle yoktur ya da açıkça
  belirtilmiştir.

Basın kitlerinde "planlanan" zamanlar ile gerçekleşen zamanlar ayrı tutulmuştur. Dinamik
sayfalar (operasyon panosu, haber listeleri) zamanla değişebilir.

## Sentinel-1 / Sentinel-1C

| Kod | Kaynak | Desteklediği bilgi | Doğrulama |
|---|---|---|---|
| S1-ESA | ESA, Sentinel-1 Facts and figures — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures | görev, 693 km, 2185 kg, boyutlar, yükleniciler, AIS, fırlatma tarihleri, ömür | Ü |
| S1-VV25 | Arianespace VV25 fırlatma kiti (yalnızca üçüncü taraf aynası) — https://sky-brokers.com/wp-content/uploads/2025/07/Arianespace-Launchkit-VV25-ESA-Copernicus-Sentinel-1c.pdf | 2186 / 2286 kg, 700 km, 98,19°, kalkış 21:20:33 UTC, olay zamanları, ömür ifadesi | Ç — resmî Arianespace kopyası bulunamadı (VV25 sayfası, medya arşivi, Wayback); 700 km ve T+1 sa 43 dk Arianespace bültenleriyle uyumlu |
| S1-ARI-RESULT | Arianespace, "Arianespace successfully launches Sentinel-1C…" — https://www.arianespace.com/news/arianespace-successfully-launches-sentinel-1c-earth-observation-satellite-for-european-unions-copernicus-program/ | 18:20 Kourou / 21:20 UTC, ~700 km, ayrılma kalkıştan 1 sa 43 dk sonra; Vega-C'nin uçuşa dönüşü | Ç |
| S1-LAUNCH | ESA, "Double win for Europe: Sentinel-1C and Vega-C take to the skies" — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Double_win_for_Europe_Sentinel-1C_and_Vega-C_take_to_the_skies | ayrılma ~00:04 CET, iletişim 00:12 CET | Ç |
| S1-ESOC-MEDIA | ESA medya daveti — https://www.esa.int/Newsroom/Press_Releases/Media_invitation_Sentinel-1C_launch_event_in_Darmstadt_Germany | "the ESA control centre in Darmstadt will take over the satellite in orbit" | Ç |
| S1-FUEL | ESA, "Sentinel-1C fuelled ahead of liftoff" (18.11.2024) — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff | 154 kg yakıt | Ü |
| S1-FUEL-IMG | ESA görseli "Fuelling Sentinel-1C" — https://www.esa.int/ESA_Multimedia/Images/2024/11/Fuelling_Sentinel-1C | "fuelling the satellite with hydrazine", "critical and hazardous" | Ç |
| S1-SP1322 | ESA SP-1322/1, *Sentinel-1: ESA's Radar Observatory Mission for GMES Operational Services* (2012) — https://esamultimedia.esa.int/multimedia/publications/SP-1322_1/offline/download.pdf | Tablo 4.1 (2300 kg / 130 kg, 7 yıl + 12 yıllık sarf, 5900 W, 324 Ah, 1410 Gbit), Tablo 5.2 (4368 W RF, 3870 W DC, 945 kg), PRIMA, itki işlevleri, itki modülü, modlar | Ü (resmî ESA kopyası; PDF s.24–30) |
| S1-WIKI | Copernicus SentiWiki, S1 Mission — https://sentiwiki.copernicus.eu/web/s1-mission | 360 Ah, SAR görev çevrimi ve tutulma mevsimi, yörünge tüpü, link hızları, tekrar döngüsü | Ç |
| S1-ISSFD | Vasconcelos vd., ESA/ESOC, ISTS-ISSFD 2017 — https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf | S1B: 3 çift OCT + 4 çift RCT, yüzeyler, 35 s / 600 s / 6 sa, 0,28 m/s, dağılım ΔV tabloları, enjeksiyon hatası düzeltmeleri, −%30 plume kaybı | Ü |
| S1-ISSFD2012 | ESA/ESOC, ISSFD 2012 FDOP2-4 — https://www.issfd.org/ISSFD_2012/ISSFD23_FDOP2_4.pdf | 60 m ölü bant; IP manevraları sürüklenme; OOP sıklığı > haftada 1 (Ay) | Ç |
| S1-SPACEOPS360 | SpaceOps 2025, ID 360 — https://publications.spaceops.org/2025/download_by_id.php?id=0360 | iki dal × yedi 1 N itici; 1B bertarafı (Şubat 2023, Nisan 2023 −Z iticiler, 29 Şubat 2024), 0,47–0,58 N, ~5 kg yakıt tahmin belirsizliği | Ü |
| S1-RAFAEL | Rafael, Space Propulsion teknik kataloğu (2024) — https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf | basılı s.9 (PEPT-590GB, ST1 A–D), s.10 (28 + 28 itici), s.41 (tank özellikleri), s.46 (1 N itici tablosu) | Ü — kapakta "proprietary" ibaresi var; yalnızca sayısal veriye atıf yapıldı |
| S1-TMM-2026 | ESA European Space Thermal Engineering Workshop 2026, Latella (Thales Alenia Space) — https://indico.esa.int/event/643/contributions/12724/ | 1C termal model korelasyonu, %98 güvende 13 → 7 °C | Ç (yalnızca özet okunabildi; sunum giriş gerektiriyor) |
| S1-CGMS-2024 | CGMS-52 ESA çalışma belgesi (Nisan 2024) — https://www.cgms-info.org/Agendas/GetWpFile.ashx?aid=6e1f8b19-7eb6-4356-bbe6-0673c950ba5e&wid=da80287c-7a8a-44bd-8c5f-3e2499564b22 | s.25: 1B bertaraf planı (< 530 km dairesel, ~12 yıl), ~88 kg kalan yakıt | Ü |
| S1B-DISPOSAL | ESA, "Sentinel-1B journeys back to Earth" (23.09.2024) — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1B_journeys_back_to_Earth | 3 Ağustos 2022 görev sonu, 12 Eylül 2024 pasivasyon, ~24 yıl | Ç |
| S1A-END | ESA, "Time to say goodbye to Sentinel-1A" (30.06.2026) — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Time_to_say_goodbye_to_Sentinel-1A | 29 Haziran 2026 görev sonu; "1C and 1D performing very well" | Ü |
| S1D-LIVE | ESA (Almanca), "Sentinel-1D ist einsatzbereit" (01.05.2026) — https://www.esa.int/Space_in_Member_States/Germany/Sentinel-1D_ist_einsatzbereit_ein_Meilenstein_fuer_Europas_Radarmission | 1D tam operasyonel | Ç |
| S1-DASH | Copernicus Operations Dashboard — https://operations.dashboard.copernicus.eu/index.html | 1D olayları 26–27 Ekim ve 28 Ekim 2026 (dinamik) | T (değişken sayfa; teslimden önce yeniden bakılmalı) |
| S1-FIRST-IMG | ESA, "Sentinel-1C captures first radar images" — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_captures_first_radar_images | ilk görüntü 56 sa 23 dk sonra, Svalbard | Ç |
| S1-NODE | CDSE, 30.01.2025 — https://dataspace.copernicus.eu/node/2103 | Ocak 2025 düğüm manevraları, devreye alma | Ç |
| S1-OPEN | CDSE, "Sentinel-1C User Data Opening from 26th of March" — https://dataspace.copernicus.eu/node/2232 | 26 Mart 2025 veri açılışı | Ç |
| S1-CAL-2025 | Copernicus, "Sentinel-1C products are now calibrated" — https://sentinels.copernicus.eu/-/sentinel-1c-products-are-now-calibrated | 19 Mayıs 2025 | Ç |
| S1-ANOM-2026 | CDSE, "Copernicus Sentinel-1A/1C/2B return to nominal operations" (20.01.2026) — https://dataspace.copernicus.eu/node/2958 | 18.01 anomalisi; 19.01 19:40 UTC nominal | Ç |
| S1-CAL-2026 | Copernicus, "Update of Sentinel-1C radiometric calibration" — https://sentinels.copernicus.eu/-/update-of-sentinel-1c-radiometric-calibration | 2025 TRM arızaları, 3 Şubat 2026 güncellemesi | Ç |
| S1-RECONF | CDSE, 28.05.2026 — https://dataspace.copernicus.eu/node/3949 | 9–23 Haziran 2026 manevraları | Ç |
| S1-CONST | Copernicus, "Sentinel-1C/1D final orbital configuration achieved" — https://sentinels.copernicus.eu/-/sentinel-1c/1d-final-orbital-configuration-achieved | 24 Haziran 2026, 6 günlük tekrar | Ç |

## Sentinel-2 / Sentinel-2C

| Kod | Kaynak | Desteklediği bilgi | Doğrulama |
|---|---|---|---|
| S2-ESA | ESA, Sentinel-2 Facts and figures — https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures | MSI 13 bant/290 km, boyutlar, 1140 kg (123 kg yakıt), 786 km, yükleniciler, Kiruna, EDRS | Ç |
| S2-VV24 | Arianespace, resmî VV24 fırlatma kiti — https://newsroom.arianespace.com/vega-flight-vv24 | s.2 toplam yük 1220 kg, 775 km, 98,57°, planlanan kalkış; s.3 1143 kg, ömür; s.4–5 kademeler ve olay zamanları (Z9 T+6:32, ayrılma T+57:33) | Ü (V2'nin kullandığı ayna eski sürümdür) |
| S2-ARI-RESULT | Arianespace sonuç bülteni — https://newsroom.arianespace.com/arianespace-successfully-launches-europes-copernicus-earth-observation-program-sentinel-2c-satellite/ | 5 Eylül 01:50 UTC, ayrılma 57 dk 27 s, 22. ve son Vega | Ü |
| S2-ESA-LEOP | ESA, "A smooth start to life in orbit for Sentinel-2C" (10.09.2024) — https://www.esa.int/Enabling_Support/Operations/A_smooth_start_to_life_in_orbit_for_Sentinel-2C | ayrılma ~04:48 CEST, ilk sinyal 05:02 CEST Troll, LEOP sonu ve devreye alma başlangıcı | Ü |
| S2-DELAY | ESA blog, "Launch update" (04.09.2024) — https://blogs.esa.int/eolaunches/2024/09/04/launch-update/ | "Owing to electrical issues… postponed"; yeni deneme 5 Eylül 03:50 CEST | Ç |
| S2-FUEL | ESA, "Sentinel-2C fully loaded" (19.08.2024) — https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/ ; görsel: https://www.esa.int/ESA_Multimedia/Images/2024/08/Fuelling_Sentinel-2C | 16 Ağustos: 133 kg hidrazin, 21 bar He, latch valfleri | Ü |
| S2-SP1322 | ESA SP-1322/2, *Sentinel-2: ESA's Optical High-Resolution Mission for GMES Operational Services* (2012) — https://sentinels.copernicus.eu/documents/247904/349490/S2_SP-1322_2.pdf | 1225 kg (70 kg marj), 117 kg hidrazin, 1700 W EOL / 7,2 m², 87 Ah, 1,4 kW, 2,4 Tbit, 560 Mbit/s 8-PSK isoflux, OCP (Tesat), sensörlerin teleskop yapısına montajı, itkinin halka içine yerleşimi, itki işlevleri | Ü |
| S2-AOCS | Wiedermann, Gockel, Winkler, Rieber, Kraft, Reggio (Airbus DS), "The Sentinel-2 Satellite Attitude Control System – Challenges and Solutions", 9th ESA GNC Conference, Porto, 2014 — kayıt: https://www.esa-gnc.eu/paper/?paper_id=830 (PDF giriş gerektiriyor); açık yazar kopyası: https://www.researchgate.net/publication/271584776 | 2 × 4 eğik 1 N itici, CHT1N-6, blow-down tank, itki modülü, yönelim gereksinimleri, AOCS donanımı, dört iticinin eşzamanlı ateşlenmesi, ömür sonu planı (653 km, 116 × 10 dk, 73 m/s, 47 kg, 25 yıl) | Ç (ResearchGate otomatik indirmeyi engelliyor; metin tarayıcıyla okunup alıntılar ayrıca karşılaştırıldı) |
| S2-AG-THR | ArianeGroup, Hydrazine Thrusters broşürü — https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf | 1 N itici: 0,32–1,1 N, Isp 220 s (200–223 s), 5,5–22 bar, 135 000 N·s, 290 g; miras listesinde Sentinel-2A/2B | Ü (s.6) |
| S2-SPACEOPS617 | Kuijper, Pol-Moreno, Affaitati, Ferrari, SpaceOps 2025, ID 617 — https://publications.spaceops.org/2025/download_by_id.php?id=0617 | 2A yörünge kontrolü (64 manevra, 1,56 cm/s ort., 24 çift CAM), ±2 km / ±10 s, 2C enjeksiyonu (−10,6 km), 6,33 + 0,57 m/s, Tablo 3–4 | Ü |
| S2-WIKI | Copernicus SentiWiki, S2 Mission — https://sentiwiki.copernicus.eu/web/s2-mission | 98,62°, LTDN 10:30, 7,25 yıl / 12 yıl, 2A'nın 36° uzaklığı | Ç |
| S2-OPS | ESA, "Sentinel-2 operations" (2017, artık güncellenmiyor) — https://www.esa.int/Enabling_Support/Operations/Sentinel-2_operations | 2300 W başlangıç, 120/123 kg, 1,6 TB çelişkisi, Kiruna, ESOC | Ç |
| S2-CDSE-DOC | CDSE belgesi, Sentinel-2 — https://documentation.dataspace.copernicus.eu/Data/SentinelMissions/Sentinel2.html | 21 Ocak 2025'te 2A'nın yerini aldı; 5/10 gün tekrar | Ç |
| S2-COMMISSION | CDSE, 08.01.2025 — https://dataspace.copernicus.eu/node/2039 | devreye alma Aralık 2024 sonunda tamamlandı | Ç |
| S2-TRANSFER | CDSE, 17.01.2025 — https://dataspace.copernicus.eu/node/2061 | planlanan ilk algılama 09:12:45 UTC | Ç |
| S2-ANOMALY | CDSE, "Sentinel-2C Unavailability on 08 August 2025" — https://dataspace.copernicus.eu/node/2594 | %36,25 veri kaybı, 8–11 Ağustos 2025 | Ç |
| S2-DECON | CDSE — https://dataspace.copernicus.eu/node/2642 | 8–9 Eylül 2025 MSI dekontaminasyonu (sayfa başlığı) | Ç |
| S2-EXT | CDSE, 09.10.2026 — https://dataspace.copernicus.eu/news/2026-10-9-sentinel-2a-extension-campaign-continue-until-end-2027 | 2A kampanyası 31 Aralık 2027'ye kadar | Ç |
| S2-OSCAR-MSI | WMO OSCAR, MSI (Sentinel-2A) — https://space.oscar.wmo.int/instruments/view/msi_sentinel_2a | 275 kg, 266 W, 450 Mbit/s (2A) | Ç |
| S2D-OSCAR | WMO OSCAR, Sentinel-2D — https://space.oscar.wmo.int/satellites/view/sentinel_2d | planlanan, ≥ 2028 | T |

## Teknoloji olgunluğu ve enkaz

| Kod | Kaynak | Doğrulama |
|---|---|---|
| TRL-ESA | ESA, Technology Readiness Levels — https://www.esa.int/Enabling_Support/Space_Engineering_Technology/Shaping_the_Future/Technology_Readiness_Levels_TRL | V2 |
| DEBRIS-ESA | ESA, yeni enkaz azaltma politikası (2023) — https://esoc.esa.int/new-space-debris-mitigation-policy-and-requirements-effect | V2 |
| DEBRIS-CHARTER | ESA, Zero Debris Charter SSS — https://www.esa.int/Space_Safety/Clean_Space/Zero_Debris_Charter_-_Frequently_Asked_Questions | V2 |

## Fotoğraflar ve lisans

| Görev | Sayfa | Kredi / lisans |
|---|---|---|
| Sentinel-1C | https://www.esa.int/ESA_Multimedia/Images/2023/01/Sentinel-1C_undergoing_tests | Thales Alenia Space/ImagIN; ESA Standard Licence, ek izin gerekebilir |
| Sentinel-2C | https://www.esa.int/ESA_Multimedia/Images/2024/07/Sentinel-2C_ready_to_be_loaded_into_its_transport_container | Airbus; "No ESA licences available" |
| Sentinel-2C (alternatif) | https://www.esa.int/ESA_Multimedia/Images/2024/08/Fuelling_Sentinel-2C | ESA-CNES-Arianespace / Optique vidéo du CSG; ESA Standard Licence |

## V2 kaynaklarından çıkarılanlar ve düzeltilenler

- **"Sentinel-1 – the future GMES C-band SAR mission" (4075 W DC):** erişilebilir kopyası
  bulunamadı; 4075 W'ı DC güç olarak veren başka kaynak yok. Rapor 3870 W DC kullanır.
- **sky-brokers VV24 aynası:** resmî kitten eski sürüm; resmî kitle değiştirildi.
- **WMO OSCAR Sentinel-1C/2C (2300/2146 kg; 1200/1080 kg):** ESA ve Arianespace verisiyle
  çelişir ve eski kayıtlardır; kütle kaynağı olarak kullanılmadı.
- **ESA blog "Sentinel-2C launched into orbit":** tarihi hatalı ("4 Eylül 03:50 CEST");
  kullanılmadı.
- **AstroBus-L adı:** kaynaklarda doğrulanamadı; çıkarıldı.

## Kaynaklar

Bu dosya kaynak dizinidir; bütün bağlantılar yukarıdadır.
