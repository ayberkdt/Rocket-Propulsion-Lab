# Sentinel-1C ve Sentinel-2C — itki odaklı görev araştırması (Revizyon V3)

UZM 508E / SASS kapsamında hazırlanan, Sentinel-1C ve Sentinel-2C görevlerinin rapor kontrol
listesine (A genel bilgi, B bütçeler, C CONOPS, D itki sistemi) göre incelenmesi.
Ders hazırlığı bağlamı: [../../SASS_PROJECT_PREPARATION.md](../../SASS_PROJECT_PREPARATION.md).

## Dosyalar

| Dosya | İçerik |
|---|---|
| [SENTINEL-1C.md](SENTINEL-1C.md) | Sentinel-1C: A–D başlıkları, kaynaklı değerler |
| [SENTINEL-2C.md](SENTINEL-2C.md) | Sentinel-2C: A–D başlıkları, gerçek manevra verisi |
| [HESAPLAR.md](HESAPLAR.md) | Tsiolkovsky, Hohmann, ömür sonu, yakıt, yanma süresi, tutulma, link, veri, rezerv |
| [KAYNAKLAR.md](KAYNAKLAR.md) | Kaynak dizini, her kaynağın doğrulama düzeyi, fotoğraf lisansları |
| [KONTROL_LISTESI.md](KONTROL_LISTESI.md) | Kontrol listesinin her maddesi için kanıt düzeyi |

## Doğrulama yöntemi

V2 paketi (29 Ekim 2026) üç bağımsız geçişle kontrol edildi:

1. **V2'nin kendi kaynak okuması.**
2. **Dört ayrı doğrulama okuması** (1C genel, 1 teknik/itki, 2C genel, 2 teknik/itki): her
   iddia birincil kaynakta arandı, kısa alıntı ve sayfa numarasıyla eşlendi.
3. **Birincil belgelerin bu çalışmada yeniden indirilmesi:** kritik sayıların geçtiği
   cümle ve tablolar doğrudan PDF'ten okundu (ESA SP-1322/1 ve /2, Rafael kataloğu,
   ArianeGroup veri sayfası, resmî VV24 kiti, ISSFD 2012/2017, SpaceOps 2025 #360 ve #617,
   CGMS-52) ve web sayfaları yeniden açıldı.

Bütün hesaplar bağımsız olarak yeniden yapıldı; V2'deki aritmetik hatasızdır.

## V2'ye göre başlıca düzeltmeler

- **VV24 kiti:** V2'nin kullandığı üçüncü taraf aynası eski sürümdür. Resmî kitte Zefiro-9
  ayrılması T+6:32 (V2: 5:52), planlanan ayrılma T+57:33; gerçekleşen 57:27'dir. Toplam yük
  1220 kg'dır.
- **"4075 W RF/DC kritik çelişkisi":** kanıtlanamadı; 4075 W'ı DC güç olarak veren kaynak
  bulunamadı. ESA SP-1322/1: 4368 W RF tepe, 3870 W DC (IW, tek veya çift polarizasyon).
- **SAR görev çevrimi:** 30 dk/tur yalnızca IW değil SM/IW/EW için geçerli; Mayıs 2026'ya kadar
  tutulma mevsiminde 25 dk idi.
- **Termal 13 → 7 °C:** uçuş ölçümünün değil, korelasyonlu ile korelasyonsuz termal model
  arasındaki belirsizlik farkıdır.
- **Rampa:** kaynaklar "Avrupa Uzay Üssü, Kourou" der; ZLV iddiası kaldırıldı.
- **LTDN 10:30 (S2):** ESA Facts and figures sayfasında değil, SentiWiki'de; eğim 98,62°.
- **"AstroBus-L":** doğrulanamadı; kaldırıldı.

## V3'te eklenen kaynaklı veri

- **S2C gerçek manevraları** (SpaceOps 2025 #617): nominalin 10,6 km altına enjeksiyon;
  ~6,33 m/s düzlem içi + 0,57 m/s düzlem dışı edinim; 6 Kasım 2024'te 314,1 s'de 0,934 m/s
  eğim manevrası; 18 Aralık 2024'e kadar toplam ~7,97 m/s. 2A rutin bakımı ~0,115 m/s/yıl.
- **S2 ömür sonu planı** (Airbus 2014): 786 → 653 km dairesel, 116 × 10 dk, 73 m/s, 47 kg,
  25 yılda yeniden giriş.
- **S2 itici** (ArianeGroup 1 N): Isp 220 s (200–223 s), 5,5–22 bar, 135 000 N·s.
- **S1 itki donanımı:** iki dal × yedi 1 N itici (SpaceOps 2025 #360); Rafael PEPT-590GB
  tankı (204 L, 4:1 blow-down'da 154 kg nominal kapasite, 27,4 → 5,5 bar).
- **S1 LEOP ΔV planlama tabloları** (ISSFD 2017) ve **1B bertarafının gerçek seyri**
  (SpaceOps 2025 #360, CGMS-52, ESA).
- **Güncel durum (29.10.2026):** 1A görevi 29 Haziran 2026'da bitti; 1C ve 1D operasyonel;
  1D'de 28 Ekim 2026'dan beri süren bir olay var; 2C'de 8 Ağustos 2025 veri kaybı; 2A
  kampanyası 2027 sonuna uzatıldı; 2D fırlatılmadı.

## Önemli uyarılar

- As-built kütle, güç, ısıl, link ve pointing bütçeleri kamuya açık değildir; bu belgeler
  aile değerlerini ve gerçek uçuş olaylarını açıkça ayırır, eksik bütçe uydurmaz.
- Senaryo hesapları (ör. Sentinel-1C rutin bakım tüketimi) ESA tahsisi değildir ve öyle
  işaretlenmiştir.
- Operasyon panosu gibi dinamik kaynaklar teslimden önce yeniden kontrol edilmelidir.
- Fotoğraf kullanımı için [KAYNAKLAR.md](KAYNAKLAR.md) lisans tablosuna bakın.

## Kaynaklar

Tam liste: [KAYNAKLAR.md](KAYNAKLAR.md). Örnek: ESA Sentinel-2 Facts and figures
<https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures>.
