# Sentinel-1C — ayrıntılı C örneği (V4 düzeltmesi)

**Ana kapsam artık seriler karşılaştırmasıdır:** [Sentinel-1](SENTINEL-1_SERISI.md) ve
[Sentinel-2](SENTINEL-2_SERISI.md). Bu dosya önceki C ayrıntılarını korur.
Kaynak kesiti 29 Ekim 2026; sayısal doğrulamanın gerçek kapsamı [HESAPLAR.md](HESAPLAR.md).

**Kanıt sınıfları:** **[D]** Sentinel-1C'ye özgü, birincil kaynakta doğrudan yazılı ·
**[M]** Sentinel-1 ailesi / A-B mirası · **[H]** bu çalışmanın hesabı
([HESAPLAR.md](HESAPLAR.md)) · **[S]** açıkça işaretlenmiş varsayım · **[Y]** kamuya açık
kaynakta bulunamadı.

## A. Genel görev bilgisi

### A1. Amaç, operatör, üretici

- **Amaç [D]:** Copernicus programı için gece/gündüz ve bulut koşulundan bağımsız
  C-bant (5,405 GHz) sentetik açıklıklı radar (SAR) görüntüleme: deniz buzu, deniz
  gözetimi, arazi deformasyonu (InSAR), afet ve tarım. Sentinel-1C ve 1D ayrıca **AIS**
  (gemi otomatik tanımlama sistemi) alıcısı taşır. [S1-ESA], [S1-VV25]
- **Operatör / müşteri [D]:** Operatör ESA; son müşteri Avrupa Komisyonu ("End-customer").
  Uçuş kontrolü ESA ESOC (Darmstadt). [S1-VV25], [S1-ESOC-MEDIA]
- **Üretici [D]:** Uydu ana yüklenicisi **Thales Alenia Space**; SAR cihazı **Airbus
  Defence and Space**. [S1-ESA]
- **Yörünge planı [M]:** Tek uydu 12 günde 175 tur (12 günlük tekrar); 180° fazlı iki uydu
  ile 6 günlük tekrar. C/D çiftinin 6 günlük deseni A/B'ye göre bir gün kaymıştır. [S1-WIKI]

### A2. Fotoğraf ve fiziksel bilgiler

![Sentinel-1C gerçek ESA fotoğrafı](figures/S1C_ESA_Photo.jpg)

Fotoğraf kaynağı, kredi ve lisans: [kaynak dizini](KAYNAKLAR.md).
Aynı fotoğraf seri raporunda açık araç kimliğiyle yer alır.

| Kalem | Değer | Kaynak ve kapsam |
|---|---|---|
| Fırlatma kütlesi | **2186 kg** | [D] VV25 fırlatma kiti s.3. Kitin resmî Arianespace kopyası bulunamadı; yalnızca üçüncü taraf aynası var (bkz. KAYNAKLAR) |
| ESA aile kütlesi | 2185 kg | [M] ESA Facts and figures |
| Tasarım kütlesi (A/B) | 2300 kg (130 kg yakıt dahil) | [M] ESA SP-1322/1 Tablo 4.1; SentiWiki aynı |
| Yüklenen itici | **154 kg hidrazin** | [D] ESA, 18.11.2024 ("fully loaded with 154 kg of fuel"); "fuelling the satellite with hydrazine" (ESA görsel sayfası) |
| Boyut | 21 × 2,5 × 4 m (açık), 2 × 10 m güneş paneli, 12 m radar anteni | [M] ESA Facts and figures |
| SAR anteni | 12,3 × 0,821 m; cihaz kütlesi 945 kg | [M] SP-1322/1 Tablo 5.2; SentiWiki |
| Güneş dizisi | 5900 W (ömür sonu, ortalama) | [M] SP-1322/1 Tablo 4.1; SentiWiki |
| Batarya | 360 Ah (A/B: 324 Ah) | [M] SentiWiki (C/D değeri A/B ile karşıtlıktan anlaşılır) |
| Operasyon yörüngesi | 693 km, güneş eşzamanlı, 98,18°, LTAN 18:00 | [M] ESA; SentiWiki |
| Enjeksiyon hedefi | ~700 km, 98,19°, LTAN 18:00 | [D] VV25 kiti s.2 (hedef; gerçekleşen yörünge elemanı değil) |
| Tasarım ömrü | "7,25 yıl için tasarlandı, 10 yıla kadar" | [D] VV25 kiti s.3 |
| Aile ömrü | "7 yıl (12 yıllık sarf malzemesi)"; ESA: "en az yedi yıl" | [M] SP-1322/1 Tablo 4.1; ESA Facts and figures |

### A3. Fırlatma ve tüm görev zinciri

| Olay | Bilgi | Kaynak |
|---|---|---|
| Kalkış | **5 Aralık 2024, 21:20 UTC** (kitin nominal saniyesi 21:20:33), Vega-C **VV25**, Avrupa Uzay Üssü (Kourou) | [D] VV25 kiti s.2; Arianespace "9:20 p.m. UTC" |
| Fırlatıcı | P120C (katı) → Zefiro-40 (katı) → Zefiro-9 (katı) → AVUM+ (sıvı, yeniden ateşlenebilir) | [D] VV25 kiti s.4–5 |
| Adaptör | "Vampire" yük adaptörü | [D] VV25 kiti; ESA "Building Vega-C flight VV25" |
| Toplam yük | 2286 kg (uydu 2186 kg; fark adaptör/ayırma donanımı) | [D] VV25 kiti s.2 |
| Planlanan olaylar | P120C ayrılması T+2:22, Z40 T+4:32, kaporta T+5:04, Z9 T+7:08; AVUM+ ateşlemeleri T+8 dk, T+53 dk, T+1 sa 40 dk; **uydu ayrılması T+1 sa 43 dk**; AVUM+ 4. ateşleme T+1 sa 49 dk | [D] VV25 kiti s.4–5 (nominal plan) |
| Gerçekleşen ayrılma | ESA: ~**6 Aralık 00:04 CET** (5 Aralık 23:04 UTC); ESA iletişimi 00:12 CET | [D] ESA "Double win for Europe" |
| Üst kademe bertarafı | ESA (fırlatma öncesi): ayrılmadan sonra son ateşleme üst kademeyi yörüngeden çıkaracak | [D] plan; uçuş sonrası teyit bulunamadı [Y] |

**Sistem sınırı:** Fırlatıcının itki ve kademe ΔV'si, uydunun 154 kg hidrazin bütçesine
dahil değildir. Uydunun ΔV bütçesi ayrılmadan sonra başlar.

### A4. Operasyonel durum (29 Ekim 2026) ve genel değerlendirme

| Tarih | Olay | Kaynak |
|---|---|---|
| 8 Aralık 2024 | İlk radar görüntüsü, kalkıştan 56 sa 23 dk sonra (Svalbard) | [D] ESA, 10.12.2024 |
| Ocak 2025 | Görev yörünge düğümüne ulaşmak için bir dizi manevra; SAR devreye alma Ocak ortası | [D] CDSE, 30.01.2025 |
| 26 Mart 2025 | Kullanıcılara veri açılışı | [D] CDSE |
| 19 Mayıs 2025 | Ürünler tam kalibre ("fully calibrated") | [D] Copernicus, 21.05.2025 |
| 18–19 Ocak 2026 | Anomali; 19 Ocak 19:40 UTC itibarıyla nominal veri | [D] CDSE |
| 3 Şubat 2026 | 2025'teki TRM arızaları sonrası anten örüntüsü ve kalibrasyon güncellemesi | [D] Copernicus, 05.02.2026 (arıza sayısı verilmiyor) |
| 9–23 Haziran 2026 | Sentinel-1C/1D yörünge yeniden yapılandırma manevraları | [D] CDSE, 28.05.2026 |
| 24 Haziran 2026 | 1C/1D son yapılandırma: 6 günlük nominal interferometrik tekrar | [D] Copernicus, 24.06.2026 |

- **Sentinel-1C:** operasyonel; ESA 30.06.2026: "Sentinel-1C and Sentinel-1D are performing
  very well". [S1A-END]
- **Sentinel-1A:** operasyonel görevi 29 Haziran 2026'da bitti (30 Haziran'da duyuruldu);
  bertaraf aşaması yaz döneminde başlayacaktı. [S1A-END]
- **Sentinel-1B:** görev sonu 3 Ağustos 2022; elektriksel pasivasyon 12 Eylül 2024;
  ~24 yılda atmosfere girmesi bekleniyor. [S1B-DISPOSAL]
- **Sentinel-1D:** 4 Kasım 2025'te fırlatıldı; 1 Mayıs 2026'dan beri tam operasyonel.
  26 Ekim 2026 16:20 UTC olayı (kurtarılamaz veri kaybı), 27 Ekim 17:32 UTC nominal;
  **28 Ekim 2026 tarihli yeni veri olayı** (ayrıntılı başlangıç/bitiş teyit edilmedi) (Operasyon Panosu,
  29.10.2026). Bu durum değişkendir; rapor teslimden önce panodan yeniden kontrol edilmeli.
  [S1-ESA], [S1D-LIVE], [S1-DASH]

**Genel değerlendirme:** C-bant SAR, bulut ve aydınlatmadan bağımsız gözlem ve InSAR için
faz kararlılığı sağlar. Bedeli yüksek elektrik yükü (SAR DC ~3,9 kW), yoğun veri akışı ve
çok sıkı yörünge kontrolüdür (interferometri için dar yörünge tüpü). İtki sistemi bu sıkı
yer izi kontrolünün ve ömür sonu bertarafının aracıdır.

## B. Bütçeler

Sentinel-1C'ye özgü, as-built alt sistem bütçeleri kamuya açık değildir. Aşağıda aile
(SP-1322/1, SentiWiki) değerleri **[M]** olarak, 1C'ye özgü olanlar **[D]** olarak verilir.

### B1. Kütle bütçesi

| Kalem | kg | Sınıf |
|---|---:|---|
| Fırlatma kütlesi (1C) | 2186 | [D] |
| Hidrazin (1C, yüklenen) | 154 | [D] |
| Fırlatma − hidrazin | 2032 | [H] yalnızca fark; doğrulanmış kuru kütle değil (basınçlandırıcı gaz dahil olabilir) |
| C-SAR cihazı (anten dahil) | 945 | [M] |
| İtici tankı (PEPT-590GB, birim) | 17,5 | [M] Rafael kataloğu s.41; uyduda tank sayısı [Y] |
| 1 N itici (birim, 1000 mm kabloyla) | ≤ 0,310 | [M] Rafael kataloğu s.46 |
| Diğer platform + marj | — | [Y] |

### B2. Güç bütçesi

| Kalem | Değer | Sınıf |
|---|---|---|
| Güneş dizisi ortalama gücü | 5900 W (ömür sonu) | [M] |
| SAR DC gücü | 3870 W (IW modu, tek veya çift polarizasyon) | [M] SP-1322/1 Tablo 5.2 |
| SAR RF tepe gücü | 4368 W (TRM çıkışlarının toplamı) | [M] SP-1322/1 Tablo 5.2 |
| SAR çalışma süresi | SM/IW/EW modlarında tur başına en çok 30 dk (A/B tasarımı: 25 dk). Mayıs 2026'ya kadar tutulma mevsiminde 25 dk ile sınırlandı; model iyileştirmeleriyle artık tüm yıl 30 dk | [M] SentiWiki |
| Valf gücü (1 N itici, birim) | 9,2 W @ 28 Vdc | [M] Rafael s.46 |
| Katalizör ısıtıcısı (birim) | iki eleman, eleman başına 257 Ω (28 V'ta ~3,05 W/eleman) | [M] Rafael s.46; [H] |

**Not (V2 düzeltmesi):** V2'deki "4075 W kritik RF/DC çelişkisi" kanıtlanamadı. SentiWiki
RF tepe satırında "4.368 kW, − 4.075 kW (IW, dual polarisations)" yazar; DC satırı 3870 W'tır.
4075 W'ı "DC güç" olarak veren bir kaynak bulunamadı (atıf yapılan makaleye
erişilemedi). Güç bütçesinde 3870 W DC kullanılır. 5900 − 3870 = 2030 W yalnızca bir farktır;
platform tüketimi, şarj ve kayıplar bilinmeden **marj değildir**. Tutulma koşulu için bkz.
[HESAPLAR.md](HESAPLAR.md) §7.

### B3. Yönelim (pointing) bütçesi [M]

- Donanım: 4 reaksiyon tekeri, 3 manyetik tork çubuğu, yıldız izleyiciler, GNSS. [S1-WIKI]
- Yönelim bilgisi 0,004°'den iyi (~70 µrad); yönelim doğruluğu ~0,01° (~175 µrad) eksen
  başına. [S1-WIKI], [S1-SP1322]
- Yörünge tüpü: SentiWiki'de "120 m çaplı (RMS)" tüp; aynı sayfada "60 m (RMS) yarıçap" ve
  tabloda "±120 m" ifadeleri de geçer (kaynak içi tutarsızlık). Uçuş dinamiği bu isteri
  60 m'lik mutlak yer izi ölü bandına dönüştürür. [S1-WIKI], [S1-ISSFD2012]
- 1C'ye özgü hata bütçesi (kovaryans, termoelastik pay) [Y].

### B4. Haberleşme linki bütçesi [M]

- X-bant: iki kanal × 260 Mbit/s; EDRS optik aktarma 520 Mbit/s. S-bant TT&C: 64 kbit/s
  yukarı, 128 kbit/s veya 2 Mbit/s aşağı. TT&C istasyonu Kiruna; SAR veri istasyonları
  Svalbard, Matera, Maspalomas, Inuvik. [S1-WIKI], [S1-SP1322], [S1-ESA]
- EIRP, G/T, kodlama ve gerekli Eb/N0 yayımlanmadığı için **gerçek link marjı
  hesaplanamaz** [Y]. [HESAPLAR.md](HESAPLAR.md) §8'deki örnek tamamen varsayımsaldır [S].

### B5. Veri bütçesi [M]

- Kayıt kapasitesi 1410 Gbit (ömür sonu) = 176,25 GB; indirme 2 × 260 Mbit/s. [S1-SP1322]
- 10 dakikalık bir geçişte iki kanalla ideal aktarım 39,0 GB [H]. Günlük üretim; mod,
  sıkıştırma ve görev planına bağlıdır, yayımlanmış günlük değer yoktur [Y].

### B6. Isıl bütçe

- Isıl dengenin ana kalemleri: SAR'ın ~3,9 kW DC tüketiminin ısıya dönüşmesi, güneş,
  albedo ve Dünya IR'si, radyatörler, MLI ve ısıtıcılar.
- **[D] 2026 termal korelasyonu:** Thales Alenia Space, 8 Ekim 2026 ESA termal mühendislik
  çalıştayında 1C'nin yaz gündönümü uçuş telemetrisiyle termal matematik modelini
  karşılaştırdı. %98 güven düzeyinde hesaplanan sıcaklık tahmin belirsizliği, korelasyonlu
  modelde korelasyonsuz modele göre **13 °C'den 7 °C'ye** (%46) indi. Bu bir sıcaklık
  sınırı veya marjı değildir. [S1-TMM-2026]
- Bileşen sıcaklık limitleri ve ısıtıcı görev çevrimi [Y]; ısıl marj hesaplanamaz.

### B7. ΔV bütçesi

Gerçek, kalem kalem 1C ΔV tahsisi yayımlanmamıştır [Y]. Kaynaklı referanslar:

| Kalem | m/s | Kaynak |
|---|---:|---|
| Fırlatma dağılımı düzeltmesi (S1B planlaması): eğim | 3,93 / 7,86 / 11,80 (1σ/2σ/3σ) | [M] ISSFD 2017 Tablo 4 |
| Dışmerkezlik | 1,78 / 3,60 / 5,43 (1σ/2σ/3σ) | [M] ISSFD 2017 §2.6 |
| Yörünge yuvasına hızlı yerleşme | 3,6–8,1 (nominal 5,3) | [M] ISSFD 2017 Tablo 3 |
| S1B gerçek enjeksiyon hatası düzeltmesi | yarı büyük eksen 1,0; dışmerkezlik 0,94; eğim 0,6 | [M] ISSFD 2017 §4.1.1 |
| 1C yörünge düğümüne geçiş (Ocak 2025), yeniden yapılandırma (Haziran 2026) | yapıldı; büyüklük yayımlanmadı | [D] / [Y] |
| Ömür sonu: 693 km → 530 km dairesel (1B'nin bertaraf hedefi) | **88,1** | [H] Hohmann; hedef [M] CGMS-52 |

Bu kalemler ile tasarım ömrü varsayımları birleştirilerek kurulan senaryo bütçesi
[HESAPLAR.md](HESAPLAR.md) §7'dedir; **ESA tahsisi değildir**.

### B8. Faydalı yük ve alt sistemler

- **C-SAR [M]:** 12,3 × 0,821 m aktif faz dizili anten; modlar SM/IW/EW/WV; RF tepe
  4368 W; DC 3870 W; kütle 945 kg; IW modunda 250 km şerit. [S1-SP1322]
- **AIS [D]:** gemi tanımlama alıcısı (yalnızca C/D). [S1-ESA]
- **Platform [M]:** PRIMA ailesi; itki modülü ayrı bir yapı olarak tüm itki donanımını
  taşır (SP-1322/1). Ayrıntılar B1–B6.

## C. CONOPS

| Evre | Zaman | Ana faaliyetler | Uydu ΔV |
|---|---|---|---|
| Yer hazırlığı | Kasım 2024 | Hidrazin dolumu (154 kg), "kritik ve tehlikeli" işlem | 0 |
| Fırlatma | 5 Aralık 2024 ~21:20 UTC → ayrılma ~23:04 UTC | Vega-C kademeleri yörüngeyi sağlar | Yalnızca fırlatıcı |
| LEOP | ilk günler | İlk temas (00:12 CET), panel açılımı, güvenli yönelim, yörünge belirleme | itici kalibrasyonu [Y] |
| Devreye alma | Aralık 2024 – Mayıs 2025 | İlk görüntü 8 Aralık; Ocak 2025 düğüm manevraları; kalibrasyon | düğüm manevraları (büyüklük [Y]) |
| Rutin | Mayıs 2025 → | IW/EW/SM/WV gözlemi, yer izi bakımı (60 m ölü bant), çarpışmadan kaçınma | düzenli düzlem içi + düzlem dışı manevralar |
| Takımyıldızı | Haziran 2026 | 1D ile 6 günlük tekrar için yeniden yapılandırma | faz manevraları (büyüklük [Y]) |
| Ömür sonu | gelecekte | Bertaraf manevraları, kalan yakıtın tüketilmesi, elektriksel pasivasyon | 1B örneği: 25 yıl yeniden giriş noktasına iniş |

**1B'nin gerçek bertarafı [M]:** Şubat 2023 RCT1 ile başlandı; plume kısıtları
nedeniyle Nisan 2023'te pitch dönüşü ve −Z iticileriyle yeni manevralar uygulandı.
SpaceOps #226 s.20'ye göre manevralar 23 Nisan 2024'te durduruldu; son yörünge yaklaşık
580×610 km, kalan hidrazin yaklaşık 70 kg, balistik yeniden giriş öngörüsü 24 yıldı.
Eski 530 km/<12 yıl hedefi ve ~88 kg ara tahmin, gerçekleşen sonuç değildir. Elektriksel
pasivasyon Eylül 2024'te yapıldı; yakıt kaldığı için tam pasivasyon sağlanamadı.
[S1-SPACEOPS226], [S1-SPACEOPS360], [S1B-DISPOSAL]

## D. İtki sistemi

### D1. Tip, mimari ve yerleşim

- **Tip [D]:** Hidrazin tek bileşenli (monopropellant), katalitik ayrışmalı kimyasal itki.
  1C'nin hidrazinle doldurulduğu ESA tarafından yazılıdır. [S1-FUEL-IMG]
- **İticiler [M]:** İki yedekli dal, **her dalda yedi adet 1 N itici** (toplam 14); ISSFD
  2017 bunları 3 çift 1 N yörünge kontrol iticisi (OCT) ve 4 çift yönelim kontrol iticisi
  (RCT) olarak tanımlar. OCT1/2/3 sırasıyla +X, −X ve −Y yüzeylerindedir. [S1-SPACEOPS360],
  [S1-ISSFD]
- **Tedarikçi [D]:** Rafael kataloğu: "28 units operational in-orbit (ST1 A & B), 28 units
  (ST1 C & D) were delivered to the customer". Uydu başına 14 sayısı ile tutarlıdır, ama
  1C'ye takılan adet doğrudan yazılı değildir. [S1-RAFAEL] s.10
- **Tank [M]:** Rafael **PEPT-590GB** (Sentinel-1 A–D için listelenmiş): Ø587 × 1170 mm,
  204 L, **4:1 blow-down oranında 154 kg nominal kapasite**, çalışma basıncı 27,4 → 5,5 bar,
  basınçlandırıcı GN₂ veya GHe, Ti-6Al-4V gövde, EPN-40 diyafram, 17,5 kg, uçuş ömrü 15 yıl.
  1C'ye yüklenen 154 kg bu tankın nominal kapasitesine eşittir; uyduda **tek tank**
  olduğu çıkarımdır [H], kaynakta yazılı değildir. [S1-RAFAEL] s.9, s.41
- **Basınçlandırma:** blow-down çalışma güçlü biçimde işaret edilir (4:1 blow-down tank,
  5,5 bar'a inen besleme aralığı, yakıt tahmininin tank basıncından yapılması). Kullanılan
  gaz (N₂/He) 1C için yayımlanmadı [Y]. [S1-SPACEOPS360]
- **Yerleşim [M]:** SP-1322/1, tüm itki donanımının ayrı bir itki modülünde olduğunu yazar.

### D2. İtici özellikleri (Rafael 1 N, katalog) [M]

| Özellik | Tip 1 | Tip 2 |
|---|---:|---:|
| Besleme basıncı | 24,5 → 9 bar | 24,5 → 5,5 bar |
| Başlangıç itkisi | 1,3–1,4 N | 1,0–1,1 N |
| Ömür sonu itkisi | > 0,2 N | > 0,2 N |
| Özgül itki (kararlı) | > 214 s @ 22 bar; > 205 s @ 5,5 bar | aynı |
| Toplam verilen impuls | 100 000 N·s | 60 000 N·s |
| Toplam darbe sayısı | 100 000 | 58 000 |
| En küçük impuls biti | 0,008 N·s @ 5,5 bar | aynı |
| Valf gücü / kütle | 9,2 W @ 28 Vdc / ≤ 310 g | aynı |
| Lüle genleşme oranı / ömür | 130 / 15 yıl (depolama + uçuş) | aynı |

Kaynak: [S1-RAFAEL] basılı s.46 (PDF s.49). Katalog **Sentinel-1'de hangi tipin uçtuğunu belirtmez**.
Tip 1 için 5,5 bar çalışma onayı çıkarılamaz; ortak Isp satırı uçuş haritası değildir.

### D3. İşlev

İtki sistemi enjeksiyon hatalarının düzeltilmesi, yörünge tüpünün (yer izinin) korunması,
çarpışmadan kaçınma, son güvenli modda yönelim kontrolü ve ömür sonu bertarafı için
kullanılır. [S1-SP1322] s.16. Fırlatma öncesi (2012) uçuş dinamiği analizine göre düzlem
içi (IP) manevralar ağırlıklı olarak atmosferik sürüklenmenin yarı büyük eksende yarattığı
azalmayı telafi eder; en yüksek enlemdeki yer izi sapmasını kontrol eden düzlem dışı (OOP)
manevraların sıklığını Güneş, Ay ve katı gelgitler belirler ve Ay etkisi nedeniyle bu sıklık
haftada birden fazla olmalıdır. [S1-ISSFD2012] s.4–5

### D4. Tsiolkovsky, yakıt, kütle debisi ve yanma süresi

$\Delta V = I_{sp} g_0 \ln(m_0/m_f)$, $\dot m = F/(I_{sp} g_0)$, $t = \Delta m/\dot m$.

- 154 kg'ın tamamı ideal kullanılırsa (2186 → 2032 kg): Isp 205 s ile **146,9 m/s**,
  214 s ile **153,3 m/s**, 220 s ile **157,6 m/s** [H]. Isp aralığı Rafael katalog
  alt sınırlarıyla kıyaslanır; 220 s yalnız varsayımsal duyarlılık noktasıdır.
- 1 N, 220 s: $\dot m$ = 0,4635 g/s; 1 m/s için 1,013 kg ve 36,4 dakikalık net 1 N
  ateşleme [H].
- Eski B hedefiyle ideal kıyas (C için kabul edilmiş EOL değil), 693 → 530 km dairesel: 88,1 m/s ve 220 s'de 87,4 kg (205 s'de 93,7 kg) [H].

**Gerçek ateşleme sınırları [M]:** S1A/S1B'de plume çarpması OCT1/2 ateşlemelerini
300 s yerine **35 s** ile sınırladı ve ~%30 performans kaybına yol açtı. OCT3, ±90° dönüşle
**600 s'ye kadar** kullanılabildi; ardından ısı ve batarya için **6 saat** normal modda
kalmak gerekti. 600 s tek dal için planlanan ΔV 0,28 m/s idi. [S1-ISSFD] §1.3, §2.3, §4.4.

### D5. Anlık mı sürekli mi?

Manevralar sürekli (elektrikli) itki değil, **ayrık kimyasal ateşleme dizileridir**. Kısa
yer izi düzeltmeleri dakikalar mertebesindedir ve yörünge periyoduna (98,6 dk) göre kısa
olduğundan anlık (impulsive) modele yaklaşır. Ömür sonu bertarafı ise onlarca saatlik toplam
ateşleme gerektirir ve çok sayıda sonlu yanmaya bölünür (1B örneği). [H], [S1-SPACEOPS360]

### D6. Rezerv ve ek ömür

1C'nin güncel kalan yakıtı yayımlanmamıştır [Y]. Bu nedenle "taşınan − gereken" farkından
ek ömür yılı **verilemez**. Sınırlı bir değerlendirme [HESAPLAR.md](HESAPLAR.md) §7'dedir.
Kalan yakıt uçuşta tank basıncından tahmin edilir; S1B için belirsizlik yaklaşık 5 kg
olarak raporlanmıştır. [S1-SPACEOPS360]

### D7. Güvenlik ve regülasyon

- Hidrazin toksik ve tehlikelidir; ESA 1C dolumunu "kritik ve tehlikeli" işlem olarak
  tanımlar. Dolum, sızıntı kontrolü, basınçlı sistem testleri ve fırlatma sahası tehlike
  analizi gerekir. [S1-FUEL-IMG]
- Bertaraf: 1B, 25 yıllık yeniden giriş ölçütünü sağlayacak yörüngeye indirilip
  pasive edildi. ESA'nın Kasım 2023 tarihli 5 yıl politikası için görev baseline'ı, geçiş ve
  uzatma değerlendirmesi gerekir; otomatik muafiyet veya uyum atanamaz. [S1-SPACEOPS360], [DEBRIS-ESA]
- Zero Debris Charter hukuken bağlayıcı bir düzenleme değildir. [DEBRIS-CHARTER]

### D8. TRL

ESA tanımında TRL 9, gerçek görev operasyonunda kanıtlanmış teknolojidir. Rafael 1 N
iticiler ve PEPT-590GB tank, Sentinel-1A/B ile uçuş mirasına sahiptir (TRL 9 düzeyinde
uçmuş tasarım). 1C bileşenleri için uçuş öncesi resmî TRL değerlendirmesi yayımlanmamıştır
[Y]; "uçuş öncesi TRL 6 → sonrası TRL 9" gibi bir geçiş belge olmadan yazılamaz. AIS gibi
yeni yükler ayrı değerlendirilir. [TRL-ESA]

## Açık kalan bilgiler

Uyduda tank sayısı ve basınçlandırıcı gaz; 1C'de uçan itici tipi; as-built kütle, güç, ısıl,
link ve yönelim bütçeleri; gerçek manevra büyüklükleri ve kalan yakıt; resmî TRL kayıtları.

## Kaynaklar

Kısa kodlar ve tam bağlantılar: [KAYNAKLAR.md](KAYNAKLAR.md). Başlıcaları: ESA Sentinel-1
Facts and figures <https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures>;
ESA SP-1322/1 <https://esamultimedia.esa.int/multimedia/publications/SP-1322_1/offline/download.pdf>;
Rafael uzay itki kataloğu <https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf>;
ISSFD 2017 <https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf>;
SpaceOps 2025 #360 <https://publications.spaceops.org/2025/download_by_id.php?id=0360>.
