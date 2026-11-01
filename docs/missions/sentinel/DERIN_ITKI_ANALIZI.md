# Sentinel-1 ve Sentinel-2: motor, besleme ve ayrıntılı itki hesabı

**SASS / UZM 508E · V5 · 1 Kasım 2026.** Karşılaştırma ilk nesil A/B/C/D serilerini
kapsar. Next Generation araçları kapsam dışında. Bu ek, [S1 seri raporu](SENTINEL-1_SERISI.md),
[S2 seri raporu](SENTINEL-2_SERISI.md) ve [genel türetimleri](HESAPLAR.md) derinleştirir.
Girdi/sonuç dosyası: [DERIN_ITKI_SONUCLARI.json](DERIN_ITKI_SONUCLARI.json).

## 1. Verilerin kanıt düzeyi

Bu çalışmada **K** yayımlanmış görev verisi, **M** aile/katalog/tasarım mirası,
**H** türetim, **V** senaryo varsayımı ve **A** bulunamayan veri anlamına gelir.
Bir üretici kataloğunun qualification tablosu, belirli uçuş biriminin kabul raporu
değildir. Bir C aracının yakıt dolumunu başka araçların kalkış kütlesiyle eşleştirmiyoruz.
Yörünge tayininden kalibre edilen ΔV, tanktan tüketilen yakıtın doğrudan ölçümü değildir.

| Eşleşen hesap girdisi | Sentinel-1C | Sentinel-2C |
|---|---:|---:|
| Yayımlanan fırlatma kütlesi [K] | 2186 kg | 1143 kg |
| ESA dolum haberindeki hidrazin [K] | 154 kg | 133 kg |
| Fark [H], sertifikalı kuru kütle değil | 2032 kg | 1010 kg |
| Görev irtifası, aile [M] | yaklaşık 693 km | yaklaşık 786 km |
| Kullanılan performans varsayımı [V] | 205 s, ayrıca 214/220 s duyarlılığı | 220 s, ayrıca 200/223 s duyarlılığı |

Kaynaklar: S1C VV25 kiti s.3 (Arianespace yazarlı üçüncü taraf aynası), S2C resmî
VV24 kiti s.3 ve [ESA S1C dolumu](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff),
[ESA S2C dolumu](https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/).
Tam URL/sürüm defteri [kaynak dizininde](KAYNAKLAR.md) bulunur.

### 1.1 A/B/C/D’ye aktarılabilecek ve aktarılamayacak bilgi

Rafael 2024 kataloğu basılı s.43, 1 N iticinin Sentinel-1 **A, B, C ve D** için
seçildiğini açıkça söyler. Dolayısıyla üretici/aile bağı doğrudandır. Bununla birlikte
Tip 1/Tip 2 konfigürasyon dağılımı ve uçuşa özgü seri numarası, thrust haritası, son
kalifikasyon sapmaları açık değildir. SpaceOps #360’ın iki dal × yedi itici tarifi,
S1B operasyon kaydıdır. Ailenin toplam 14 donanımı aynı yönde 14 N olarak toplanamaz.

Sentinel-2 için Airbus yazarlarının 2014 GNC makalesi Tablo 5'te **CHT1N-6 / EADS ST**
yazar. Aynı metin iki soğuk yedekli dal × dört eğik itici, blow-down tank ve modül
entegrasyonunu anlatır (s.5,9,12). ArianeGroup broşürünün s.6 miras listesinde 2A/2B
vardır. Bu kanıt, her C/D iticisinin aynı alt sürümde olduğunu tek başına kanıtlamaz.
S2C için hidrazin ve 21 bar helyum dolumu doğrudan ESA haberiyle doğrulanır.

## 2. Kullanılan motorların ayrıntılı karşılaştırması

| Özellik [M] | Rafael 1 N ailesi, Sentinel-1 bağı | ArianeGroup 1 N ailesi, S2 A/B tasarım bağı |
|---|---|---|
| Kimyasal tip | Katalitik hidrazin monopropellant | Katalitik hidrazin monopropellant |
| Valf | Solenoid, çift bobin/çift seat, yayla NC | Tek gövdede seri iki NC monostabil valf |
| Giriş basıncı | Tip 1: 24,5–9 bar abs; Tip 2: 24,5–5,5 bar abs | 5,5–22 bar |
| Kuvvet | BOL Tip 1: 1,3–1,4 N; Tip 2: 1,0–1,1 N; EOL >0,2 N | 0,320–1,1 N; nominal 1 N |
| Özgül impuls | >214 s @22 bar; >205 s @5,5 bar, SSF | nominal 220 s; aralık 200–223 s |
| Debi | Eşleşmiş F/Isp’den hesaplanır; tam harita açık değil | nominal 0,44 g/s; aralık 0,142–0,447 g/s |
| Minimum impuls biti | 0,008 N·s @5,5 bar ve 0,02 s/1 s duty | 0,010–0,043 N·s |
| Nozul alan oranı | 130 | 80 |
| Birim kütle | ≤310 g, 1000 mm kabloyla | 290 g, valflerle |
| Toplam impuls | Tip 1: 100000; Tip 2: 60000 N·s | 135000 N·s |
| Darbe/çevrim | Tip 1: 100000; Tip 2: 58000 | 59000 çevrim |
| Diğer ömür verileri | 15 yıl depolama+uçuş, tabloda | 67 kg throughput, 12 saat tek yanma, 50 saat toplam, <10°C’de 10 soğuk başlangıç |
| Elektrik | FCV 9,2 W @28 V; 23–36 V | Ürün metninde yedekli katalizör ısıtıcısı; bu incelemede valf/ısıtıcı W açık |
| Isıl/filtre | çalışma 5–90°C; 15 µm mutlak filtre | Isı bariyeri ve yalıtım; ayrıntılı uçuş limitleri açık |

Birincil kaynaklar: [Rafael teknik katalog](https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf),
basılı s.43–46; [ArianeGroup hidrazin broşürü](https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf),
s.5–6; [üretici 1 N sayfası](https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html).
İnternet sayfası minimum impuls bitinin birimini N yazıyor, PDF s.6 **N·s** diyor.
İmpulsun fiziksel boyutuyla tutarlı PDF birimini kullanıyoruz. Broşürde 2B’nin 2016
miras tarihi eski plandır; gerçek 2017 fırlatma tarihi yerine geçirilmez.

### 2.1 Katalog nominal değerleri neden aynı nokta sayılamaz?

`F = ṁ g0 Isp`, `g0=9,80665 m/s²`.
1 N ve 220 s eşleşmesi `ṁ=0,463507 g/s` gerektirir. Buna karşılık katalogdaki
0,44 g/s ve 220 s çarpımı `F=0,949284 N` verir. Fark yaklaşık %5,07’dir.
Bu bir fizik ihlali ilanı değildir: nominal isimlendirme, yuvarlama ve farklı
referans koşulları aynı tabloda olabilir. Bu nedenle bir yanda 1 N, diğer yanda
0,44 g/s kullanıp aynı yanmada farklı toplam impulslar üretmek kabul edilmez.

Hesaplarda ya F/Isp ikilisi seçilir ve debi türetilir ya da ölçülmüş eşleşen
F/ṁ/Isp haritası kullanılır. Basınç, yatak sıcaklığı, pulse width ve duty haritası
bulunmadan uç noktalar arasında doğrusal thrust haritası oluşturulmadı.

### 2.2 Geometri ve entegrasyon

ArianeGroup'un açık boyut çizimindeki referans uzunlukları valf tarafında 97 mm,
montaj yüzeyinden sıcak uca 75 mm, ısıtıcı bölümünde 42 mm’dir. Bunlar resimdeki
tanımlı referanslar arasındadır, bütün uçuş kablo zarfı değildir. Kesit/arayüz
görselinde Ø30 mm, Ø41 mm ve Ø4,3 mm delik bilgileri bulunur. Genel ürün çizimi,
S2C’ye ait imzalı ICD kabul edilmez. Çizimler yerel kaynak paketindedir.

S1 tankı platform merkezinde, iticiler yönelim/yörünge işlevlerine göre farklı
yönlerdedir. S2 tank/boru/itici grubu fırlatıcı arayüz halkası içinden entegre
edilebilen bağımsız modüldür. Hat uzunlukları ve moment kolları olmadığı için tam
yerleşim, hat basınç kaybı ve tork sayıları uydurulmaz.

## 3. Katalitik süreç ve motor hesabının fiziksel sınırı

Hidrazinin temsilî ayrışma yolu `3 N2H4 → 4 NH3 + N2` ve ardından
`4 NH3 → 2 N2 + 6 H2` şeklinde yazılabilir. Amonyağın ayrışan oranı sıcaklık ve
katalizör koşuluna bağlıdır. Sabit ürün molekül ağırlığı ve sabit γ her durumda
doğru değildir. Bu nedenle katalog Isp’sini doğruladığımızı iddia eden hayalî
oda sıcaklığı/oda basıncı/nozul boğaz alanı hesabı yapılmadı.

Genel eşitlikler `F=ṁ ve+(pe−pa)Ae`, `c*=pc At/ṁ`, `CF=F/(pc At)` ve
`Isp=CF c*/g0` olur. Bilinen F/Isp ve nozul alan oranı, tek başına pc, At,
katalizör sıcaklığı, verim ve genleşme durumunu tekil olarak belirlemez.
Nozul oranı 80 ile 130’u yalnız karşılaştırıp Isp üstünlüğü ilan edemeyiz.
Katalizör başlangıcı, ayrışma, ısı transferi ve pulse kayıpları da etkilidir.

Genel hidrazin özellikleri için [NASA/TM-2013-216560](https://ntrs.nasa.gov/api/citations/20140009582/downloads/20140009582.pdf),
Tablo 1 (basılı s.5 / PDF s.8): 300 K’de yoğunluk **1004 kg/m³**. Bu yoğunluk
varsayımı, gerçek tankın bütün sıcaklık/basınç geçmişini temsil etmez.

## 4. Vektörel kuvvet, kütle debisi ve tork

Her iticinin yön birimi ui, konumu ri, açık kalma oranı di olsun:

```text
Fnet = Σ di Fi ui
τnet = Σ (ri − rCG) × (di Fi ui)
ṁ = Σ di Fi / (g0 Ispi)
ce,etkin = (Fnet · hedef_yön) / ṁ
```

Eşit Isp, simetrik eğim ve sıfır yanal bileşen için `η=cosθ` ve
`ce,etkin=g0 Isp cosθ`. Aynı duty bütün motorlarda yalnız süreyi uzatır.
Eşitsiz off-modulation ise kuvvet yönünü ve torku da değiştirir, tek ortak duty
modeline indirgenemez. Yakıt, bütün iticilerin skaler debisinden gider.

| 1 m/s örneği [V/H] | S1C tek 1 N, Isp205 s | S2C dört 1 N, Isp220 s |
|---|---:|---:|
| Cant 0° yakıt | 1,08710 kg | 0,52967 kg |
| Cant 0° süre | 2185,46 s | 285,68 s |
| Cant 20° yakıt | 1,15684 kg | 0,56365 kg |
| Cant 20° süre | 2325,68 s | 304,01 s |
| Cant 30° yakıt | 1,25522 kg | 0,61158 kg |
| Cant 30° süre | 2523,45 s | 329,87 s |

Bu açılar uçuşta ölçülen yerleşim açıları değildir. Birim kuvvet 0,5 N olursa
aynı sabit Isp varsayımında yakıt değişmez, süre iki katına çıkar. Gerçek blow-down’da
Isp de değişebileceğinden bu sonuç koşulludur. S1’in bir saatten uzun tek motor
örneği, yayımlanmış plume/yanma kısıtına uyum göstermez. Operasyonel bir öneri değildir.

### 4.1 Küçük darbelerin çözünürlüğü

Küçük darbe için `δv≈Ibit/m`. S1C’de 0,008 N·s → **0,00366 mm/s**;
S2C’de 0,010–0,043 N·s → **0,00875–0,03762 mm/s** tek itici, başlangıç kütlesi
ve η=1 varsayımıyla. Dört iticiyi eşzamanlı ateşleme, cant, valf gecikmesi,
yatak başlangıcı ve kontrol örneklemesi uygulanabilir net çözünürlüğü değiştirir.
Katalog MIB değeri, uçuşta garanti edilen ΔV doğruluğu değildir.

## 5. Tsiolkovsky, değişen Isp ve olay defteri

Momentum dengesiyle `dv=−ce,etkin dm/m`. Sabit ce halinde:

```text
ΔV = g0 Isp η ln(mi/mf)
mp = mi [1 − exp(−ΔV/(g0 Isp η))]
mf = mi − mp
ṁ = Σ diFi/(g0Isp)
t = mp/ṁ
```

Küçük ΔV’de sayısal anlamlılık için kod `−expm1(−ΔV/ce)` kullanır. Motor kuvveti
ideal roket denkleminin ΔV kapasitesini doğrudan belirlemez, süreyi belirler.
Isp’nin kütleyle değişmesi halinde `ΔV=∫mf^mi ce(m)/m dm`; aritmetik ortalama Isp
genellikle doğru eşdeğer değildir. Ayrık, her biri sabit ce'li segmentler için
`ΔV=Σ ce,j ln(mj/mj+1)`. Güvenilir harita olmadan bu integral uçuş tahminiyle kapatılamaz.

S1C bütün 154 kg'ı η=1 ile kullanabilse 205/214/220 s varsayımlarında
146,86/153,31/157,61 m/s; S2C bütün 133 kg için 200/220/223 s ile
242,63/266,89/270,53 m/s verir. Bunlar ideal kapasite. Bertaraf rezervi,
kontrol tüketimi, kullanılamayan yakıt ve motor ömür tahsisi düşülmemiştir.

### 5.1 S2C yayımlanmış manevralarının rekonstrüksiyonu

[SpaceOps 2025 #617](https://publications.spaceops.org/2025/download_by_id.php?id=0617),
s.15–16 ve Tablo 3: dört büyük IP toplamı yaklaşık 6,33 m/s, 29 Ekim ince edinim
0,012 m/s ve iki ilk OOP toplamı yaklaşık 0,57 m/s. OOP metni planlama dili de
içerir. Bunlara Tablo 3'ün sekiz olayı eklenir. Toplam **yaklaşık 7,97 m/s**.
Aşağıdaki basamaklar yalnız tekrar üretim içindir, kaynak doğruluğu bu kadar yüksek değildir.

| Olay / tarih 2024 | Kalibre ΔV bileşeni [m/s] | Süre [s] | İdeal yakıt defteri |
|---|---|---:|---|
| Büyük IP/OOP ve 29 Ekim toplamı | 6,33 + 0,57 + 0,012 | yayımlanmamış | Sıralı kütle güncellemesi |
| 31 Ekim | 0,0146 | 5,3 | aynı |
| 6 Kasım | 0,9344 / 0,0537 | 314,1 | ortogonal norm 0,9359418 |
| 8 Kasım | 0,0025 | 0,9 | aynı |
| 23 Kasım ilk CAM | −0,0207 | 7,3 | büyüklük 0,0207 |
| 23 Kasım dönüş CAM | +0,0206 | 7,3 | büyüklük 0,0206 |
| 27 Kasım | 0,0309 | 11,0 | aynı |
| 29 Kasım | 0,0029 | 1,0 | aynı |
| 18 Aralık | 0,0283 | 10,0 | aynı |

6 Kasım tek yanmasında bileşenler normla birleştirilir. Ayrı zamanda gerçekleşen
iki zıt CAM’ın yakıt maliyeti ise **0,0413 m/s** olur. Net işaretli −0,0001 m/s'yi
yakıt hesabına koymak yaklaşık 413 kat küçük değer üretir. 31 Ekim ikinci kez eklenmez.

220 s, η=1 modeli **4,2138 kg** tüketim ve **128,7862 kg** kalan verir.
Bu ölçüm değildir. Ek olarak 1140/1143/1146 kg, 200/220/223 s ve 0/20/30°
ile 27 köşe/ızgara senaryosu **4,1463–5,3636 kg** verir. Bu aralık bir güven aralığı,
Monte Carlo olasılığı veya fiziksel bütün ihtimallerin sınırı değildir. Girdiler
ve örneklenen aralıklar varsayımdır. Eksik olaylar/safe mode tüketimi bu aralıkta yoktur.

### 5.2 Süreden tersine kuvvet ve tanımlanabilirlik

6 Kasım olayı öncesi model kütlesi 1139,336 kg; ideal tüketim 0,494155 kg.
`Feksen=mp g0Ispη/t` ile η=1 için **3,39421 N** ortalama eksensel kuvvet çıkar.
İlk yaklaşım `m ΔV/t` yakın sonuç verir, tam denklem kütle değişimini korur.

Dört eş motor için `di Fi cosθi` toplamı gözlenir. Yaklaşık 3,394 N;
dört 1 N ve θ≈32° ile de, motor başına yaklaşık 0,849 N ve θ≈0° ile de
elde edilebilir. Duty ayrıca bilinmeyendir. Tek bir ΔV/süre kaydı cant açısını,
motor basıncını ve valf görev oranını ayrı ayrı belirleyemez. Isp belirsizliği
tüketimi değiştirir, küçük ΔV’de çıkarılan kuvveti çok az etkiler.

## 6. CONOPS’tan ΔV bütçesi

Fırlatıcı alt/üst kademelerinin ürettiği Dünya yörüngesine erişim ΔV'si uydu
hidrazin bütçesine yazılmaz. Uyduya ayrılma anındaki enjeksiyon kovaryansı girdi olur.
S1C Vega-C/AVUM+, S2C Vega/AVUM; A/B’lerde sırasıyla Soyuz/Fregat ve Vega/AVUM,
S1D’de Ariane 6/Vinci zincirleri ayrı fiziksel sistemlerdir.

| Evre | Kullanılacak veri / denklem | Belirsizlik |
|---|---|---|
| İlk yönelim/güvenli mod | Valf zaman geçmişi, açısal impuls ve moment kolları | Net translasyon ΔV’sinden çıkarılamaz |
| Enjeksiyon düzeltme | Gerçek durum vektörü, hedef yörünge ve kovaryans | Kit hedefi ölçülen enjeksiyon değildir |
| Faz/yer izi edinimi | Gauss varyasyonları ve yörünge yayılımı | Zaman ve yörünge toleransları |
| Düzlem içi bakım | Sürüklenme/eksantriklik kontrolü | Yoğunluk, CdA/m, güneş etkinliği |
| Düzlem dışı bakım | Eğim ve yerel saat/J2 kontrolü | Manevra düğüm/faz seçimi |
| CAM ve geri kazanım | Her olayın impuls normu, eş/dönüş ateşlemeleri | Olay frekansı ve performans hatası |
| EOL | Bertaraf hedefi, sonlu yanma, basınç ve ömür | Ulaşılan yörünge/yeniden giriş modeli |

### 6.1 Eğim düzeltmesi

Saf düzlem değişimi `ΔV=2v sin(Δi/2)`. 0,01° düzeltme için S1’de **1,3104 m/s**,
S2’de **1,3019 m/s** gerekir. Bu, irtifa kaynaklı hız farkının küçük olduğunu,
düzlem kontrolü ihtiyacının yıllık drag'dan ayrı hesaplanması gerektiğini gösterir.
Bu 0,01° değeri yayımlanmış yıllık hata değildir. J2 düğüm kayması
`Ωdot=−(3/2)J2 n(RE/p)² cosi` ile yerel saat kontrolüne bağlanır. SSO uyumunu
her yıl 0,01° varsayarak hesaplamak yeterli değildir.

### 6.2 Drag duyarlılığı

Sabit yoğunluklu örnek `aD=½ρ Cd A vrel²/m` ve `ΔVyıl≈aD tyıl`.
Yerine dönmeyen atmosfer yaklaşımında vrel≈yörünge hızı, Cd=2,2, A=10 m²:

| Varsayılan yoğunluk [kg/m³] | S1C [m/s/yıl] | S2C [m/s/yıl] |
|---|---:|---:|
| 10⁻¹⁴ | 0,0895 | 0,1690 |
| 10⁻¹³ | 0,8951 | 1,6898 |
| 10⁻¹² | 8,9515 | 16,8977 |

A=20 m² olursa sonuç iki katıdır. İki irtifada aynı yoğunluk **yalnız duyarlılık
karşılaştırmasıdır**, gerçek atmosfer tahmini değildir. Gerçek modelde ortak güneş/
jeomanyetik sürücüler altında irtifaya ve zamana bağlı NRLMSISE/JB2008 gibi model,
atmosfer dönüşü, tutum ve etkin alan gerekir. Birbirinden bağımsız yoğunluk ve
alan uçlarını olasılıksal %95 marj diye sunamayız.

### 6.3 Bertaraf geometrisi

Vis-viva `v²=μ(2/r−1/a)`. Dairesel r0’dan rp×ra’ya ilk olarak r0×rp transferi,
ardından rp’de hedef elipse geçiş yapılır. İki retrograde impuls büyüklüğü toplanır.

| Geometri [H] | İdeal iki impuls |
|---|---:|
| S1 693→530 km dairesel, eski B hedefi | 88,0576 m/s |
| S1 693→580×610 km, B son durumuyla kıyas | 52,5822 m/s |
| S2 786→653 km dairesel, 2014 aile planı | 70,2157 m/s |

B’nin son durum geometrisi, B’de gerçekte harcanan ΔV’nin ölçümü değildir.
Manevralar arasında sürüklenme, başarısız/yarım yanmalar ve yakıt harcayan yönelim
işlemleri vardır. S2’nin yayımlanmış 73 m/s planı ideal geometriden %3,965 fazladır.
Kaynak kayıpları ayırmadığı için bu farkın tamamına sonlu yanma kaybı denmez.

## 7. Sonlu yanma: sayısal yörünge deneyi

İmpuls yaklaşımı yanma boyunca konumu sabit sayar. Gerçek 600–700 s yanmalar,
yaklaşık 100 dakikalık LEO turunda önemli açı tarar. Bu nedenle ek kod iki boyutlu
iki cisim hareketini RK4 ile, değişen kütle ve hızın tersine sürekli yönlenen itkiyle çözer:

```text
r′=v
v′=−μr/|r|³ − (F/m) v/|v|
m′=−F/(g0Isp)
W′=a_itki · v
H′=(r × a_itki)z
```

Senaryo [V]: 693 km dairesel, 2100 kg, Isp205 s, skaler ∫F/m dt=0,7 m/s.
2100 kg, S1B’nin ölçülen o günkü kütlesi değildir. 2,18 N,
SpaceOps #360 Tablo 3’teki dört kuvvetin toplamına eşit bir senaryo büyüklüğüdür.
Yerleşim/yön/eşitsiz momentler bu modelde yoktur.

| Sabit kuvvet [N] | Süre [s] | Yakıt [kg] | Son oskülatör rp×ra [km] |
|---:|---:|---:|---|
| İmpuls referansı | 0 | 0,73108 eşdeğer | 690,36354 × 693,00000 |
| 1,00 | 1469,74 | 0,73108 | 690,49328 × 692,87021 |
| 2,18 | 674,19 | 0,73108 | 690,39150 × 692,97203 |
| 4,00 | 367,44 | 0,73108 | 690,37188 × 692,99165 |

2,18 N yanması turun **%11,393**’ünü kaplar. Aynı skaler roket ΔV’si, impuls
modelindeki apsislerle aynı sonucu vermez. 2,18 N örneğinde perigee yaklaşık 28 m
daha yüksek, apogee yaklaşık 28 m daha düşüktür. Bunu doğrudan “%x yakıt kaybı”
olarak yorumlayamayız. Hedef yörünge için yanma merkezi, yön ve süre yeniden optimize
edilmelidir. Görülen değerler sayısal eğitim deneyidir, gerçek B manevra rekonstrüksiyonu değildir.

1 s ve 0,5 s adımlarıyla rp/ra farkı 10⁻⁷ km'nin altındadır. Türetilen son enerji
değişimi ile entegre özgül iş, açısal momentum değişimi ile entegre itki momenti ve
yakıt defteri ayrı kontrol edilir. Adım küçültme ve denge kapanışı, eksik J2/atmosfer/
attitude modelini doğrulamaz. İşlem uçuş propagatörü veya çarpışma taraması için kullanılamaz.

## 8. Tank ve blow-down: erişilebilir yakıt

Rafael PEPT-590GB [M], basılı s.41: 204 L, 154 kg nominal kapasite @4:1,
Ø587×1170 mm, 17,5 kg, Ti-6Al-4V gövde, EPN-40 diyafram. İşletme 27,4→5,5 bar;
proof 41,1, burst 54,8 bar. Bunlar görev dolum basıncı veya motor girişi değildir.
S1C’nin gerçek basınç/ısı geçmişi bulunamadı. S2 tank hacmi bu S1 tankından kopyalanmadı.

İdeal gaz ve sabit sıvı yoğunluğu varsayımıyla:

```text
Vg0 = Vtank − mload/ρ
P Vg^n = P0 Vg0^n
Pcut,tank = Pmin,inlet + ΔPloss
Vg,cut = min[Vtank, Vg0 (P0/Pcut,tank)^(1/n)]
mres = ρ (Vtank−Vg,cut)
maccessible = mload−mres
```

204 L, 154 kg ve ρ=1004 kg/m³ birleşimi yaklaşık **50,614 L** ilk gaz hacmi
verir. S1 C yükü+aile tankı eşleştirmesi bir duyarlılık modelidir.
Pmin=5,5 bar ile Tip 2 zarfı seçilmiştir, Tip 1’in 9 bar altına onay olduğu söylenmez.
Hat/filtre/valf kaybı toplam **0,5 bar sabit varsayım**, modelle hesaplanmış değer değildir.

| P0 [bar abs, V] | n [V] | Basınç kesiminde kalan [kg, H] | Erişilebilir [kg, H] |
|---:|---:|---:|---:|
| 22,0 | 1,0 | 18,491 | 135,509 |
| 22,0 | 1,2 | 54,769 | 99,231 |
| 22,0 | 5/3 | 94,010 | 59,990 |
| 24,5 | 1,0 | 0 | 154 |
| 24,5 | 1,2 | 40,689 | 113,311 |
| 24,5 | 5/3 | 86,618 | 67,382 |

n=1 izotermal; n=5/3 ideal tek atomlu helyumun ısı alışverişsiz genişleme sınırıdır.
Aralarında aylarca yeniden ısınma bulunan görev, bütün yakıt boyunca adyabatik
modelle temsil edilemez. Bu tablo **S1C erişilebilir yakıt tahmini veya güven aralığı değildir**.
27,4 bar tank üst sınırını 0,5 bar kayıpla kullanırsak ilk motor girişi 26,9 bar olur,
24,5 bar katalog üst sınırını aşar. Bu nedenle 27,4 barı eşleşmiş çalışma noktası gibi kullanmadık.
27,4/5,5 basınç oranı ile katalog 4:1 hacim blow-down ifadesi aynı oran değildir.

Gerçek modelde Z, gaz/sıvı sıcaklığı, buhar basıncı, diyafram hacmi, ısı geçişi,
gaz çözünmesi ve ölçüm kalibrasyonu gerekir. `Pgiriş=Ptank−ΔPhat−ΔPfiltre−ΔPvalf`.
Darcy–Weisbach `ΔP=(fL/D+ΣK)ρu²/2`; fakat L,D,pürüzlülük, K ve debi haritası
olmadan 0,5 bar bir hesap sonucu diye sunulamaz.

### 8.1 Basınçtan kalan yakıt ölçümü

İzotermal, gaz kütlesi sabit ve Z≈1 ise `mres=ρ(Vtank−mg Rg T/P)`.
Yaklaşık `δm≈ρVg(δP/P−δT/T−δmg/mg)` olur. Basınç, sıcaklık ve gaz kütlesi
birbiriyle ilişkiliyse bağımsız RSS kullanılamaz. SpaceOps #360 s.3'teki yaklaşık
**5 kg** kalan yakıt belirsizliği belirli S1B operasyon bağlamıdır, C için ±5 kg
kabul test sınırı değildir. Geçmiş valf sürelerini sabit debiyle çarpmak da blow-down
nedeniyle sistematik hata üretebilir (#226 s.14).

## 9. Motor ömrü ve yakıt rezervi aynı bütçe değildir

Bir iticiye tahsis edilen `Ij=∫Fjdt`, `mj=∫ṁjdt`, toplam açık süre, çevrim sayısı ve
soğuk başlangıçlar ayrı ayrı sınanır. Yedek dal, otomatik olarak iki kat ömür sağlamaz.
Bir ortak tank tüm yakıtı taşısa da her itici tüm yakıtı işleyemez.

| İdeal yük paylaşımı [V/H] | Motor başına impuls | Katalog eşiğine oran |
|---|---:|---:|
| S1C 154 kg, 205 s, tek Rafael Tip 2 | 309596 N·s | 5,160 × 60000 |
| Aynı yük 4 iticiye eşit | 77399 N·s | 1,290 × 60000 |
| S2C 133 kg, 220 s, tek AG 1 N | 286943 N·s | 2,126 × 135000 |
| Aynı yük 4 iticiye eşit | 71736 N·s | 0,531 × 135000 |

S1’de bu tablo bütün motorların Tip 2 olduğunu veya uçuşta limit aşıldığını
kanıtlamaz. Tip 1 için aynı oranlar 100000 N·s’ye göre farklıdır. S2’de 33,25 kg/
motor, 67 kg katalog throughput’un yaklaşık %49,6’sıdır, ancak çevrim/ısı/tek yanma
sınırı hâlâ ayrıca gerekir. Ayrıca 67 kg×220g0≈144549 N·s ile 135000 N·s aynı anda
tek bir sabit nokta değildir. Kalifikasyon satırları tek eşzamanlı çalışma köşesi gibi kullanılamaz.

2014 S2 bertaraf planı 116×600 s = **69600 s = 19,333 saat** toplam yanma
öngörür. Dört iticinin hepsi her yanmada açıksa bu süre motor başına da 19,333 saattir,
4'e bölünmez. 47 kg eşit paylaşılırsa 11,75 kg/itici ve 116 ateşleme başlangıcı olur.
Off-modulation valf çevrimlerini artırır. Tek yanma 600 s, katalog 12 saat sınırından
kısa olsa da gerçek uzay aracı plume/termal limiti daha kısıtlayıcı olabilir.

### 9.1 Rezervden koşullu ek yıl

`p=mmevcut−mEOL−memniyet−merișilemeyen` pozitifse ve mevcut toplam kütle m,
sabit yıllık ΔVyıl ve etkin egzoz hızı ce varsayılırsa:

`Nyakıt = ce ln[m/(m−p)]/ΔVyıl`.

Örnek [V]: m=2100 kg, tankta 70 kg, toplam ayrılan/tutulan 60 kg, rutin için
10 kg, Isp205 s, η=1, yılda 1 m/s → **9,596 ideal yıl**. Bu, S1B’nin görev
uzatma tahmini değildir. 70 kg tarihsel son durumla aynı seçilmiş olsa da 2100 kg,
60 kg tahsis ve yıllık 1 m/s bu örnekte varsayımdır. Gerçek ek yıl
`N≤min(Nyakıt, Nbasınç, Nmotor, Nçevrim, Ntermal, NEPS, Nyük, Ngüvenlioperasyon)`.
Kalan fiziksel yakıt, erişilebilir ve güvenle kullanılabilir yakıtla aynı değildir.

## 10. S1B: doğrulanmış operasyon dersleri

[SpaceOps #360](https://publications.spaceops.org/2025/download_by_id.php?id=0360)
Tablo 3'te simülatöre girilen dört kuvvet **0,57/0,56/0,47/0,58 N**'dur.
Toplam 2,18 N; maksimum–minimum farkının ortalamaya oranı **%20,18**.
Bu oran bir motorun bozulma yüzdesi değildir. Moment için montaj konumları ve
yönler gerekir, yalnız toplam kuvvet torku belirlemez.

2017 edinim kaydı, normal düzlem içi iticilerde plume kaynaklı süre kısıtını ve
yaklaşık %30 etkili performans kaybını anlatır. 2023 bertarafında −Z iticiler için
pitch dönüşüyle farklı çalışma modu kullanılır. Eşitsiz performans teker hızlarını
yükseltip safe mode doğurur. 2023’te geliştirilen Stop Burn FDIR, yalnız teker hızına
bakmakla kalmayıp tutum ve açısal hız eşiklerini de kapsar.

Makaledeki eşikler (Tablo 6) o konfigürasyon için ±250 rad/s teker hızı,
±0,01 rad tutum hatası ve ±0,001 rad/s hız hatasıdır. Çift motor arızası simülasyonu
yapılmadığı açıkça belirtilir. Bu eşikleri başka sistem için güvenlik kuralı diye kopyalamıyoruz.

[SpaceOps #226](https://publications.spaceops.org/2025/download_by_id.php?id=0226)
s.20 ve ESA kapanış haberi: manevra kampanyası Nisan 2024'te yaklaşık 580×610 km’de,
yaklaşık 70 kg yakıtla sonlandırıldı; yeniden giriş yaklaşık 24 yıllık model öngörüsüydü.
Elektriksel pasivasyon Eylül 2024'te ayrı gerçekleşti. Tankı boşaltacak özel valf
olmadığı ve iticilerle boşaltmanın öngörülemez yörünge değişimi yaratacağı için tam
yakıt/presürizasyon pasivasyonu mümkün olmadı. Bu durum yeni tasarımda EOL
pasivasyonunu ilk günden gereksinim yapmanın somut gerekçesidir.

## 11. Güç, ısıl, pointing ve diğer bütçelerle kapanış

Rafael FCV 9,2 W @28 V ve her 257 Ω eleman için `P=V²/R=3,0506 W`.
Üç valf + her birinde tek eleman açık örneği 36,7518 W. Bunu 14 iticiye sürekli
yaymak doğru değildir. Soğuk yedek dal ve ısıtıcı duty/ön ısıtma süresi ayrı tutulur.
AG 1 N için Rafael’in elektrik değerini kopyalamıyoruz.

Isıl bütçe, elektik ısıtmasının yanında katalizör/gövde iletimi, plume, gölge,
Güneş/albedo/Dünya IR ve radyasyon içerir. `Qin−Qout=dU/dt` ve `Prad=εσA(T⁴−Tenv⁴)`.
Yerleştirme/view factor/ısıl direnç ve yatak kabul limitleri olmadığından toplam
modül ısıl kapanışı gösterilemez. Yakıt donmasıyla kalan basınç ve valf sızıntısı
aynı risk değildir.

S2 2014 AOCS gereksinimleri: pointing ≤1200 µrad/eksen (%99,7), knowledge
≤10 µrad/eksen (%95,5), rate ≤20 µrad/s (%99,7). Bunların olasılık seviyeleri
farklıdır, tek 3σ RSS bütçesinde doğrudan toplanamaz. İtki darbesinin yön/tork
hatası, optik görüntüleme sırasında normal tekerli kontrol koşuluyla karıştırılmaz.

S1 SAR güç/mod profili, S2 MSI gündüz gözlemi, optik iletişim pointing’i ve
downlink fırsatları manevra zamanlamasını sınırlar. Önceki seri raporlarındaki
kütle, güç, link, veri ve ısıl bütçe açıkları korunur. Katalog W bilgisi, uçuş EPS
bütçesi değildir. ΔV’yi hesaplamak bütün sistem bütçelerini kapatmaz.

## 12. Güvenlik, regülasyon ve TRL

Hidrazin toksisitesi, katalitik ayrışma, yüksek basınç, valf sızıntısı, plume ve
pasivasyon riskleri tasarım doğrulama konusudur. Katalog proof/burst değerlerini
uçuş basıncıyla karıştırmıyoruz. Basınç altında dolum/ateşleme işlemi için bu rapor
bir uygulama prosedürü sunmaz.

Eski görevlerin 25 yıllık bertaraf tasarımı ve ESA’nın 2023'te yürürlüğe giren
5 yıllık yaklaşımı, farklı tarihli gereksinimlerdir. Projeye uygulanacak sürüm ve
kurumsal kapsam doğrulanmalıdır. [ESA 2023 uzay enkazı yaklaşımı](https://www.esa.int/Space_Safety/Space_Debris/ESA_leads_the_way_towards_a_Zero_Debris_future).
Güncel esaslar ve görev özelindeki uyum kanıtı olmadan “tam regülasyon uyumlu” denmez.

TRL birimi belirtilmelidir: motor, tank, entegre modül veya uzay aracı.
Uçuşta başarılı kullanılan donanım/konfigürasyon için TRL9 yorumu yapılabilir.
Önceki geliştirme evresinin resmî TRL’si açık raporda yoksa otomatik TRL6 atanmaz.
Sentinel mirası, kendi motor+besleme tasarımımıza TRL5 sağlamaz. İlgili ortamda
temsilî bileşen/breadboard testleri, ölçüm belirsizliği ve sistem arayüz kanıtları gerekir.

## 13. Güncel durum düzeltmesi ve seriler arası kararlar

V4’te 28 Ekim S1D olayının kapanışı açık bırakılmıştı. 1 Kasım taramasında
[30 Ekim CDSE bülteni](https://dataspace.copernicus.eu/news/2026-10-30-copernicus-sentinel-1d-data-unavailability-28-october-2026)
28–29 Ekim aralığını doğruluyor; [operasyon panosu](https://operations.dashboard.copernicus.eu/index)
29 Ekim 15:29 UTC’de nominale dönüş bildiriyor. Kaybolan veri geri gelmiş sayılmaz.
Bu olay için itki kaynaklı neden yayımlanmamış, öyle nitelendirilmedi.

S1A bilim görevi sonu ve S1B elektriksel arızası, hidrazin tükenmesiyle aynı olay
değildir. S2 A/B/C/D için donanım veri sürümleri ayrı tutulur; D’nin yerde depolama
durumu uçuşta TRL9 kanıtı sayılamaz. Her iki seride de yedeklilik, darbe tekrarlanabilirliği,
yatak/hat ısısı, nominal görev ile EOL geometrisi ve kalan basınç birlikte önemlidir.

Yeni tasarıma aktarılacak somut gereksinimler: görev bazlı ΔV ve manevra zamanı,
net kuvvet/tork tahsisi, erişilebilir yakıt, her iticinin impuls/çevrim/ısı defteri,
yakıt ölçüm belirsizliği ve EOL pasivasyonu. Görev hedefi ve kütle henüz seçilmediği
için bu incelemeden kendi sistemimize kesin motor veya tank boyutu seçilmedi.

## 14. Tekrar üretim ve bağımsız kontroller

```powershell
python scripts/sentinel_series_budget.py --output docs/missions/sentinel/SERI_SONUCLARI.json
python scripts/sentinel_propulsion_deep_dive.py --output docs/missions/sentinel/DERIN_ITKI_SONUCLARI.json
python -m unittest discover -s tests -p "test_sentinel*.py"
```

V5 ekinin kontrolleri: PVⁿ ve kütle kapanışı, duty/kuvvet/cant ayrımı,
enerji–iş ve açısal momentum–itki momenti kapanışı, adım yarılama, impuls sınırına
yaklaşım ve hatalı girdilerin reddi. Önceki 26 test ayrı korunur. Sayısal kontroller,
uçuş yakıt ölçümü veya motorun temsilî ortamda denenmesi yerine geçmez.

## 15. Kaynak taramasının kapsamı ve açık kalan belgeler

Temel kaynaklar: ESA görev/operasyon/dolum sayfaları, Arianespace araç bazlı kitler,
ESA SP-1322/1 ve /2, ISSFD 2012/2017, SpaceOps 2025 #226/#360/#617,
Airbus 2014 GNC yazar makalesi, Rafael 2024 teknik kataloğu ve ArianeGroup ürün
broşürü/çizimleri. NASA genel hidrazin özelliği ve itki durum raporu çapraz kontrol sağlar.
İndirilen dosyalar paket içinde URL, erişim zamanı ve SHA-256 ile kayıtlıdır.

2014 GNC yazar metni web üzerinden okunabildi. Resmî ESA GNC PDF erişimi giriş
gerektiriyor, ResearchGate otomatik PDF/HTML indirmesi başarısız. İndirilemeyen belge
pakette varmış gibi gösterilmedi. Yeni ESA LPS2025 C/D sunumunun indirmesi de başarısız
oldu; ayrıntılı motor verisi için ona dayanılmadı. Kamuya açık kataloglarda kullanım
hakları saklıdır; özellikle Rafael belgesindeki proprietary notu korunur, çizimleri
sunuma yeniden basılmadı. Kaynak kopyaları akademik inceleme arşividir.

**Kapanmayan uçuş verileri:** A/B/D gerçek yükleme/kuru kütle defterleri, C/D itici
alt sürümleri ve yerleşim açıları, her birimin basınç-sıcaklık-darbe performans haritası,
S2 tank hacmi/ICD, uçuş hat/filtre kayıpları, gerçek tank P/T ve kalan gaz kütlesi,
tam motor ömür defteri, C/D onaylı EOL tahsisleri ve güncel uzay aracı ısıl/güç bütçeleri.
Bu belgeler olmadan hesaplar ayrıntılı ve tekrar üretilebilir olabilir, ancak tam
uçuş tüketim/ömür tahminine dönüştürülemez.
