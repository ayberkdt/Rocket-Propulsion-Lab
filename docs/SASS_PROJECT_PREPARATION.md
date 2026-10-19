# UZM 508E / SASS: propulsion projesi hazırlığı

Araştırma tarihi: **18 Ekim 2026**. Bu belge bir proje hazırlığıdır; öğretim
üyesinin yayımladığı ödev yönergesi değildir. Sayısal örnekteki görev, kütle,
marjin, güç ve tercih ağırlıkları özellikle **sentetik tasarım varsayımı**
olarak işaretlenmiştir.

## 1. Doğrulanabilen ders bilgisi

- **UZM 508E — Spacecraft Architecture and Subsystems Design** dersinin
  koordinatörü Alim Rüstem Aslan; ders dili İngilizce. Resmî Ninova kaydının
  açık amaç/tanım alanları yalnızca “Sassd” içeriyor. Bu kayıt ayrıntılı bir
  haftalık içerik veya proje isterlerini kanıtlamıyor. [ITU-UZM508]
- ÖBS'nin erişilen ders grubunda ders 3 kredi, 7.5 AKTS, haftada 3 teorik
  saatle listeleniyor. Grup adı **Zorunlu Seçmeli Ders III**; bundan dersin
  bütün program öğrencilerine tek başına zorunlu olduğu sonucunu çıkarmıyoruz.
  Kullanıcının kendi ders/proje yükümlülüğü esas alındı. [ITU-UZM-GROUP]
- Aynı bölümün **lisans** UZB 419E ders kataloğu; görev/isterler, sistem
  bütçeleri, itki, alt sistem etkileşimi ve takım tasarım raporu içeriyor.
  Bu, kapsam tahmini için ilgili bağlamdır; **UZM 508E izlencesinin yerine
  geçirilmemiştir**. [ITU-UZB419]

## 2. Önceden hazırlanması mantıklı proje dalları — çıkarım

Aşağıdakiler ders adından, ilgili bölüm tasarım dersinden ve sistem tasarımı
iş akışından çıkarılan hazırlık seçenekleridir; kesin ödev beklentisi veya
öğretim üyesinin istediği teslimatlar değildir. [ITU-UZM508; ITU-UZB419;
NASA-SE-APP]

| Olası proje dalı | Hazırlanan analiz | Sonradan gereken görev verisi |
|---|---|---|
| Küçük uydu için itki alt sistemi seçimi | Δv kalemleri, yakıt/rezerv, kütle-güç-zaman bütçesi, kaynaklı seçenek karşılaştırması | Gerçek yörünge, ömür, manevra listesi, yük/kütle ve operasyon sınırları |
| Kimyasal ve elektrikli itki karşılaştırması | Ani Hohmann ve düşük itki spiral bütçeleri ayrı; mevcut elektrikli itki modeli ve EPS arayüzü | Eşleşmiş itki-Isp-güç haritası, güç geçmişi ve görev süresi |
| Basınç beslemeli/soğuk gaz sistemi | Tank hacmi, kullanılabilir gaz, basınçlandırıcı envanteri; mevcut dinamik besleme modellerine geçiş | Depolama koşuluna ait yoğunluk/gaz verisi, gerçek donanım ve sıcaklık |
| Hassas manevra / ADCS ile ortak itki | Mevcut pulse-mode, gimbal, motor yerleşimi, kuvvet-moment ve kütle türevleri | İtki darbesi ölçümleri, yerleşim, kütle merkezi ve kontrol isterleri |

Kaynak dayanağı: görev/alt sistem tasarım yaklaşımı için [ITU-UZB419;
NASA-SE-APP; NASA-SE-6]; manevra ayrımı için [MIT-522-L5; MIT-522-L6];
itki seçenekleri ve bileşen/sistem ayrımı için [NASA-SOA26]. Olası projelerin
öncelik sırası bu hazırlığın mühendislik tercihidir, istatistiksel olasılık değildir.

## 3. Yeni modüller ve kaynak bağı

Paket: `rocket_propulsion.studies.sass`. Çalışma kopyası:
`.worktrees/advanced-propulsion-sidera`. Ana klasörün daha eski `src/` kopyasına
dosyalar kopyalanmadı; geliştirme güncel dal üzerinde yapıldı.

| Modül / fonksiyon | Hesap / davranış | Kaynak ve kapsam |
|---|---|---|
| `references`, `Citation`, `CitedScalar` | Kaynak kimliği, alanın gerçek konumu, erişim tarihi, yayımlanmış/varsayım/türetilmiş ayrımı | [NASA-SE-APP] Ek D; [NASA-SE-6] 6.6/6.8. Veri şeması bu yazılımın tasarımıdır. |
| `delta_v_budget` | Her kaleme `count × nominal × (1+margin)` | [NASA-SE-6] tahsis/varsayım yönetimi. Marjin yüzdesi kaynak tarafından zorunlu kılınmaz; kullanıcı seçer. |
| `hohmann_transfer` | İki ayrı, işaretli ani teğetsel hız değişimi ve yarım elips süresi | [MIT-522-L6] Denklem 7; [MIT-522-L5] yörünge enerjisi. Yarıçaplar merkezden ölçülür. |
| `circular_spiral_delta_v` | Yakın dairesel, düşük teğetsel itki için hız harcaması yaklaşımı | [MIT-522-L6] Denklem 4–8. Kaçış, gölge, yönlendirme ve güç değişimi çözülmez. |
| `perigee_lowering_delta_v` | Dairesel yörüngeden yalnızca ilk geri yönlü yanma | [MIT-522-L5] vis-viva; [MIT-522-L6] ilk Hohmann impulsu. Yeniden giriş süresi veya mevzuata uygunluk hesaplanmaz. |
| `burn_duration_screen` | Yanma süresi/yörünge periyodu, kullanıcının seçtiği sınırla kıyas | [MIT-522-L6] ani/düşük itki model ayrımı; [NASA-SE-6] kriter seçimi. Evrensel bir eşik uydurulmaz. |
| `size_mission` | Farklı Isp ve yön kaybıyla sıralı yakıt, impuls, süre ve ıslak kütle | [NASA-ROCKET; NASA-ISP; BIPM-SI]. Aşağıdaki rezerv denklemi bu kaynaklardan cebirsel olarak türetilmiştir. |
| `bipropellant_split` | Toplam kütleyi açık O/F ile iki bileşene ayırma | [ARIANE-10N] O/F tanımı; kütle korunumu. Artık yakıtların dağılımı ayrıca varsayılmalıdır. |
| `liquid_tank_volume` | `V_l=m/rho`, `V_t=V_l/(1-u)` | [NASA-GAS] yoğunluk tanımı; [NASA-SP8112] ullage bağlamı. Yoğunluk kullanıcıdan kaynaklı gelir. |
| `ideal_cold_gas_storage` | İki basınç sınırı arasındaki ideal, izotermal kullanılabilir gaz | [NASA-GAS; NASA-SP8112]. Yüksek basınçta gerçek gaz geçerliliği doğrulanmış sayılmaz. |
| `ideal_regulated_pressurant` | Aynı sıcaklık varsayımıyla tüp-ullage envanter dengesi | [NASA-GAS; NASA-SP8112]. Başlangıç tank gazı ve regülatör kayıpları ayrıca gerekir. |
| `power_resources` | Ayrık ateşleme/bekleme/ön ısıtma evrelerinde `E=sum(P t)` | [BIPM-SI] W=J/s; [NASA-SE-6] EPS arayüzü. Her evre gücü kullanıcıdan kaynaklı gelir. |
| `check_requirements` | Kaynaklı eşiklere karşı artı/eksi marjin, eksik veride `unknown` | [NASA-SE-APP] Ek D. Sonuç yalnızca ön analiz kontrolüdür. |
| `weighted_trade` | Önce zorunlu sınırlar, sonra sabit uçlu ve kırpılmış doğrusal fayda | [NASA-SE-6] 6.8.1.2.3–4. Bu doğrusal ölçek uygulamanın açık seçimidir; ağırlıklar kullanıcı tercihidir. |
| `design_sensitivity` | Aynı tabandan açık Δv ve platform kütlesi senaryolarını yeniden çözme | [NASA-SE-6] 6.2.1.2.4/6.8. Olasılıksal dağılım veya katalog performans ölçeklemesi değildir. |
| `analyze_project`, JSON/HTML raporu | Tam girişler, katalog varyantı, sayısal sonuç, eksik bilgiler, kaynakça ve çıktı kimliği | [NASA-SE-APP] karar raporu; [NASA-SE-6] veri/konfigürasyon. HTML fizik hesaplamaz. |

Üretim paketindeki **özel yardımcılar ve sınıf metotları dahil her tanımlı
fonksiyon/sınıf**, `References` alanıyla kayıtlı bir kaynağa bağlanır. AST testi
bu koşulu denetler. Katalogdaki her sayısal alan ayrıca kendi birim ve kaynak
konumunu taşır; açıklayıcı katalog bilgileri de kayıt düzeyinde kaynaklıdır.

## 4. Rezerv ve kütle kapanışı

`D`: tüm kuru donanım dahil kuru araç kütlesi; `P`: yüklenecek toplam yakıt;
`f`: **yüklenen yakıtın** kullanılamayan kısmı; `R`: buna ek, uçuş sonunda
tutulacak sabit rezerv. Son kütle `D+fP+R` olur. Her yanmanın yararlı hız
katkısı için `c_i=g0 Isp_i eta_i` ve `S=sum(dv_i/c_i)` kullanıldığında:

```text
D + P = exp(S) (D + f P + R)
P = [D (exp(S)-1) + exp(S) R] / [1 - f exp(S)]
```

Bu, [NASA-ROCKET] kütle oranının açık rezerv tanımına cebirsel uygulanmasıdır.
Payda pozitif değilse sonlu çözüm yoktur. Küçük yanmalarda basamak kaybını
azaltmak için `expm1` kullanılır. Debi `[NASA-ISP]` üzerinden `F/(g0 Isp)`;
yön kaybı yararlı kuvveti düşürürken tank debisini gizlice düşürmez.

Sistem tasarım örneğinde kuru kütle `platform + diğer itki donanımı + kurulu
itici sayısı × katalog bileşen kütlesi` olarak toplanır. Diğer donanım tahsisi
tank, hat, elektronik, montaj ve tutulacak basınçlandırıcı envanterini
kapsamalıdır. Bir itici başlığının kütlesi bütün sistemin kütlesi olarak
sunulmaz. [NASA-SOA26; NASA-SE-6]

ADCS momentum boşaltma, sızıntı ve diğer net Δv üretmeyen sarfiyatlar bu
translasyon örneğine uydurma bir Δv olarak eklenmemelidir. Mevcut
`burns.vectoring`, `vehicle_geometry`, `pulse_mode` ve akış sözleşmesi bu
işler için kullanılabilir; tam görev bütçesinde ayrı tüketim olayı olarak
entegre edilmelidir. Bu örnek çalışma bunları içerdiğini iddia etmez.
[NASA-ISP; NASA-SE-6; bkz. `VEHICLE_GEOMETRY.md`, `PULSE_MODE.md`]

## 5. Katalog disiplini ve açık veri problemleri

Katalog beş ayrı donanım kaydı içerir: Moog 058E151 soğuk gaz, ArianeGroup
1 N hidrazin, S10-13 çift bileşenli, ECAPS 1 N HPGP ve Busek BIT-3. Kaynak
ve varyantlar kodda ve raporda birlikte saklanır.
[MOOG-COLD; ARIANE-1N; ARIANE-10N; NASA-SOA26; BUSEK-BIT3]

- Moog ve iki Ariane kaydı yalnızca etiketli nominal itki/Isp ile sabit
  çalışma noktası oluşturabilir. Bu noktaların sürdürülmesi ve görev
  koşullarına uygulanabilirliği ayrıca doğrulanmalıdır. [MOOG-COLD;
  ARIANE-1N; ARIANE-10N]
- Ariane 1 N için yayımlanan nominal `F`, `Isp` ve debi tam kapanmıyor.
  Kaynak debi ayrı tutulur; `F/(g0 Isp)` ile kapatılmış çalışma debisi ve
  farkı raporlanır. Fark, ölçüm belirsizliği diye adlandırılmaz.
  [ARIANE-1N; NASA-ISP]
- ECAPS kaydındaki aralıkların uçları eşleşmiş çalışma noktası değildir;
  otomatik ortalama veya en iyi uç seçilmez. [NASA-SOA26]
- Busek'in `/bit3` ve `/rf-ion-thrusters` sayfaları farklı üst performans
  sınırları veriyor. Tarihli kaynaklar ayrı tutulur; yeni bir ortak harita
  uydurulmaz. Bu nedenle elektrikli aday örnek karşılaştırmada `incomplete`
  kalır. Bu, elektrikli itkinin elenmesi anlamına gelmez.
  [BUSEK-BIT3; BUSEK-RF]
- Katalogda belirtilmeyen tam sistem kütlesi/gücü, yoğunluk, maliyet, TRL,
  itki haritası ve donanım ömrü hesaplanmış gerçek gibi doldurulmaz.
  [NASA-SE-6] Belirli rakam istenirse kaynağı eklenir veya sayının bir
  tasarım varsayımı olduğu açıkça yazılır.

## 6. Ekstrapolasyon politikası

**Yeni pakette katalog veya malzeme verisinin ekstrapolasyonu kapalıdır.**
Bir üst sınır, aralık ortası veya farklı ürünün değeri çalışma noktası diye
seçilemez. JSON'a `extrapolate` alanı veya `extrapolated` veri statüsü eklemek
hata verir. Mevcut performans yüzeyi modülleri de kendi geçerlilik sınırlarını
korur; SASS paketi bu kontrolleri gevşetmez.

İdeal gaz hesabı, açık varsayımlı bir analitik modeldir; gerçek gaz tablosunun
uzatılması değildir. [NASA-GAS] Deterministik senaryoda görev Δv'sini değiştirmek
de donanım performansı ekstrapolasyonu değildir. [NASA-SE-6]

İleride ekstrapolasyon gerekirse önce **asıl veri kaynağı ve geçerlilik
aralığı**, sonra **o veri türü/koşul için uzatma yöntemini destekleyen bağımsız
kaynak**, uygulama gerekçesi ve hata/geçerlilik sınırı sağlanmalı; bunlar
olmadan özellik açılmamalıdır. Bu, kullanıcının bu çalışma için koyduğu
gerekliliktir; genel bir NASA zorunluluğu olarak sunulmaz.

## 7. Çalıştırma ve teslimatlar

Komutlar güncel çalışma kopyasının kökünde çalıştırılır:

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) "src"
python -m rocket_propulsion.studies.sass examples/sass_mission.json `
  --output output/sass/study.json --report output/sass/study.html
python -m unittest tests.test_sass_studies -v
```

`examples/sass_mission.json` gerçek bir görev için kullanılacak hazır veri
değildir. Öğretim üyesinin isterleri geldikçe kaynak/varsayım zarfları içindeki
girişler değiştirilir. Program:

1. Üç nominal kimyasal/soğuk gaz adayını aynı görev tahsisleriyle hesaplar.
2. Eşleşmiş haritası olmayan iki adayı açık eksik veri durumunda tutar.
3. Kütle, süre ve açık güç tahsislerini eşiklerle kıyaslar.
4. Kullanıcı ağırlıkları ve sabit ölçeklerle koşullu karşılaştırma verir.
5. Girişlerin tamamını, kaynakları, yöntemleri ve varsayım soyunu JSON/HTML'ye yazar.

`design_sensitivity` ile farklı görev tahsisleri karşılaştırılabilir. Bu
senaryolar ağırlıklı bir olasılık dağılımı değildir. Belirsizlik çalışması
istenirse mevcut hot-fire/ensemble modülleri gerçek deney verisiyle
beslenmelidir; rastgele yüzdeler kalibrasyon sonucu gibi gösterilmemelidir.
[NASA-SE-6; bkz. `HOT_FIRE_CALIBRATION.md`]

## 8. Doğrulama ve gerçek isterler gelince

Yeni testler: bağımsız bisection ile rezervli/çok Isp'li kütle kapanışı,
roket denklemi, impuls/debi kimliği, düşük Δv duyarlılığı, ters Hohmann,
tek impulsla perigee düşürme, gaz uç envanterleri, EPS evre dengesi, kaynak
denetimi, atıfsız/bozuk veri reddi, eksik adayların sıralanmaması, senaryo
yeniden üretimi ve Python/CLI/JSON/HTML tutarlılığı.

Nihai proje yönergesi geldiğinde hâlâ öğrenilmesi gerekenler: görev ve yörünge,
ömür, kütle/güç/zaman sınırları, istenen doğruluk, izin verilen teknolojiler,
güvenilirlik/yeniden ateşleme beklentisi, donanım/prototip talebi ve teslim
biçimi. Bunlar araştırmayla uydurulmamış ve eksiklerinin altı çizilmiştir.
[ITU-UZM508; NASA-SE-APP]

Mevcut ayrıntılı motor, nozül, elektrikli itki, dinamik besleme, gimbal ve
Sidera aktarım modülleri bu hazırlığın devamında kullanılabilir. Tam yerel
Sidera tüketicisi veya bağımsız NASA CEA karşılaştırmasının tamamlandığı
iddia edilmez; ilgili mevcut geçerlilik sınırları sürer.

## References

Tam başlık, yayıncı, erişim tarihi ve kullanılan konumlar
`studies/sass/references.py` içinde; bütün URL'ler `docs/REFERENCES.md`
indeksinde saklanır. Ana kaynaklar:

- [ITU-UZM508: resmî Ninova ders kimliği](https://ninova.itu.edu.tr/en/courses/faculty-of-aeronautics-and-astronautics/37923/uzm-508e/)
- [ITU-UZM-GROUP: resmî ÖBS ders grubu](https://obs.itu.edu.tr/public/DersPlan/_DersGrupSearch?grupId=13265)
- [ITU-UZB419: ilgili lisans tasarım dersinin kataloğu](https://obs.itu.edu.tr/public/DersKatalog/DersKatalogBilgiBransDersKodu?bransKodu=UZB&dersNo=419)
- [NASA-SE-APP: sistem mühendisliği el kitabı ekleri](https://www.nasa.gov/reference/system-engineering-handbook-appendix/)
- [NASA-SE-6: teknik yönetim ve karar analizi](https://www.nasa.gov/reference/6-0-crosscutting-technical-management/)
- [NASA-ROCKET: ideal roket denklemi](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/)
- [NASA-ISP: özgül itki ve impuls](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/)
- [MIT-522-L6: ani ve düşük itki manevraları](https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/7f725e54b9be201164d56ebbd5e08023_MIT16_522S15_Lecture6.pdf)
- [NASA-GAS: ideal gaz ve yoğunluk tanımı](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/)
- [NASA-SP8112: basınçlandırma sistemleri](https://ntrs.nasa.gov/citations/19760015212)
- [NASA-SOA26: uzay içi itki derlemesi](https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/)
- [MOOG-COLD: üretici nominal koşulları](https://www.moog.com/products/propulsion/space-propulsion/spacecraft-propulsion/thrusters/cold-gas-thruster.html)
- [ARIANE-1N: üretici teknik karakteristikleri](https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html)
- [ARIANE-10N: üretici teknik karakteristikleri](https://www.space-propulsion.com/spacecraft-propulsion/bipropellant-thrusters/10-bipropellant-thrusters.html)
- [BUSEK-BIT3: tek ürün sayfası](https://www.busek.com/bit3)
- [BUSEK-RF: ayrı aile sayfası](https://www.busek.com/rf-ion-thrusters)
