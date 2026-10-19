# Motor ve besleme sistemi: doğrulama ve olgunluk dosyası

Bu dosya NASA danışmanlığı, NASA onayı veya sertifikasyon belgesi değildir.
Kullanıcının tanımladığı test nesnesi motor ve besleme donanımını birlikte kapsıyor.
Kullanıcı henüz tasarım seçmediğini bildirdi. Bu nedenle donanımın TRL seviyesi
**değerlendirilmedi**; TRL 5 başarılmış değildir. Yazılım testlerinin sonucu bu
durumu değiştirmez. Ders hazırlığı çıktıları yalnızca ön tasarım taraması içindir.

## Kabul temeli ve açık kapılar

NASA ESTO'nun donanım TRL 5 tanımında gerçekçi destek donanımıyla ilgili ortamı
temsil eden koşullarda çalışan bir brassboard ve kritik performansın gösterimi
istenir. Çıkış ölçütü analitik öngörülerle belgelenmiş uyumu ve sonraki aşamalar
için ölçekleme gereklerini içerir. Bu, bir Python test paketine verilen sertifika
olarak yorumlanamaz. [NASA-TRL]

NASA-STD-7009B'nin burada kullanılan sürümü PDF kapağına göre 2024-03-05 tarihli
B revizyonudur. Modelin kod doğrulaması, gerçek sistemle geçerlemesi, bunların
geçerli olduğu alanlar ve belirsizlik kaydı birbirinden ayrılmalıdır. Kabul
ölçütleri proje tarafından belirlenir; aşağıdaki çalışma NASA adına uygunluk
kararı vermez. [NASA-7009B, 4.1.1.5, 4.2.3–4.2.7, 4.3.8]

| Kapı | Şimdiki durum | Kapatmak için gereken kanıt |
|---|---|---|
| Konfigürasyon | Açık | Seçilmiş motor/akışkanlar, P&ID, BOM, arayüzler, revizyon |
| İlgili ortam | Açık | Görevden türetilen basınç/sıcaklık/vakum/darbe profilleri ve eşzamanlı yükler |
| Kabul ölçütleri | Açık | Onaylı performans gereksinimleri, hata payları, belirsizlik ve karar kuralı |
| Test nesnesi | Açık | Kimliği ve üretim kayıtları belli motor + besleme brassboard'u |
| Metroloji | Açık | Kalibrasyon izleri, DAQ zaman eşlemesi, sensör/işleme belirsizlik bütçesi |
| Fiziksel gösterim | Açık | Entegre nesnenin ilgili ortamda ölçümlü çalışması ve ham kayıtlar |
| Model geçerlemesi | Açık | Kalibrasyonda kullanılmamış veriyle, sınırları belirtilmiş karşılaştırma |
| Ölçekleme | Açık | Analiz-deney uyumu ve sonraki tasarıma aktarım gerekçesi |
| Teknik inceleme | Açık | Bağımsız inceleme, uygunsuzlukların kapatılması, olgunluk değerlendirmesi |

Bu tablo bizim proje kapı listemizdir; resmî NASA gereksinim matrisi veya eksiksiz
standart uygunluk incelemesi değildir. Sayısal test basınçları, ateşleme reçetesi
ve evrensel başarı yüzdesi tasarım seçilmeden uydurulmamıştır.

## Bulunan ve giderilen gerçek kusurlar

Serbest uzay hesabındaki pozitif Δv'nin sıfıra dönüşmesi, yanıltıcı sıfır
hata farkı ve tablolu integral/moment taşması
[sayısal çekirdek incelemesinde](FREE_SPACE_TABULATED_NUMERICS_REVIEW.md)
ele alınır. Yönlü sayısal aralıklar ölçüm belirsizliği veya fiziksel geçerleme
değildir. Kısa segmentlerin göreli sahipliği, kaydırılmış koordinatı ve kanonik epoch
olayları [zaman incelemesinde](TABULATED_TIME_CONTRACT_REVIEW.md) ele alınır.
Ham belge sözleşmeleri, native tüketim/rezerv aritmetiği ve bütün fiziksel
geçerleme bununla kapanmış sayılmaz.

Yanma kayıtlarının scalar/enum/bayrak doğrulaması, değişebilir alt listeleri,
JSON ham türleri ve doğrulanmamış sonuç yazımı
[kayıt sözleşmesi incelemesinde](BURN_RECORD_CONTRACT_REVIEW.md) ele alınır.
Aynı inceleme, fonksiyon envanterinin iç kontrol akışlarını ve lambda'ları
atlayan taramasını da onarır. Tam fiziksel geçerleme sonucu çıkarılmaz.

Doğrudan tank envanteri yardımcılarındaki eksik eşleme, geçersiz zaman,
akış toplama ve rezerv kararı kusurları
[ayrı envanter yardımcıları incelemesinde](DIRECT_INVENTORY_HELPERS_REVIEW.md)
ele alınır. L0/L1 bu modülün yalnız kapsam denetimini kullanır; süre/kütle
yardımcılarıyla çözücü hesaplarının çağrı yolları ayrı belirtilir.

Native Lunaris dürtü normu/kütle güncellemesindeki aralık kusuru için
[ayrı yama ve kanıt incelemesi](NATIVE_LUNARIS_IMPULSE_REVIEW.md) tutulur.
Aday kaynağın ve canlı uygulamanın durumları ayrıdır; dosya hash'i, bağımsız
referanslar ve tam koşu olmadan canlı düzeltme tamamlanmış sayılmaz.

Sidera sabit/ensemble ve Lunaris dürtü adaptörlerinde büyük/küçük yönlerden
sıfır veya tanımsız vektör üreten aritmetik, ayrıca ensemble yeniden açılışında
koordinat/hash kayması [adaptör yön incelemesinde](ADAPTER_DIRECTION_NUMERICS_REVIEW.md)
ele alınır. Native kontroller sabit yanma/dürtü sınırıdır; ileri weighted
provider yayılımı veya fiziksel yönlendirme geçerlemesi sayılmaz.

L1 profil yolundaki pozitif hedefe sıfır çıktı verip başarı bildirme, komut
terslemesi, moment ve kesme/envanter hataları
[profil incelemesinde](PROFILE_BURN_NUMERICS_REVIEW.md) ele alınır. Bütün
çıktılar için evrensel hata garantisi veya fiziksel geçerleme çıkarılmaz.

L0 sabit yanma çözücüsünde küçük Δv terslemesi/ileri hesabı, kesme sınırını
aşan hedef toleransı ve yuvarlanmış tank durumlarından kaybolan tüketim
[sabit yanma incelemesinde](CONSTANT_BURN_NUMERICS_REVIEW.md) ele alınır.
L1 profil çözümü bu incelemeyle doğrulanmış sayılmaz.

Tank envanteri yardımcılarında artık iptali ve yuvarlanmış artığın toleransla
karşılaştırılmasından doğan yanlış kapanışlar ayrıca giderildi. Üç modelin
ortak sayısal defteri, JSON'da korunabilen tam artık oranı ve bağımsız
kontroller [envanter kapanışı incelemesinde](INVENTORY_CLOSURE_NUMERICS_REVIEW.md)
açıklanır. Bu sonuçlar donanımın kütle bütçesini ölçümle doğrulamaz.

İlk sekiz regresyon testi eski uygulamada başarısız oldu. Komşu yörüngeler testi
alt kontrolleriyle birlikte on başarısızlık bildirdi. Düzeltmelerin amacı:

| Kusur | Etki | Düzeltme / saldırgan kontrol |
|---|---|---|
| Komşu yarıçaplarda iki büyük hızın çıkarılması | Sıfır/yanlış işaretli küçük manevra | Rasyonelleştirilmiş fark; 75 basamaklı Decimal oracle ve komşu float girdiler |
| Yakıt tüketimini başlangıç kütlesinden çıkarma | Büyük kütle oranında son kuru kütle sıfırlanıyor | Son kütleyi üstel bağıntıyla ayrı hesaplama; pozitiflik ve kapanış denetimi |
| Küçük oksitleyici miktarını toplamdan çıkarma | Pozitif bileşen sıfır oluyor | İki bileşenin doğrudan oranla hesaplanması |
| Küçük gaz artığını toplamdan çıkarma | Depoda kalan gaz kayboluyor | Uç basınçlardan doğrudan artık hesabı |
| Dondurulmuş sonuç içinde değişebilir listeler | Hesap sonrası sonuç değiştirilebiliyor | İç içe kayıtları tuple'a dönüştürme |
| Varsayımı published olarak yeniden etiketleme | Kaynak niteliği yanlış sunulabiliyor | Varsayım soyunu ve salt kayıt-yöntemi atıflarını published için reddetme |
| Yinelenen JSON alanı | Son alan önceki kanıtı sessizce eziyor | Her derinlikte yinelenen alanları reddetme |
| Devasa donanım adedi | Sonraki hesapta denetimsiz taşma | Hesaba girmeden sonlu temsil denetimi |
| Tek darbe hesabının gereksiz yere periyot hesabına bağlanması | Periyot taşınca sonlu darbe hesabı da reddediliyor | Tek darbe bağıntısını doğrudan hesaplama; yapay uç ölçekli Decimal oracle |

Son satırdaki dokuzuncu kusur daha sonraki denetimde ayrıca yeniden üretildi.
Bu uç sayısal durum gerçek bir görev önerisi değildir.

Sonraki kütle incelemesi ayrıca uygulanabilirlik sınırında %50 yakıt hatası,
son çıktılar sonlu iken ara kütle oranı taşması ve pozitif yanmanın sıfır
raporlanması kusurlarını yeniden üretti. Hassasiyet kontrollü hesap, ileri
yayılım oracle'ı ve artık/koşulluluk kanıtları
[kütle kapanış incelemesinde](MASS_CLOSURE_REVIEW.md) açıklanır.

Faz enerjisinde ayrıca yakın sürelerin çıkarımından doğan bekleme enerjisi
hatası ve ayrı çarpımların alt taşmasından kaybolan toplam enerji yeniden
üretildi. [Enerji defteri incelemesi](POWER_LEDGER_REVIEW.md) düzeltmeyi ve
bağımsız kontrolleri açıklar.

Depolama incelemesi de sonlu hacim/envanteri bozan ara taşma ve alt taşmayı,
pozitif ullage'ın sıfır raporlanmasını yeniden üretti. Başlangıç/son gaz
envanterinden bağımsız hacim kökü ve son çıktı temsil sınırı
[depolama incelemesinde](STORAGE_NUMERICS_REVIEW.md) açıklanır.

Karar zincirinde eşit ağırlıklı toplamların farklı puan alması ve pozitif
görev ölçeklemesinin sıfır tahsise dönüşmesi ayrıca yeniden üretildi.
[Karar hesabı incelemesi](DECISION_NUMERICS_REVIEW.md) düzeltmeyi açıklar;
fiziksel belirsizlik altında sıralama dayanıklılığı henüz değerlendirilmedi.

Deney karşılaştırmasının raporladığı hata uçlarının tam aralığı içeri
daraltabildiği de yeniden üretildi. [Aralık incelemesi](INTERVAL_VALIDATION_REVIEW.md)
dışa yuvarlamayı, bağımsız köşe oracle'ını ve sayısal karar korumasını açıklar.
Sayısal karşılaştırma pass'i fiziksel geçerleme veya kalibrasyon onayı değildir.

Ek olarak rapor yazımı dosya takma adlarını/sabit bağları denetler, iki çıktıyı
önce hazırlar ve olağan yazma hatasında önceki dosyaları geri yüklemeyi dener.
İkinci dosyanın değiştirilmesinde hata enjekte edilerek ilk dosyanın korunduğu
kontrol edilir. Farklı dosya sistemleri arasında çift dosyalı bir işlem için
güç kesintisine karşı atomiklik iddia edilmez; eşzamanlı CLI yazarları desteklenmez.
Yalnızca içerik özeti kullanmak dosyanın yazarını veya doğruluğunu kanıtlamaz.

## Fonksiyon ve model inceleme sınırları

[Çok doğrusal harita incelemesi](MULTILINEAR_NUMERICS_REVIEW.md), kuvvet/toplam
debi ve bağımsız akışkan kanallarının ağırlık, gradyan ve köprü defterini kapsar.
Sayısal harita tutarlılığı fiziksel çalışma alanı geçerlemesi değildir.

[PCHIP ve kübik harita incelemesi](PCHIP_NUMERICS_REVIEW.md), ortak çekirdeği
ve iki tüketicide değer/türev aritmetiğini kapsar. Çok doğrusal harita yolu
ve fiziksel performans geçerlemesi bu sınırlı yazılım sonucundan çıkarılamaz.

[Gaz kelebeği debi incelemesi](THROTTLE_FLOW_NUMERICS_REVIEW.md) tüketim,
türevler, özgül itki ve kütle kapanışını verilen eğri örnekleri üzerinden
ele alır. Eğri çekirdeğinin tamamı için doğruluk iddiası değildir.

[Yön geometrisi incelemesi](DIRECTION_NUMERICS_REVIEW.md) üç tüketicide ölçek,
işaret ve kestirim türevlerini kapsar; donanım yöneltme geçerlemesi sağlamaz.

`output/assurance/audit.json` `src/rocket_propulsion` altındaki Python kaynak
dosyalarının SHA-256 özetini, açık Python
fonksiyon/metot/sınıf/async-function/lambda bildirimlerinin konumlarını ve her
öğenin inceleme kapsamını içerir. Kontrol akışı ve imza içleri de taranır.
Dinamik/üretilmiş callable'lar, bu ağaç dışındaki script'ler ve diğer diller bu AST listesinde yer almaz.
Listede görünmek incelenmiş olmak anlamına gelmez; işaretli kapsam dışında
yeniden denetim veya bütün kütüphane doğruluğu iddia edilmez.
Atıf varlığı, formülün veya donanım uyumunun kanıtı değildir.
Standarttaki seçilmiş maddelere karşı açık/eksik kanıt eşlemesi
`docs/assurance/nasa7009b_crosswalk.json` içindedir. Tasarım seçildikten sonra
tamamlanacak kayıt alanları `examples/hardware_validation_manifest.template.json`
dosyasında boş bırakılmıştır; bu şablon deney kanıtı sayılmaz.

| Fonksiyon / grup | Kabul edilen kullanım | Kalan fiziksel / teknik sınır |
|---|---|---|
| `hohmann_transfer` | İki cisim, eşdüzlemli dairesel uç yörüngeler, ani itki | Sonlu yanma/pertürbasyon geçerlemesi yok |
| `circular_spiral_delta_v` | Küçük teğetsel ivmeli yaklaşık dairesel spiral bütçesi | Güç, gölge, uçuş süresi ve optimal kontrol çözümü yok |
| `perigee_lowering_delta_v` | Tek retrograd ideal darbe | Atmosferik ömür veya bertaraf uygunluğu yok |
| `burn_duration_screen` | Kullanıcının seçtiği süre/orbit oranı sınırına tarama | Bu oran bir yörünge hata sınırı veya evrensel kabul ölçütü değil |
| `delta_v_budget`, `DeltaVItem` | Açık adet ve marjla kaynak tahsisi | Gerçek istasyon koruma tarihçesi üretilmiyor |
| `size_mission`, `Maneuver` | Verilen sabit sistem Isp/F ve eksenel verimle sıralı kütle bütçesi | Isp çalışma noktası ölçülmeli; valf/ısı/besleme geçişleri çözülmüyor |
| `bipropellant_split` | Sabit toplam boşaltım O/F oranı | Ayrı tank artıkları ve zamanla değişen karışım çözülmüyor |
| `power_resources` | Ayrık çalışma/boşta/ısınma fazlarının enerji defteri | EPS, PPU, termal ve batarya geçerlemesi yok |
| `liquid_tank_volume` | Verilen koşula ait yoğunluk ve açık ullage | Tank yapısı, PMD, çalkantı ve yoğunluk tablosu yok |
| `ideal_cold_gas_storage` | Sabit T, ideal gaz ve mutlak basınçla envanter | Yüksek basınç gerçek gazı ve blowdown ısı geçişi doğrulanmadı |
| `ideal_regulated_pressurant` | Aynı sabit sıcaklıkta ideal şişe/ullage envanteri | Başlangıç ullage gazı, regülatör geçişi/kayıpları ve çözünme yok |
| `nominal_operating_point` | Üretici nominal F/Isp'sinden açıkça varsayımsal sabit çalışma hesabı | `require_validated=True` tüm mevcut kayıtları reddeder |
| `check_requirements`, `weighted_trade` | Deterministik ön tasarım taraması ve tercihler | Belirsizlik içeren donanım uygunluğu veya mimari seçimi kanıtı değil |
| `analyze_project`, JSON okuyucular | İzlenebilir varsayımlar ve yeniden üretilebilir senaryo | Tekrarlı manevra adetleri toplam bütçeye birleştirilir; darbe/başlatma ömrü denetlenmez |
| `design_sensitivity` | Bağımsız açık senaryoları yeniden hesaplama | Olasılık dağılımı veya güven aralığı çıkarımı yok |
| Kanıt/katalog/sunum yardımcıları | Birim, kayıt, girdi ve çıktı bütünlüğü | Kaynak sahibinin beyanını doğrulayan bir otorite değiller |
| `compare_validation_sample`, `assess_validation` | Kaydedilmiş bağımsız deney ve öngörü aralıklarını karşılaştırma | Dosya kimliği denetlenir; deneyin gerçekliği ve onay kayıtları otomatik doğrulanmaz |

Kütle oracle'ı 250 sabit tohumlu senaryoda 75 basamaklı Decimal aritmetiğiyle
bağımsız geri yayılım ve kök araması yapar; üretimdeki kapalı formülü kopyalamaz.
Uygulanamaz artık kapanışları da ret davranışıyla sınanır. Toleranslar **sayısal
doğrulama toleransıdır**, fiziksel model doğruluğu veya deney toleransı değildir.
SASS mutasyon aracı on altı seçilmiş anlamsal kusuru bellekte enjekte eder; yalnızca
assertion başarısızlığını yakalanmış kusur sayar. İlişkisiz istisna başarı değildir.
Bu çalışma bütün olası mutasyonları taramaz.

## Katalog kayıtlarının yeniden değerlendirilmesi

| Kayıt | Denetim sonucu | Endüstriyel kullanım için açık konu |
|---|---|---|
| Moog 058E151 | Sütun/model, N2, F/Isp, mutlak giriş basıncı, bileşen kütlesi ve ayrı açma/tutma gücü kontrol edildi | Ortak çalışma koşulları, ölçüm belirsizliği ve entegre besleme kanıtı yok |
| Ariane 1 N | Kaynakta 1 N, 220 s, 0.44 g/s ayrı nominal değerler olarak bulunuyor | F/(g0 Isp) ile yayımlanan debi arasında yaklaşık %5.34 fark; belirsizlik yok, uzlaştırılmadı |
| Ariane S10-13 | Tek valf kütlesi çift valfle birleştirilmedi; O/F aralık olarak tutuldu | Küçük debi tutarsızlığı da otomatik belirsizlik kabul edilmedi |
| ECAPS 1 N HPGP | NASA tablosundaki ayrı aralıklar çalışma noktasına dönüştürülmedi | Birlikte ölçülmüş F/Isp/güç ve revizyonu belli üretici verisi gerekli |
| Busek BIT-3 | `/bit3` ve aile sayfasındaki çelişki korundu | Maksimumlar birlikte erişilebilir nominal nokta değil; revizyon/çalışma haritası açıklaması gerekli |

JSON ve HTML çıktıları açıkça `hardware_readiness=not_assessed`,
`experimental_validation=not_performed`, `qualification_evidence=false` taşır.
Kaynak alıntısı veya varsayım etiketinin bir saldırgan tarafından tamamen
yeniden yazılmasını bu veri modeli tek başına engelleyemez; dış veri onayı ve
kontrollü kaynak deposu gerekir. `published` durumu otomatik doğrulanmış demek değildir.

## Bağımsız veri araştırması

İndirilen dört PDF'nin orijinal byte özetleri
`output/assurance/sources/manifest.json` içinde tutulur. Büyük kaynak PDF'leri
paketin çalışma zamanı bağımlılığı değildir. Deney özeti ve kabul edilmeme
gerekçeleri `docs/assurance/research_sources.json` içindedir.

- **SAMSON 2019:** s.10–11, Şekil 12, azotla atmosferik ölçümlerin ortalamaları
  60 psig için 22.08 mN, 30 psig için 9.06 mN. Vakumdaki kripton için yaklaşık
  20 mN değeri aynı koşulda doğrudan ölçüm değil, düzeltme hesabıdır. Ham zaman
  serisi ve sayısal toplam ölçüm belirsizliği bu incelemede bulunmadı. Metindeki
  “43.5 psig (2 bar)” eşlemesi de açıklama gerektiriyor. Bunlar bizim motorun
  bağımsız geçerleme verisi olarak kabul edilmedi. [SAMSON-2019]
- **GFSSP v6:** NASA/TP-2016-218218, 884 PDF sayfası. Örnek 29, basılı s.293–301,
  MHTB hidrojen tankı basınçlanmasını deneyle karşılaştırıyor. Isı yalıtımı
  bozunma katsayısının aynı testle ayarlanması bağımsız holdout kanıtı yerine
  geçmez. GFSSP'nin çalıştırıldığı veya bizim düşük mertebeli modele karşı
  doğrulama tamamlandığı iddia edilmiyor. [NASA-GFSSP6]
- **NASA/JSC 20170000809:** indirilen dosya yalnızca tek sayfalık özet.
  Deneylerin yapıldığını bildiriyor; kullanılabilir ham ölçüm serisi sağlamıyor.
  Arama aracının otomatik PDF başlığı yanıltıcıydı; yerel kapağa göre kimlik
  kontrol edildi. Karşılaştırma için kabul edilmedi.
- **DLR DX2025/LUMEN:** Bölüm 4.3'te sağlanan veri, deneylerle sınanmış bir
  simülatörün çıktısı olarak açıklanıyor. Gerçek sistem adıyla verilen ayrı
  simülatör trajektorileri hot-fire ham ölçümü değildir. Algoritma karşılaştırma
  adayıdır; bizim motorun deney kanıtı değildir.

## Yeniden üretme ve kayıtlar

PowerShell, bu worktree kökünde:

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) 'src') + ';' + (Get-Location)
python -m unittest tests.test_sass_studies tests.test_sass_adversarial tests.test_sass_validation
python scripts/run_sass_mutations.py
python scripts/assurance_audit.py
python scripts/run_bounded_regression.py --include-deferred
python -m ruff check src tests examples scripts --no-cache
python -m rocket_propulsion.studies.sass examples/sass_mission.json --output output/sass/study.json --report output/sass/study.html
```

Geniş regresyon koşusunun ayrıntısı JUnit XML'de tutulur. Tam paket geçti iddiası
için atlanan/çalıştırılmayan kontrollerin de kapanması gerekir. Önceden uzun süre
tamamlanmayan iki gerçek Lunaris testi ve iki rastgele pulse planlayıcı testi
ayrı izlenir; dışarıda bırakıldıklarında yeşil sayılmazlar. Windows sandbox'ında
sabit bağ oluşturma izni yoksa gerçek hardlink deneyi atlanır; `samefile` ret
dalı ayrı hata enjeksiyonuyla kontrol edilir. Rapor HTML'sinin yapısal/kaçış
kontrolleri vardır; bu turda tarayıcıda görsel doğrulama iddia edilmiyor.

Son regresyon aracı her modülü ayrı süreçte, en çok 120 saniye çalıştırır;
`--include-deferred` önceki dört uzun kontrolü de kapsar. Num­ba önbelleği açıkça
çalışma alanındaki `output/assurance/regression/numba-cache` dizinine yönlendirilir.
Yerel HTTP/sabit bağ işlemleri sandbox tarafından engellenirse aynı kontroller
izinli yerel koşuda yürütülür; izin hatası model hatası veya başarılı test sayılmaz.
İlk 60 saniyelik sınırda pulse modülünün dokuz testi 58.04 saniyede tamamlanmış,
ancak süreç kapanışı süreyi aşmıştır; daha uzun sınırlı son koşu ayrı kayıttır.
Güncel sonuç `output/assurance/regression/summary.json` içindedir. Kaynaklar
koşu sırasında değişirse `complete_suite_passed` doğru olamaz.

Bir C-tabanlı yığın izi zamanlayıcısı çöküşü ve eksik JUnit'in eski yeşil
özeti bırakması ayrıca yeniden üretildi. [Tanılama/kanıt incelemesi](DIAGNOSTIC_EVIDENCE_REVIEW.md)
kontrol alt süreçlerini ve düzeltmeyi açıklar. İz toplama Python callback'ine
taşındı; 120 saniyelik sert modül sınırı ve bütün testler korunur. Eksik
kanıt güncel başarısız rapor ve sıfır olmayan komut çıkışı üretir.

Enerji düzeltmesinden sonraki ilk geniş koşuda 466 metot geçti, HTTP yetenek
sorgusu iki saniyelik soket sınırını aştı. Başarısız çıktı
`output/assurance/http_parallel_failure.txt` ve `verification_parallel_attempt.json`
olarak korundu. Aynı Sidera yolu ve aynı süre sınırıyla tek başına HTTP modülünün
altı testi geçti (`http_isolated.xml`). Nedeni yalnız bu iki gözlemle kesinleşmez.
İzolasyona rağmen kusur daha sonra tekrarlandı (`http_isolated_failure.*`);
bu önceki protokol kök sorunu çözmüş sayılmaz. Son düzeltme isteğe bağlı
paket incelemesini HTTP yanıt yolundan çıkarır ve `checking` durumunu açık
gösterir. Son koşu HTTP dahil bütün modülleri dört işçiyle çalıştırır.
[HTTP incelemesi](HTTP_CAPABILITY_PROBE.md) API sürümünü, hata durumlarını ve
testte işin tamamlanmasını gözleme penceresini açıklar. İstek başına iki
saniyelik soket sınırı korunur. Genel yük altındaki HTTP gecikme dağılımı
raporda `not_assessed` kalır.

İlk SASS incelemesinin tarihsel koşusunda 444 test metodu toplandı;
temel yapılandırmada 443 geçti, bir isteğe bağlı yerel
Sidera kontrolü atlandı. `D:\sidera\src` yalnızca testin PYTHONPATH'ine açıkça
eklenerek yapılan ek koşuda beş Sidera testi de geçti; atlanan yöntem böylece
ayrıca çalıştırıldı. Bu yerel yol üretim koduna eklenmedi. Dolayısıyla 444
yöntemin tümü için ayrı koşularda başarılı kanıt toplandı. Güncel sonuçlar
`output/assurance/verification_summary.json` içinde JUnit ve kaynak özetlerinden
yeniden üretilir; tarihsel sayılar güncel koşunun yerine kullanılmaz.
Sidera sonucu sabit yanma public adaptörüyle sınırlıdır;
gelişmiş yerel provider tüketicisinin tamamlandığını kanıtlamaz.
JUnit toplamı alt testleri de saydığından 444 yöntem sayısıyla karıştırılmamalıdır.

SASS için statement kapsamı 889/956 (%92.99), branch kapsamı 227/288 (%78.82).
Kapsanmayan yollar ve deneysel geçerleme eksikliği açık sınırlardır. Yeşil
regresyon veya bu yüzdeler donanım TRL düzeyi belirlemez.

Sonraki motor/besleme modeli incelemesi, ek kusurlar ve korunum oracle'ları:
[Besleme fiziği incelemesi](FEED_PHYSICS_REVIEW.md). Görev gereksinimleri de
henüz tanımlanmamıştır; donanım ve ilgili ortam seçimine ilişkin kapılar açıktır.

Bağlı modelin negatif olay denemelerinde tank türevini kırpıp yanlış sıfır
kütle kapanışı raporlaması ayrıca yeniden üretildi ve düzeltildi:
[olay uzantısı incelemesi](COUPLED_FEED_REVIEW.md). İşaretli matematiksel devam,
fiziksel çözümün sıfır-akış olayında durması şartını kaldırmaz. Darbe
planlayıcısının aynı kapsamlı rastgele deneyi önceki paralel koşuda 120 saniyeyi
aştı. Profil, genişlik aramasında gereksiz çift eğri değerlendirmesini gösterdi;
yalnız kütle değerlendirmesi, model başına en çok 512 skaler ters çözüm ve
çağrı ömrüyle sınırlı darbe önbelleği kullanılarak
iyileştirildi. Test girdileri, oracle'lar, toleranslar, 120 saniyelik modül
sınırı ve dört paralel işçi değiştirilmedi. Güncel geniş koşu ayrıca kaydedilir;
tekil profil süresi genel yük altında gecikme garantisi değildir.

Aynı bağlı modelde ortak küçük debi ölçeklerinin karışım oranı türevini
sıfırlaması ve ham bölme hatası üretmesi ayrıca düzeltildi. Yakın iptal ve
temsil edilebilir çıktı kontrolleri [bağlı model incelemesinde](COUPLED_FEED_REVIEW.md)
açıklanır. Genel yazılım durumu artık geçen pytest sonuçlarının yanında bütün
tanımlı seçilmiş mutasyonların da yakalanmasını ister; durum ve kimlik
tutarlılığı [kanıt raporu incelemesinde](DIAGNOSTIC_EVIDENCE_REVIEW.md) açıklanır.
Bu daha sıkı genel sonuç, fiziksel geçerleme veya test yeterliliği kararı değildir.

Köprüdeki kanal/tank debi defteri ayrıca [ayrıntılı olarak incelendi](BRIDGE_FLOW_NUMERICS_REVIEW.md).
Temsil edilebilir küçük toplamlar ve karşıt terimler korunur; temsil
edilemeyen ihraç edilen türev açık hata verir. Bu kapsam, bütün kuvvet/yön
hesaplarının veya propagatörün tamamının incelendiği anlamına gelmez.

## References

- [NASA-TRL](https://esto.nasa.gov/trl/)
- [NASA-STD-7009B, 2024-03-05](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [SAMSON-2019](https://www.mdpi.com/2226-4310/6/8/91)
- [GFSSP v6, NASA/TP-2016-218218](https://ntrs.nasa.gov/citations/20160005401)
- [NASA/JSC özet kaydı](https://ntrs.nasa.gov/citations/20170000809)
- [DLR DX2025, Bölüm 4.3](https://elib.dlr.de/219953/1/DX2025benchmark.pdf)
- Formül, katalog ve birim kaynakları: [kaynak dizini](REFERENCES.md),
  [ilk model ve varsayım açıklamaları](SASS_PROJECT_PREPARATION.md).
