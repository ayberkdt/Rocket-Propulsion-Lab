# Motor ve besleme bağlantısında olay uzantısı

Bu inceleme, tek ve çift hatlı bağlantıda sıfır-akış olayının iki tarafındaki
sayısal deneme durumlarını kapsar. Kullanılan tank, motor haritası ve akışkan
değerleri sentetiktir. Donanım geri akış davranışı veya deneysel geçerleme
kanıtı değildir.

## Yeniden üretilen kusur

Dinamik hat negatif akışı olay kökü araması için tutarken, bağlantı tank
talebini `max(0, q)` ile sıfırlıyordu. Hat depolama türevi ise `q - q_engine`
olarak kalıyordu. Böylece döndürülen tank ve hat türevleri birlikte korunumu
sağlamıyordu. Çift hatlı bağlantının hata raporu gerçek tank türevi yerine
yeniden `-q` yazdığı için bu açığı sıfır raporluyordu.

Oksitleyici hattında `q=-0.1 kg/s` örneğinde rapor sıfır, döndürülen türevlerin
defteri `-0.10000000000000009 kg/s` idi. Motor açık/kapalı, tek/iki negatif
bacak ve sıfır-akış kontrollerinden oluşan sekiz kaydın altısında uyuşmazlık
vardı. Tarihsel kaynak hash'leriyle kayıt:
`output/assurance/coupling_reverse_flow_probe.json`. Bu dosya düzeltilmiş
kaynağın başarılı testi diye kullanılmaz. İlk yeni regresyon yönteminin
eski kodda verdiği gerçek kütle-türevi assertion hataları ve yeni alan henüz
bulunmadığı için oluşan istisnalar `coupling_before_fix.xml` içindedir;
bunlar ayrı ayrı fiziksel kusurlar diye sayılmaz.

## Düzeltmenin sözleşmesi

`RegulatedFeedSystem.evaluate(..., for_event_location=True)` sıfır statik
direnç şartıyla negatif akışın matematiksel devamına izin verir. Varsayılan
çağrı negatif talebi reddetmeye devam eder. İşaretli akış hem tank sıvı
türevi `-q` hem `dV_ullage/dt=q/rho` için kullanılır. Negatif deneme durumu
böylece gazı sıkıştırır; sınır işi aynı işaretle enerji denklemine girer.

Sonuçta `reverse_flow_event_extension` alanı, tek/çift bağlantıda da aynı
adlı özellik, negatif durumları açıkça ayırır. Fiziksel çözüm sıfır-akış
kesme olayında durmalıdır. Bayrak kendi başına olayı uygulayan bir
propagatör değildir ve ters çalışan bir donanım modeli oluşturmaz.
Tank manifesti bu politikayı içerir; bağlı model hash'i de değişir.
Manifest hash'i çalıştırılabilir kaynak kodunun hash'i değildir; test kanıtı
kaynak hash'lerini ayrıca tutar.

Çift hatlı kütle kapanışı, döndürülen tank ve basınç türevleri ile motorun
döndürdüğü kütle debileri üzerinden hesaplanır. Her hattın depolaması
`C * dp_manifold/dt` ile yeniden kurulur. Oksitleyici ve yakıt motor
debilerinin toplamı, iki hattın ayrı kütle dengeleri ve araç dengesi dört
ayrı kontroldür. Kabul, verilen binary64 sayıları temsil eden tam rasyonel
değerler üzerinde yapılır. Toleransa eşit artık kabul edilir. Raporlanan
artık daha sonra binary64'e yuvarlanır.

### Kapanış kararının sayısal sınırları

Sentetik aktif motorda 1,3 ve 1,5 MPa manifold basınçlarıyla bileşen debileri
2,16 ve 0,8625 kg/s, toplam debi 3,0225000000000004 kg/s oldu. Eski
toplama işlemi farkı sıfır gösterdi. Sağlanan değerlerin tam farkı
yaklaşık 2,22 × 10⁻¹⁶ kg/s olduğundan 10⁻¹⁶ kg/s toleransında
bileşen kontrolü durmalıdır.

Tanklardan hatlara 0,11 ve 0,30 kg/s giderken eski araç artığı
2,91 × 10⁻¹⁶ kg/s idi. Basınç türevlerinden hesaplanan araç artığı
8,40 × 10⁻¹⁷ kg/s. Her hat da 2 × 10⁻¹⁶ kg/s sınırını sağlar. Buna
karşılık motor kapalı, hat akışları +0,10 ve −0,10 kg/s iken araç artığı
sıfır olsa da her bacağın mutlak artığı yaklaşık 3,46 × 10⁻¹⁸ kg/s.
1 × 10⁻¹⁸ kg/s toleransında bu durum reddedilir. Negatif akış burada
sıfır-akış olayını aramak için kullanılan bir deneme durumudur.

Bağımsız Decimal denetimleri motor açık ve kapalıyken döndürülen
türevlerden artığı hesaplar. Toleransa eşit ve komşu float sınırları,
her bacağın ayrı reddi ve araç toplamının ayrı reddi sınanır. Bir
sentetik durumda tam artık, ekranda gösterilen toleransı yaklaşık
10⁻¹⁰⁰ kg/s aşar. Rapor ve tolerans float olarak eşit görünürken kabul
kararı tam değeri kullanır.

Manifold basınç Jacobian'ındaki dört motor-basınç terimi ve tank
basıncından hat debisine olan duyarlılıklar tam bölümden sonra
aktarılır. Ölçek testinde sıfır olmayan bir duyarlılık binary64
aralığının dışında kalırsa alan adıyla hata döner. Manifest
`mass_closure_numerics` alanında kapanış politikasını belirtir.
Bu sayısal sınırlar ölçüm toleransı olarak kullanılmaz. [NASA nodal
besleme modellemesi](https://ntrs.nasa.gov/citations/20240003493)

### Kriyojenik tank deneyinde ölçütlerin ayrılması

NASA TN D-7638'deki 1,52 m çaplı metan tankının altı statik helyum
koşusu, 2025 tarihli GFSSP çalışmasının 3–6 numaralı tablolarında
yeniden kullanılmıştır. Aynı çok düğümlü modelin basılı değerlerinden
yeniden hesaplanan helyum tüketimi ortalama mutlak yüzde hatası %4,43.
İlk basınçlandırmada metan yoğuşması ortalama %73,73 düşük, boşaltmada
net buharlaşma ortalama %120,32 yüksek tahmin edildi. Bekletmede 36 ve
40 numaralı koşuların raporlanan faz değişimi işareti modelle ters.
Çalışma bu iki küçük ölçüm için belirsizlik tahmini vermez.

Bir kriyojenik konfigürasyon seçilirse, toplam basınçlandırıcı kütlesi
ile basınçlandırma, bekletme ve boşaltma evrelerindeki faz değişimi
ayrı kanallarda karşılaştırılır. Bu kıyas başka bir tesisteki GFSSP
sonuçlarını yeniden hesaplar. Bu depodaki regülatör modeline veya
henüz seçilmemiş donanıma deney geçerlemesi atamaz.
[Deney raporu](https://ntrs.nasa.gov/citations/19740019160)
[Model karşılaştırması, Tablolar 3–6](https://ntrs.nasa.gov/citations/20240014846)

## Bağımsız kontroller

- Her bacakta depolama `C * dp_manifold/dt` üzerinden yeniden kurulur;
  tank türevi, bu depolama ve motor debisi ayrı ayrı kapatılır. Normal,
  sıfır, negatif, karşıt işaretli ve çok küçük negatif akışlar; motor açık
  ve kapalıyken sınanır. Raporlanan araç hatası gerçek defterle karşılaştırılır.
- Kapalı adyabatik gazın birinci yasası, negatif olay uzantısında
  `cv * m_gas * dT/dt = -p * dV/dt` ile denetlenir. Kütleyi veya sıkıştırma
  işini sıfırlamak ve uzantı bayrağını gizlemek ayrı assertion'larla yakalanır.
- Gerçek `solve_ivp` adımlarıyla sıfır-akış kökü aranır. Testin negatif
  deneme adımı kullandığı, tek bir durdurucu köke ulaştığı ve motor kapalıyken
  kabul edilen adımlarda `m_tank + C*p_manifold` değişmezini koruduğu doğrulanır.
  Bu kontrol, negatif denemeyi tamamen yasaklayan bir düzeltmeyi de reddeder.
- Tek hatlı bağlantı aynı bağımsız defter ve alan işaretiyle ayrıca sınanır.

Bu dört yeni test yöntemi mevcut besleme kontrolleriyle birlikte çalıştırılır.
Üç ek seçilmiş kusur enjeksiyonu tank kütlesini kırpmayı, negatif ullage
işini kırpmayı ve olay bayrağını gizlemeyi hedefler. İlişkisiz istisnalar
başarılı kusur yakalama sayılmaz. Bu kontroller, fiziksel belirsizlik,
kalibrasyon veya bütün giriş uzayında doğruluk kanıtı değildir.

Güncel kaynak bağlı sonuçlar `output/assurance/verification_summary.json`
ve `feed_mutations.json` içindedir. [Donanım açıkları](PROPULSION_ASSURANCE.md)
ayrı kalır.

## Yapılandırma ve kesme olayı girdi sınırları

Ek incelemede oran sınırının `True` olarak kabul edilmesi, eksik adın ham
`AttributeError` üretmesi ve eksik tank/hat kayıtlarının kurulabilmesi yeniden
üretildi. Bilinmeyen olay türü reddedilmiyor; tanınan minimum türünün metin
biçimi bile maksimum sınırı gibi hesaplanıyordu. `stops_burn="false"` metni de
kabul ediliyordu. Bunlar fiziksel kalibrasyon sorunu değil, API sözleşmesi ve
olay seçimi kusurlarıdır.

Yapılandırma artık sayısal ve metin sınırlarında tür denetimi yapar; bool
sayısal değer değildir. Tamsayı ve float ile verilen aynı sayısal
yapılandırmanın model kimliği eşitlenir. Tank/hat ve bağlantı kayıtları
kuruluşta doğrulanır. Olay türü yalnız kayıtlı enum veya onun tanınan metin
değeri olabilir; sınır pozitif sonlu sayı, durdurma bayrağı bool olmalıdır.
Aktif oran girdisi de sonlu pozitif sayı olarak denetlenir. `None` mevcut
inaktif-motor sözleşmesini korur.

Dört yeni yöntem bozuk ad/sayıları, bool sınırlarını, eksik alt kayıtları ve
olayın minimum/maksimum işaretini sınar. Eski kodun başarısızlıkları
`bipropellant_boundary_before_fix.xml` içindedir; alt kontrol hataları ayrı
test yöntemleri veya ayrı fiziksel kusurlar diye sayılmaz. İki ek seçilmiş
mutasyon bool sayıları kabul etmeyi ve olay türünü doğrulamadan kullanmayı
hedefler. Bu aşamadaki besleme/API kümesi on üç seçilmiş kusurdu; ilk sekiz ve
olay uzantısındaki sonraki üç kusur tarihsel alt kümelerdir.

SciPy olay araması adım içindeki işaret değişimine dayanır. Buradaki tek kök
deneyi, bir adımda olabilecek bütün çoklu kökleri yakalama garantisi vermez.

## Karışım oranı türevinde sayısal ölçek ve iptal

Ortak itki/debi katsayısı karışım oranından ve onun basınç türevinden
analitik olarak elenir. Eski uygulamada `1e-160` ortak katsayısında yaklaşık
`5.1724e-7 1/Pa` türev sıfır oluyordu; daha küçük katsayılar ham
`ZeroDivisionError` üretiyordu. Karşı örnek
`output/assurance/bipropellant_scale_screening.json` içinde tarihsel kaynak
hash'iyle tutulur. `mixture_before_fix.xml` mevcut kontrollerin bu kusuru ve
sonuç kayıtlarındaki sonlu-olmayan değer kabulünü yakaladığını gösterir.

`d(O/F)/dx = (O'F - OF') / F²` artık verilen binary float değerlerini tam
rasyonel sayılara çevirerek hesaplanır. Ara çarpımlar ve fark yuvarlanmaz;
yalnız son çıktı float'a çevrilir. Bu, kayıp ölçüm basamaklarını geri getirmez
ve fiziksel belirsizliği azaltmaz. Float çıktıya sığmayan büyüklük veya sıfıra
yuvarlanan tam sıfır olmayan sonuç açık `DomainError` verir; gerçek sıfır
türev kabul edilir. Motor/hat oranları da bu çıktı denetiminden geçer.
Manifest `mixture_ratio_numerics` politikasını içerir; kaynak hash'leri ayrıca
gereklidir.

Altı yeni test yöntemi şu bağımsız dayanakları kullanır:

- Gerçek bağlı modelde `1e-280`–`1` ortak ölçekleri, katsayısı önceden
  analitik olarak elenmiş sentetik doğrusal basınç haritasıyla karşılaştırılır.
- Türev çekirdeğinde iki kuvvetleriyle tam ölçeklenen faktörlü akış yolları
  kullanılır; `2^-1000`–`2^1000` aralığı bütün bağlı sistemin fiziksel veya
  sayısal çalışma alanı iddiası değildir.
- Birbirine çok yakın bağıl eğimler için üretim formülü kopyalanmaz:
  bileşik doğrusal-kesirli akış yolunun iki tarafı Decimal ile değerlendirilir.
  `1e-50` yol adımı ve 160/240 basamaklı iki hesap aynı float referansta
  anlaşmalıdır. Bu seçilmiş örneklerde bir yakınsama kontrolüdür; bütün giriş
  uzayında hata sınırı veya deneysel geçerleme değildir.
- Temsil edilemeyen sonuçlar, tam sıfır, sonlu olmayan bileşik Jacobian ve
  oran çıktıları; boyut, koordinat benzersizliği ve değişebilir girdi
  matrisinden sonuç kaydına veri sızması ayrıca sınanır.

Dört yeni seçilmiş mutasyon ara çarpımları yuvarlamayı, bağıl eğimleri erken
yuvarlamayı, sıfıra kaybolan çıktıyı kabul etmeyi ve Jacobian denetimini
atlamayı hedefler. Bu aşamadaki besleme/API/sayısal küme on yedi seçilmiş kusurdu.
İlişkisiz istisna kusurun yakalandığı kanıtı sayılmaz.

Daha düşük ölçekli ek taramada, `1e-320` ortak katsayısında bileşen debi
türevlerinin köprüde önceden sıfıra yuvarlandığı ayrıca görüldü
(`mixture_upstream_subnormal_screening.json`). Tam-rasyonel O/F çekirdeği
kendisine ulaşmadan kaybolan bu bilgiyi geri getiremez. Sonraki
[köprü debi incelemesi](BRIDGE_FLOW_NUMERICS_REVIEW.md) bu temsil/çıktı kaybını
ele alır: tam döndürülen bileşen türevi float'a sığmıyorsa çağrı açık alan
hatasıyla kapanır; birleşerek temsil edilebilir toplam oluşturan terimler
korunur. Bu yeni davranış bütün pozitif float ölçeklerinde başarılı dönüş
garantisi veya fiziksel geçerleme iddiası değildir.

## References

- [Python 3.12 Fraction: float değerinin tam rasyonel karşılığı](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal: hassasiyet bağlamı ve float dönüşümü](https://docs.python.org/3.12/library/decimal.html)

- [SciPy 1.14.1 solve_ivp: events, terminal ve direction](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.solve_ivp.html)
- [NASA SP-8112, Pressurization Systems for Liquid Rockets](https://ntrs.nasa.gov/citations/19760015212)
- [NASA liquid-rocket feedline dynamics](https://ntrs.nasa.gov/citations/19740028545)
- [NASA-STD-7009B: doğrulama, geçerleme ve kanıt kapsamı](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
