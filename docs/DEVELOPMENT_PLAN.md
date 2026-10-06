# Geliştirme Planı

Bu belge, Rocket Propulsion Lab'in `0.3.0` sürümünden doğrulanabilir bir
`1.0.0` mühendislik çalışma aracına ilerlemesi için uygulanacak sırayı,
bağımlılıkları ve kabul ölçütlerini tanımlar. Ayrıntılı fizik açıkları
[`GAP_ANALYSIS.md`](GAP_ANALYSIS.md), mevcut mimari kuralları ise
[`ARCHITECTURE.md`](ARCHITECTURE.md) içinde tutulur.

## 1. Ürün hedefi ve kapsam sınırı

`1.0.0` hedefi; roket itki sistemi ön boyutlandırmasında aşağıdaki analiz
zincirini tek bir çalışma dosyası içinde kurabilmektir:

```text
Akışkan / karışım seçimi
        ↓
Termokimyasal oda durumu
        ↓
Boğaz ve tasarım dışı nozül akışı
        ↓
Nozül geometrisi ve kayıplar
        ↓
İtki, Isp, kütle debisi ve görev karşılaştırması
        ↓
Grafik, tablo, belirsizlik ve dışa aktarma
```

Araç ön tasarım ve eğitim amaçlı olacaktır. Üç boyutlu CFD, yanma
kararsızlığı, ayrıntılı enjektör tasarımı, yapısal sonlu eleman analizi ve uçuş
sertifikasyonu `1.0.0` kapsamına dahil değildir. Bu sınır, "eksiksiz" ifadesini
belirsiz bırakmamak için bilinçli olarak konmuştur.

## 2. Geliştirme ilkeleri

1. Mevcut sabit-`gamma` çözücüler korunacak; yüksek doğruluk modelleri onların
   yerini sessizce almayacak, açıkça seçilebilir olacaktır.
2. Fizik kodu saf Python fonksiyonları ve değişmez `dataclass` sonuçları
   üretecek; API ve tarayıcı katmanı denklem içermeyecektir.
3. Her sonuç model adı, varsayımlar, birimler, uyarılar ve veri kaynağı
   bilgisini taşıyacaktır.
4. SI birimleri iç model olacaktır. Arayüzdeki birim dönüşümleri giriş/çıkış
   sınırında yapılacaktır.
5. Yeni bir fizik özelliği; referans değer testi, sınır testi, ters ilişki testi
   ve en az bir uçtan uca API testi olmadan tamamlanmış sayılmayacaktır.
6. Ağır termokimya bağımlılıkları isteğe bağlı tutulacak; ideal-gaz çekirdeği
   kurulum gerektirmeden çalışmaya devam edecektir.
7. Çalışma dosyaları sürümlü bir şemaya sahip olacak ve eski kayıtlar için
   göç kodu bulunacaktır.

## 3. Önceliklendirme

| Öncelik | İş paketi | Neden önce geliyor |
|---|---|---|
| P0 | Sayısal ve veri sözleşmesi temeli | Sonraki bütün modeller aynı birim, hata ve sonuç yapısına dayanır. |
| P0 | Sıcaklığa bağlı tür/karışım özellikleri | Roket oda sıcaklıklarında sabit `cp` ve `gamma` ana doğruluk sınırıdır. |
| P0 | Denge kimyası ve CEA doğrulaması | Gerçekçi `Tc`, molekül ağırlığı, `gamma`, `c*` ve Isp için gereklidir. |
| P0 | Tasarım dışı nozül rejimleri | Deniz seviyesi/irtifa davranışı ve iç şoklar mevcut nozül modelinin en önemli açığıdır. |
| P1 | MOC/Rao geometrisi | Geometriyi yalnızca görsel yaklaşımdan hesaplanmış kontura taşır. |
| P1 | Kayıplar, ısı transferi ve iki faz | İdeal performansı gerçek motor tahminine yaklaştırır. |
| P1 | Çalışma yönetimi ve belirsizlik | Mühendislik kararlarını tekrarlanabilir ve karşılaştırılabilir kılar. |
| P2 | Motor çevrimi ve görev optimizasyonu | Sağlam oda/nozül altyapısından sonra güvenilir anlam kazanır. |

## 4. Sürüm planı

Efor aralıkları tek geliştirici için odaklı çalışma günü tahminidir; takvim
taahhüdü değildir. Her faz bir öncekinin kabul kapısını geçtikten sonra başlar.

### Faz 0 — `0.3.1`: mühendislik temeli (4–6 gün)

**Durum:** Uygulandı. Ortak birim/meta-veri/hata sözleşmeleri, yakınsama
kanıtları, kaynak etiketli referans vakalar ve Python 3.11/3.12 CI matrisi eklendi.

**Teslimatlar**

- `core.units`: boyut adı, SI dönüşümü ve ekranda gösterim hassasiyeti.
- Ortak `ModelMetadata`, `WarningMessage` ve `ValidityRange` veri tipleri.
- JSON hata sözleşmesi: kod, alan, kullanıcı mesajı ve teknik ayrıntı.
- Sayısal çözücülerde yakınsama raporu, iterasyon sınırı ve açık hata türleri.
- Referans test verisini koddan ayıran `tests/reference_data/` yapısı.
- CI kapısı: birim/API testleri, Ruff ve paket kurulum denemesi.

**Kabul ölçütleri**

- Mevcut tüm hesaplar aynı sonuçları tolerans içinde üretir.
- Her API yanıtı model/varsayım/birim meta verisi taşır.
- Alan hataları HTTP 400, çözücü yakınsama hataları HTTP 422 olarak ayrılır.
- Temiz Python 3.11 ve 3.12 ortamlarında test ve paket kurulumu geçer.

### Faz 1 — `0.4.0`: değişken özellikli termodinamik (8–12 gün)

**Durum:** Uygulandı. NASA Glenn verili 10 tür, donmuş ideal karışımlar,
entalpi/entropi ters çözücüleri, API ve sabit özellik karşılaştırmalı üç grafik
eklendi. Kimyasal denge bilerek Faz 2 sınırında tutuldu.

**Teslimatlar**

- `thermochemistry/species.py`: tür kimliği, mol kütlesi ve sıcaklık aralığı.
- `thermochemistry/nasa_polynomials.py`: NASA 7 katsayılı `cp(T)`, `h(T)` ve
  `s(T)` değerlendirmesi.
- Başlangıç veri kümesi: H2, O2, H2O, OH, H, O, N2, CO, CO2, CH4 ve temel
  hidrokarbon yanma türleri; her kayıt kaynak ve aralık bilgisiyle.
- `Mixture` modeli: mol/kütle kesri dönüşümü, karışım mol kütlesi, `R`, `cp`,
  `cv`, `gamma`, entalpi ve entropi.
- Entalpi ve entropiden sıcaklık bulan sağlam ters çözücüler.
- Kalorik mükemmel gaz ile termal mükemmel gazı karşılaştıran UI bölümü ve
  `cp(T)`, `h(T)`, `gamma(T)` grafikleri.

**Kabul ölçütleri**

- Saf tür özellikleri yayımlanmış NASA polinom değerlerini seçili sıcaklıklarda
  tanımlı toleransla yeniden üretir.
- Kesirlerin toplamı, negatif oranlar ve polinom aralığı ihlalleri açık hata
  veya uyarı üretir; sessiz ekstrapolasyon yapılmaz.
- `h -> T -> h` ve `s -> T -> s` gidiş-dönüş testleri tüm veri aralıklarında
  geçer.
- Sabit özellikli eski API geriye dönük uyumlu kalır.

### Faz 2 — `0.5.0`: yanma dengesi ve CEA uyumluluğu (12–18 gün)

**Durum:** Çekirdek uygulandı; dış CEA kabul kapısı açık. TP/HP/UV sağlayıcı
sözleşmesi, deneysel yerel Gibbs çözücüsü, CEA subprocess adaptörü, O/F taraması
ve UI tamamlandı. Depoya CEA ikili dosyası gömülmediği için dokuz noktalı harici
CEA karşılaştırması kurulum sahibi ortamda çalıştırılmayı bekliyor.

**Teslimatlar**

- `thermochemistry/equilibrium.py`: element dengeli Gibbs serbest enerjisi
  minimizasyonu için sağlayıcı arayüzü.
- Önce NASA CEA ile dosya/subprocess adaptörü; ardından aynı arayüze uyan yerel
  çözücü için ayrı deneysel sağlayıcı.
- Sabit basınç/entalpi yanma, sabit iç enerji/hacim ve verilen `T-p` denge
  problemleri.
- Donmuş ve denge halinde nozül genişlemesi karşılaştırması.
- O/F taraması: adyabatik alev sıcaklığı, karakteristik hız `c*`, ortalama mol
  kütlesi, `gamma`, itki katsayısı ve vakum Isp grafikleri.
- CEA girdi/çıktı dosyaları için ayrıştırıcı, önbellek ve kaynak izi.

**Kabul ölçütleri**

- LOX/LH2, LOX/RP-1 ve LOX/CH4 için en az dokuz referans noktasında oda
  sıcaklığı, `c*` ve Isp sonuçları CEA ile önceden ilan edilmiş toleranslarda
  uyuşur.
- Element ve kütle dengesi kalıntıları her sonuçta raporlanır ve toleransın
  altında kalır.
- CEA kurulu değilse ideal çekirdek çalışır; kullanıcıya kurulabilir özellik
  eksikliği anlaşılır biçimde bildirilir.
- Adaptör çıktısı ham metne değil sürümlü Python veri modellerine dönüşür.

### Faz 3 — `0.6.0`: tasarım dışı sıkıştırılabilir nozül akışı (10–15 gün)

**Durum:** Sabit özellikli çekirdek uygulandı. Rejim eşikleri, iç normal şok,
irtifa taraması, Summerfield ayrılma ekranı, artıklar, API ve grafikler
tamamlandı. Değişken özellik/denge sağlayıcısının aynı rejim çözücüsüne doğrudan
bağlanması sonraki entegrasyon kapısında tutuluyor.

**Teslimatlar**

- Geri basınca göre rejim sınıflandırması: boğulmamış, yeni boğulmuş, nozül
  içinde normal şok, çıkışta şok, aşırı genişlemiş, tasarım noktası ve eksik
  genişlemiş.
- Verilen geri basınç için iç şok konumu ve şok sonrası subsonik çıkış çözümü.
- İrtifa/ortam basıncı taraması ve itki–irtifa, Isp–irtifa, basınç oranı
  grafikleri.
- Ampirik ayrılma başlangıcı uyarısı; korelasyon adı ve geçerlilik aralığıyla.
- Değişken özellik ve donmuş/denge genişleme sağlayıcılarını kullanabilen ortak
  nozül arayüzü.
- Şok, ayrılma uyarısı ve akış rejimini nozül konturu üzerinde gösteren grafik.

**Kabul ölçütleri**

- Rejim sınırları analitik limitlerde süreklidir ve her rejim için referans test
  vardır.
- İç şok konumunda kütle debisi korunur; toplam basınç kaybı ve çıkış basıncı
  artıklarla birlikte raporlanır.
- Geri basınç azaldıkça rejimlerin fiziksel sırası bozulmaz (özellik tabanlı
  test).
- UI, modelin ayrılma tahmini yaptığını fakat ayrıntılı viskoz çözüm olmadığını
  açıkça belirtir.

### Faz 4 — `0.7.0`: hesaplanmış nozül geometrisi (12–18 gün)

**Durum:** Uygulandı. Konik ve `preliminary` Hermite sözleşmesine Rao/TOP ve
düzlemsel uyumluluklı, eksenel-simetrik alan eşlemeli MOC eklendi. Ağ inceltme,
sınır artıkları, ortak karşılaştırma, karakteristik ağ grafiği ve CSV/SVG/DXF/STL
dışa aktarma testlerle kapatıldı.

**Teslimatlar**

- `geometry/moc.py`: düzlemsel ve eksenel simetri düzeltmeli karakteristik ağ.
- Minimum uzunluklu nozül üretimi, karakteristik kesişim kontrolleri ve ağ
  yakınsama raporu.
- Rao tabanlı parabolik bell seçeneği; boğaz/çıkış açıları ve yüzde uzunluk
  parametreleri.
- Mevcut Hermite bell modelinin `preliminary` olarak yeniden adlandırılması.
- Kontur karşılaştırması: konik, preliminary bell, Rao ve MOC için uzunluk,
  yüzey alanı, tahmini diverjans kaybı ve çıkış düzgünlüğü.
- CSV ve SVG dışa aktarma; doğrulama tamamlandıktan sonra DXF, ardından STL.

**Kabul ölçütleri**

- MOC duvar sınırında akış teğetliği ve merkez hattı simetrisi kalıntıları
  raporlanır.
- Ağ inceltmeyle çıkış Mach sayısı ve kontur uzunluğu kararlı değere yakınsar.
- DXF/STL çıktısı yeniden içe alındığında boğaz ve çıkış alanları tolerans
  içinde korunur.
- Aynı tasarım girdisinde bütün kontur türleri ortak bir `NozzleContour`
  sözleşmesi üretir.

### Faz 5 — `0.8.0`: kayıplar ve ısıl yükler (15–25 gün)

**Durum:** Uygulandı (ön tasarım/deneysel sınırlar açık). Açılıp kapanabilir
kayıp bütçesi, belirsizlik bandı, Bartz duvar dağılımı, tek fazlı deneysel
rejeneratif kanal dengesi ve yoğuşmuş faz bağlaşım bandı API ve UI ile eklendi.
Eksik soğutucu verisi varsayılmıyor; alan ihlalleri hata veya kaynak etiketli
uyarı üretiyor.

**Teslimatlar**

- Deşarj, diverjans, sınır tabakası ve yanma verimi kayıp bütçesi.
- Bartz temelli gaz tarafı ısı akısı; korelasyon ve özellik seçiminin açık
  gösterimi.
- Duvar boyunca adyabatik duvar sıcaklığı, ısı akısı ve toplam ısı yükü.
- Basit rejeneratif kanal enerji dengesi ve soğutucu basınç kaybı için ayrı,
  deneysel modül.
- Yoğuşmuş faz kütle kesri ve parçacık sürükleme kaybı için önce uyarı/bant
  modeli, doğrulama verisi geldikçe genişleyen iki faz sağlayıcısı.
- İdeal, kayıplı ve belirsizlik bantlı performansın yan yana gösterimi.

**Kabul ölçütleri**

- Her kayıp açılıp kapatılabilir ve toplam farkın hangi terimden geldiği izlenir.
- Korelasyonlar kaynak, birim, geçerlilik aralığı ve sınır dışı uyarısı taşır.
- Enerji dengesi ve toplam kayıp bütçesi artıklarla doğrulanır.
- Termal hesaplar, eksik malzeme/soğutucu verisini varsayarak gizlice devam
  etmez.

### Faz 6 — `0.9.0`: çalışma alanı ve sistem analizi (12–20 gün)

**Durum:** Uygulandı. Sürümlü çalışma şeması ve v1→v2 göçü, vaka
karşılaştırması, sabit tohumlu LHS/Monte Carlo, duyarlılık sıralaması, kısıtlı
kademe optimizasyonu, üç ön çevrim enerji dengesi ve JSON/CSV/SVG/HTML rapor
paketi API ve tarayıcı iş akışına bağlandı.

**Teslimatlar**

- Sürümlü `.rplab.json` çalışma dosyası; kaydetme, açma, çoğaltma ve göç.
- Birden fazla vakayı aynı grafikte karşılaştırma ve sonuç fark tablosu.
- Parametre taraması, Latin hypercube/Monte Carlo belirsizlik yayılımı ve
  duyarlılık sıralaması.
- Çok kademeli delta-v, yapısal oran, faydalı yük oranı ve kademe optimizasyonu.
- Basınç beslemeli, gas-generator ve staged-combustion çevrimleri için önce
  bileşen enerji dengesi; pompa/türbin haritası sonraki alt sürüme bırakılır.
- CSV/JSON/SVG rapor paketi ve tüm girdileri/varsayımları içeren yazdırılabilir
  HTML özet.

**Kabul ölçütleri**

- Kaydedilen çalışma tekrar açıldığında aynı sürümde sayısal olarak aynı sonucu
  üretir.
- Sabit tohumlu belirsizlik analizi tekrarlanabilir ve örnekleme bilgisi sonuçta
  bulunur.
- Optimizasyon, kısıt ihlallerini çözüm olarak sunmaz; yakınsama ve aktif
  kısıtları raporlar.
- Dışa aktarılan rapor model sürümlerini, veri kaynaklarını ve bütün uyarıları
  içerir.

### Faz 7 — `1.0.0`: doğrulama ve kararlı yayın (8–12 gün)

**Durum:** Yerel yayın kapıları uygulandı; kararlı sürüm kapısı açık. Tolerans
matrisi, kullanıcı/API kılavuzu, örnek çalışma, klavye/print/küçük ekran CSS'i,
Windows paket betiği ve 3.11/3.12 CI hazırdır. `1.0.0` etiketi, harici CEA
kurulumunda dokuz referans noktasının kapatılmasına kadar bilinçli olarak
verilmemiştir; mevcut paket sürümü `0.9.0`'dır.

**Teslimatlar**

- Altın referans vakalar: ideal gaz, CEA karşılaştırması, tasarım dışı nozül,
  MOC kontur ve kayıplı performans.
- Sayısal tolerans matrisi ve model bazlı doğrulama raporu.
- API şeması, Python kullanım örnekleri ve kullanıcı kılavuzu.
- Erişilebilirlik, klavye kullanımı, küçük ekran ve yazdırma denetimi.
- Windows paketleme, örnek çalışma dosyaları ve sürüm notları.

**Kabul ölçütleri**

- Bütün P0/P1 özellikler doğrulama tablosunda kaynak ve toleransla kapalıdır.
- Kritik fizik modüllerindeki dallar referans veya sınır testleriyle kapsanır.
- Temiz makinede kurulum, ilk hesap, dosya kaydetme/açma ve dışa aktarma uçtan
  uca geçer.
- Bilinen sınırlamalar arayüzde ve belgede aynı ifadeyle yayımlanır.

### Sonraki genişleme — bağımsız burn analizi ve Sidera bağlantısı

Yörünge mekaniğini bu projede yeniden kurmadan; bağımsız Python/API/CLI/UI
kullanımı, propellant envanteri, toplam impuls, değişken performans profilleri,
motor kümeleri ve state-dependent propulsion sağlayıcılarını kapsayan burn
katmanı planlanmıştır. Sidera mevcut sabit finite-burn modeli için opsiyonel
adaptördür; burn çekirdeğinin çalışması Sidera kurulumuna bağlı değildir.

Ayrıntılı sahiplik sınırı, veri modelleri, denklemler, sayısal yöntemler,
Sidera uyumluluk matrisi, fazlar ve kabul kapıları
[`SIDERA_BURN_INTEGRATION_PLAN.md`](SIDERA_BURN_INTEGRATION_PLAN.md) içinde
tanımlanmıştır. Sürüm numarası, bu çalışma başlamadan önce dış CEA `1.0`
kapısının kapanıp kapanmamasına göre `0.x` veya sonraki `1.x` feature minor
olarak seçilecektir.

## 5. Mimari evrim

Hedef paket yapısı aşağıdaki gibi olacaktır:

```text
src/rocket_propulsion/
├── api/                 HTTP sözleşmesi ve şema dönüştürme
├── compressible/        1B akış, şoklar, kanal ve nozül rejimleri
├── core/                birimler, meta veri, uyarılar, nümerik çözücüler
├── geometry/            conical, preliminary bell, Rao, MOC, dışa aktarma
├── propulsion/          itki, görev ve motor çevrimleri
├── studies/             vaka, tarama, belirsizlik ve optimizasyon
├── thermochemistry/     tür verisi, karışımlar, denge ve CEA adaptörü
├── thermodynamics/      süreçler, enerji dengeleri ve ısı transferi
└── web/                 sunum, giriş doğrulama ve grafikler
```

Katman bağımlılığı tek yönlü kalır: `web -> api -> domain -> core`.
`thermochemistry` bir sağlayıcı protokolü sunar; `compressible` veya
`propulsion` belirli bir CEA kurulumunu doğrudan içe aktarmaz. Grafikler için
Python örneklenmiş seri ve meta veri üretir, JavaScript yalnızca etkileşim ve
çizim yapar.

## 6. Test ve doğrulama stratejisi

Her iş paketi aşağıdaki beş kanıt türünü sağlamalıdır:

| Kanıt | Amaç | Örnek |
|---|---|---|
| Denklem testi | Tek ilişkinin doğruluğu | NASA polinomunda `cp(T)` |
| Limit testi | Fiziksel sınır davranışı | `M -> 1`, `p_b -> p_e` |
| Ters/özellik testi | Geniş giriş uzayında tutarlılık | `A/A* -> M -> A/A*` |
| Referans vaka | Harici güvenilir araç/veriyle uyum | CEA LOX/LH2 vakası |
| Uçtan uca test | UI/API/model zinciri | çalışma aç, hesapla, dışa aktar |

Toleranslar test içine gelişigüzel gömülmeyecek; referans kaynağı, mutlak/bağıl
tolerans ve gerekçesiyle veri dosyasında tutulacaktır. Rastlantısal testlerde
tohum sabitlenecek, başarısız örnek yeniden üretilebilir olacaktır.

## 7. Veri ve kaynak yönetimi

- Tür katsayıları makinece okunur veri dosyalarında tutulur; lisans, yayın,
  sürüm ve sıcaklık aralığı her kayıtla birlikte gelir.
- CEA karşılaştırma dosyaları girdi, ham çıktı, ayrıştırılmış sonuç ve kullanılan
  CEA sürümünü birlikte saklar.
- Ön tanımlı itici çiftleri yalnızca ad değil; varsayılan O/F aralığı, faz,
  basınç aralığı ve kaynak izi taşır.
- Kaynak değiştirilirse referans sonuç farkı otomatik üretilmeden veri güncellemesi
  kabul edilmez.
- Kullanıcının kaydettiği çalışma dosyasına büyük ham tablolar gömülmez; veri
  kümesi kimliği ve sürümü kaydedilir.

## 8. UI ve grafik planı

Arayüz tek uzun hesap makinesi sayfasından çalışma odaklı düzene evrilecektir:

1. **Girdiler:** akışkan, oda, nozül, ortam ve model seçimi.
2. **Özet:** itki, Isp, `c*`, debi, rejim ve kritik uyarılar.
3. **Detay:** termokimya, istasyon tablosu, kayıp bütçesi ve yakınsama bilgisi.
4. **Grafikler:** ortak yakınlaştırma, imleç değeri, seri aç/kapat, birim seçimi
   ve SVG/CSV indirme.
5. **Karşılaştırma:** seçilen vakaları aynı eksenlerde ve fark tablosunda gösterme.

Her grafikte eksen adı ve birimi zorunlu olacak; kullanılan model, uyarı ve
geçerlilik sınırı grafiğin yanında görülecektir. Kritik rejim değişimleri çizgi
rengiyle değil, etiket/işaretçiyle de aktarılacaktır.

## 9. İlk üç sprintlik uygulanabilir sıra

### Sprint 1 — sözleşme ve altyapı

- Ortak sonuç meta verisi ve uyarı tipleri.
- Birim dönüşüm katmanı.
- API hata şeması ve regresyon testleri.
- CI ve referans veri klasörü.

**Sprint çıkışı:** `0.3.1`; mevcut özelliklerde davranış kaybı olmadan yeni
fizik modüllerinin kullanacağı kararlı sözleşme.

### Sprint 2 — saf tür özellikleri

- NASA polinom ayrıştırma ve değerlendirme.
- İlk tür veri kümesi ve kaynak doğrulaması.
- `cp/h/s/gamma` eğrileri, aralık uyarıları ve ters sıcaklık çözücüleri.

**Sprint çıkışı:** değişken özellikli saf gaz hesapları, API ve grafiklerle
uçtan uca çalışır.

### Sprint 3 — karışım ve nozül entegrasyonu

- Mol/kütle kesirli `Mixture` modeli.
- Karışım entalpi/entropi ve adyabatik süreç hesabı.
- Mevcut nozül çözücüsüne sabit veya değişken özellik sağlayıcısı seçimi.
- Sabit-`gamma` ile değişken özellik sonuçlarının karşılaştırma görünümü.

**Sprint çıkışı:** `0.4.0`; denge kimyasına geçmeden önce veri ve sağlayıcı
arayüzü gerçek bir akış zincirinde doğrulanmış olur.

## 10. Karar kapıları ve başlıca riskler

| Karar / risk | Kapı | Azaltma yaklaşımı |
|---|---|---|
| CEA dağıtımı ve lisanslama | Faz 2 başlamadan | Çalıştırılabilir dosyayı depoya gömmek yerine kullanıcı kurulumu ve adaptör kullanmak. |
| Yerel denge çözücü kapsamı | CEA adaptörü doğrulandıktan sonra | Önce sağlayıcı sözleşmesi; yerel çözücüyü deneysel etiketlemek. |
| Tür verisi kaynağı ve sürümü | Faz 1 veri eklenmeden | Kaynak/provenans şemasını önce tamamlamak. |
| Ayrılma korelasyonu seçimi | Faz 3 UI uyarısından önce | Birden fazla korelasyonu sağlayıcı olarak tutmak, geçerlilik alanını göstermek. |
| MOC eksenel düzeltme doğruluğu | DXF/STL yayınından önce | Ağ yakınsaması ve yayımlanmış kontur vakası olmadan CAD dışa aktarmayı kararlı saymamak. |
| Kapsamın CFD'ye kayması | Her faz gözden geçirmesi | `1.0` sınırını korumak; yüksek doğruluk ihtiyacını harici çözücü adaptörüne yönlendirmek. |

## 11. Tamamlanma tanımı

Bir iş maddesi ancak kod, belge ve arayüz birlikte güncellendiğinde; bütün yeni
girdiler doğrulandığında; kaynak ve model sınırı gösterildiğinde; birim,
referans, sınır ve API testleri geçtiğinde; eski çalışma/API davranışı ya
korunduğunda ya da belgeli göç yolu sağlandığında tamamlanır.

`1.0.0` ise yalnızca özelliklerin eklenmesiyle değil, Faz 7 doğrulama raporunda
her P0/P1 kabiliyetin dış referans, tolerans ve bilinen sınırlamayla izlenebilir
olmasıyla tamamlanmış sayılır.
