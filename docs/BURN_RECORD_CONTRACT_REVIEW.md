# Yanma kayıtları ve artifact sınırı

Bu inceleme `models.py` içindeki 14 kayıt türünün kurucularını, enumları,
ilgili iç doğrulayıcıları ve v1/v2 yanma sonucu JSON giriş/çıkışını kapsar.
Sınırlı sayısal karşı örnekler ayrıca çalışma noktası Isp/karışım oranı ve
başlangıç kütle bütçesi denetimine uygulanır. Bütün fizik modelleri veya bütün
değişim biçimleri bu incelemeyle doğrulanmış sayılmaz.

## Değişen davranış

| Önceki kusur | Güncel sözleşme |
|---|---|
| Boolean değerler kütle, hedef veya sayaç olabiliyordu | Fiziksel büyüklükler sonlu gerçek sayılara, sayaçlar boolean olmayan tam sayılara normalize edilir; hatalar alan bağlamı taşır |
| Rol/hedef/hedef-kesme/destinasyon metinleri kimlik karşılaştırmalarını bozabiliyordu | Bilinen enum yazımları enum değerine çevrilir; tanımsız değerler reddedilir |
| Dondurulmuş kayıt içindeki listeler sonradan değişebiliyordu | Akış, tank, uyarı, profil ve isim/kütle çiftleri tuple olarak kopyalanır |
| JSON `"false"` değeri `True`, kesirli sayaç `1` olabiliyordu | Okuyucu ham türleri kuruculara iletir; `str/float/int/bool` ile anlam değiştiren ön dönüşümler kaldırılır |
| Üç elemanlı çiftin son elemanı veya tanımsız alanlar hash dışında kalabiliyordu | Çift biçimi ve benzersiz kimlikler sıkı denetlenir; her kayıt tam olarak şemanın alanlarını taşır |
| Yinelenen JSON alanının yalnız son değeri kullanılıyordu | Her derinlikte yinelenen ad reddedilir; NaN/Infinity tokenları kabul edilmez |
| Yeniden hesaplanmış hash, değiştirilmiş sonucu dışa aktarmaya yetiyordu | Yazıcı da okuyucu gibi tür, hash ve ilan edilen girdilerden yeniden üretim denetimi yapar |
| Ara çarpım/oran/toplam taşması geçerli Isp'yi reddediyor veya yanlış kütle/karışım onayı veriyordu | Verilen binary girdilerin rasyonel oran ve toplamları karar bitene kadar korunur |

Örneğin geçerli bir artifact içindeki `target_achieved: true` alanını
`target_achieved: "false"` yapmak eski okuyucuda aynı geçerli hash ile kabul
edilebiliyordu. Yeni okuyucu bunun boolean olmadığını bildirir. Tanımlı
sayısal alanın `1` veya `1.0` olarak verilmesi aynı fiziksel değerin tek
normalize kimliğini üretir; `true` veya `"1"` bu sayısal alan için kabul edilmez.
Serbest Unicode metni ve emoji korunur; kanonik UTF-8'e kodlanamayan
surrogate karakteri alanı belirtilmiş hata üretir.

Kurucular boyutsal alanların temel aralıklarını, türlerini ve kayıt yapısını
denetler. Sonuçların her çapraz bağıntısı kurucuda yeniden çözülmez: çözücü,
kendi hash'i hesaplanmadan önce geçici `result_hash=""` kaydı oluşturabilir.
Bu geçici kayıt yazıcıdan çıkamaz. Sonucu elle değiştirip hash'i yeniden
hesaplamak da modelden yeniden üretim denetimini atlatmaz. Hash, gönderenin
kimliğini veya ölçümün doğruluğunu onaylayan imza değildir.

Sayaçlara keyfî bir binary64 üst sınırı konmaz. Sıfır adım/segment kanıtı,
sıfır süreli gerçek kesme, boş uyarı listesi, sıfır debili izlenmeyen akış ve
boş izlenmeyen tank kimliği için geçerli kontroller korunur. Fiziksel skaler
çıktılar binary64 olduğundan bu alanların normalize değerleri sonlu olmalıdır.

Çalışma noktasının önceki akış/ISP/karışım kapanış eşikleri ve tanımın kütle
bütçesi eşiği bu turda sıkılaştırılarak karşı örnekler gizlenmedi. Aynı
eşik ifadeleri, ara taşma ve `inf > inf` yanlış kararından korunur. Bu
eşikler deneysel belirsizlik bütçesi veya donanım kabul ölçütü sayılmaz;
L0/L1'in toleransla kabul edilmiş akış uyuşmazlığını ayrı raporlaması korunur.

## Kapsam ve karşı kanıt

`test_burn_record_contracts.py` içindeki 25 yöntem; gerçek çözücü sonuçları,
geçerli v1/v2 alışverişleri, değişebilir liste saldırıları, ham tür/alan
saldırıları, tekrar eden JSON adları, sahte yeniden hash'lenmiş sonuçlar,
Unicode ve tipik/sıfır/uç ölçek kontrollerini içerir. Isp örneğinde bağımsız
Decimal100/180 referansı, sonlu oranı bozabilecek büyük ara çarpıma rağmen
aynı kabul edilebilir sonucu verir. Kütle bütçesi bağımsız Decimal500 ile
kontrol edilir. Geçersiz fizik girdisini yalnız beklenmedik bir hata üreterek
reddetmek test başarısı sayılmaz.

Önceki doğrudan envanter testleri, yeni kurucuların artık reddettiği bozuk
kayıtları yardımcı sınırına ulaştırabilmek için açıkça kopyada alan bozarak
çalışır. Yardımcıdaki savunmalar veya assertion'lar kaldırılmadı.

Son test dosyası arşivlenmiş eski kaynak üzerinde tekrar çalıştırılır;
paketin gerçekten o kaynak kökünden yüklendiği ve kaynak/test hash'leri
doğrulanır. İlk keşif test dosyası ile nihai dosya ayrı tutulur. Seçilmiş
14 kayıt/artifact kusuru ve bir envanter-kapsam kusuru bellek içinde
enjekte edilir. Yalnız assertion başarısızlığı geçerli yakalamadır; yanlış
enjeksiyon konumu veya ilgisiz çalışma zamanı hatası sayılmaz.

Yazıcıdaki tekrar çözüm bir maliyettir: dışa aktarma artık salt `json.dumps`
değildir. Bu doğrulama çözücünün kendi hash hesabına eklenmez; döngüsel
yeniden çözüm oluşmaz. Aynı model ve girdilerden tekrar sonuç elde etmek,
bağımsız bir fizik modeline veya deneye karşı geçerleme değildir.

## İnceleme envanterinin onarılması

Bağımsız `ast.walk` taraması önceki koleksiyonun 26 dosyada üç iç fonksiyon
ve 44 lambda'yı atladığını gösterdi. Kontrol akışı gövdeleri, imzalar ve
bütün AST çocukları artık ziyaret edilir. Lambda adında statik kapsayıcı,
satır ve sütun bulunur; bu ad çalışma zamanı `__qualname__` iddiası değildir.
`test_source_review_inventory.py` içindeki iki yöntem hem iç içe bir saldırı
örneğini hem bütün kaynak ağacındaki AST konum kümesini karşılaştırır.
Toplam envanter sayısının artması daha fazla incelemenin tamamlandığı anlamına
gelmez; yeni görünür öğeler ayrıca kapsam etiketi taşır.

Bu liste açık Python class/function/async-function/lambda bildirimlerini
kapsar; kaynak ağacı `src/rocket_propulsion` ile sınırlıdır. Dekoratörün
ürettiği yöntemler, dinamik callable'lar, bu ağaç dışındaki script'ler, JavaScript,
katalog verileri ve teknik iddialar için eksiksiz inceleme kanıtı değildir.
Yeni sayılar eski eksik sayılarla doğrudan ilerleme yüzdesi olarak karşılaştırılmaz.

Güncel kaynak bağlı sonuç `output/assurance/verification_summary.json` ve
inceleme durumu `objective_audit.json` içindedir. Ayrı serbest uzay
kuadratür/tabelalı moment aritmetiğinde yeni uç ölçek karşı örnekleri
`free_space_metrics_followup_screening.json` ile kaydedildi; bu turda
onarılmış sayılmaz. Dokuz donanım kanıtı kapısı açık kaldığından TRL 5
sağlandığı veya sertifika alınabileceği sonucu çıkarılmaz.

## Shared normalization follow-up

[CORE_VALIDATION_CONTRACT_REVIEW.md](CORE_VALIDATION_CONTRACT_REVIEW.md) records
the later shared scalar/text/container review. Its branch audit strengthened
the non-depleting-stream Unicode case, corrected manual literal normalization,
and then retired the now redundant model-local Unicode deletion variant from
the active curated set. Earlier checkpoint results remain historical evidence;
the retired variant is not reported as a new kill.

## References

- [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259): bölüm 4, 6, 8.1–8.2. Tekil ad için SHOULD ile projenin daha sıkı reddetme politikası ayrı tutulur.
- [Python 3.12 JSON](https://docs.python.org/3.12/library/json.html): object_pairs_hook, parse_constant ve varsayılan dönüştürme davranışları.
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html): verilen float değerinin tam rasyonel karşılığı.
- [Donanım kanıtı sınırları](PROPULSION_ASSURANCE.md)
