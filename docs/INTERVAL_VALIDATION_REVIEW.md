# Deney karşılaştırmasının aralık ve kabul incelemesi

Bu modül, sağlanan sayısal aralıkları karşılaştırır. Deneyin gerçekliğini,
sensör kalibrasyonunu, modelin dosyadan yeniden çalıştırılmasını veya donanım
olgunluğunu doğrulamaz. Onaylı fiziksel kabul sınırları henüz mevcut değildir.

## Yeniden üretilen kusur

Eski kabul kuralındaki sayısal koruma, incelenen sınır örneğinde hatalı bir
geçişi engelliyordu. Ancak `maximum_error` adıyla raporlanan binary sayı,
tam aritmetikteki gerçek üst hata sınırından küçük olabiliyordu. Benzer
biçimde alt sınır da içeri doğru yuvarlanabiliyordu. Başka bir tüketicinin
bu sayıları güvenilir aralık uçları olarak kullanması doğru değildi.

Örnek: nominal fark `3.0941676495717276`, iki sınır
`5.050647017735236` ve `7.246650282258111` iken eski üst hata
`15.391464949565073` raporlanıyordu. Tam toplam daha büyüktür. Bunlar
sentetik yazılım girdileridir; gerçek motor ölçümü değildir.

Başarısız ilk koşu `output/assurance/interval_before_fix.xml`, sonraki
hedefli koşu `interval_after_fix.xml` içinde saklanır. Alt kontrol
başarısızlıklarının sayısı ayrı fiziksel kusur veya test metodu sayılmaz.

## Aralık ve sayısal karar politikası

Ölçülen m, öngörülen p, mutlak ölçüm sınırı Um ve mutlak öngörü sınırı Up
için tam aralık:

`lower = max(0, abs(p-m)-Um-Up)`

`upper = abs(p-m)+Um+Up`.

Girdiler binary sayıların tam rasyonel değerleriyle hesaplanır. Raporlanan
alt uç aşağı, üst uç yukarı yuvarlanır; böylece tam aralık içeri daraltılmaz.
Üst sınır için sonlu bir dışa yuvarlanmış sayı yoksa `DomainError` üretilir.
Nominal `discrepancy` değeri ise en yakın binary gösterime yuvarlanan bir
tanıdır; dışa yuvarlanmış uç olduğu iddia edilmez.

Sayısal inceleme koruması, raporlanan üst hata ve izin verilen hata
ölçeğinin büyüğünün 16 float epsilon ile çarpımının yukarı yuvarlanmış
değeridir. Bu, bizim sayısal politika seçimimizdir; NASA'nın belirlediği
donanım toleransı değildir. Bildirilen sınır ve koruma değerleri üzerinden
karşılaştırmalar tam rasyonel aritmetikte yapılır:

- `maximum_error + numerical_guard < allowable_error`: sayısal karşılaştırma pass.
- `minimum_error > allowable_error + numerical_guard`: sayısal karşılaştırma fail.
- Diğer durumlar: indeterminate.

Çok küçük subnormal bir koruma, bir binary nicem kadar yukarı yuvarlanabilir.
Bu durumda sınıflama daha temkinli olur; gizli sıfır koruma raporlanmaz.
Eksik sınır veya dışarıda kalan doğrulama alanı unassessed kalır. Bir başarısız
holdout örneği toplam karşılaştırma sonucuna baskın gelir.

## Bağımsız kontroller

Oracle ölçüm ve öngörü aralıklarının dört köşesindeki imzalı farklardan en
büyük mutlak hatayı ve sıfırı içerme durumundan en küçük hatayı hesaplar.
Üretimdeki merkez farkı/artı-eksi toplam formülünü kopyalamaz. 256 sabit
tohumlu örnekte raporlanan uçların bu tam köşe aralığını kapsaması, pass'in
tam üst hatayı izin içinde bırakması ve fail'in tam alt hatayı izin dışında
bırakması denetlenir. Üç karar durumunun da oluşması aranır; sürekli
indeterminate dönen bir uygulama testi geçemez. Bu örnekleme fiziksel bir
olasılık dağılımı veya bütün girdiler için doğruluk ispatı değildir.

Ek sınamalar sonlu float üst sınırının aşılması, subnormal koruma ve ilan
edilmiş inceleme payının kaldırılmasını kapsar. İki ek mutasyon üst ucu
içeri yuvarlamayı ve sayısal korumayı kaldırmayı enjekte eder. Yalnızca
assertion başarısızlığı yakalanmış kusur sayılır.

## Fiziksel geçerleme için kalanlar

Um/Up bu arayüzde mutlak hata sınırlarıdır; standart sapma veya güven
aralığı otomatik olarak bunların yerine kullanılamaz. Kapsam ve kabul
kuralı proje tarafından gerekçelendirilip onaylanmalıdır. Şablon bu kaydı
boş bırakır. Dosya özeti yalnızca incelenen byte'ları bağlar; sayısal
transkripsiyonu veya onay yetkisini doğrulamaz. Model kimliği ve özeti
sağlanan beyan olarak kaydedilir, bağımsız olarak denetlenmiş sayılmaz.

Sonuçta `empirical_validation_status=not_established_by_this_numeric_check`
ve `qualification_evidence=false` korunur. Görev, donanım, ilgili ortam,
metroloji ve bağımsız fiziksel holdout kanıtı seçilmeden TRL kararı verilmez.

## References

- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [NASA hardware TRL definitions](https://esto.nasa.gov/trl/)
- [BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure)
- [Açık sistem kanıt kapıları](PROPULSION_ASSURANCE.md)
