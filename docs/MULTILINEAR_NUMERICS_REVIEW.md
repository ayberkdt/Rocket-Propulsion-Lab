# Çok doğrusal harita ve kanal defteri incelemesi

İnceleme `RectilinearPerformanceSurface` çok doğrusal yolunu ve
`RectilinearStreamPerformanceSurface` bağımsız akışkan kanallarını birlikte
kapsar. Hücre geometrisi, katsayılar, gradyanlar ve köprüye aktarılan
kuvvet/tank tüketimi denetlenir. Veriler sentetiktir; fiziksel motor haritası,
ölçüm belirsizliği veya motor/besleme donanımı geçerlenmiş değildir.

## Yeniden üretilen kusurlar

Sonlu eksen uçlarının farkı float'a sığmayınca hücre ortası ilk düğüm gibi
hesaplanıyor ve temsil edilebilir gradyan sıfırlanıyordu. Bir başka örnekte
`1e-200 × 1e-200` köşe ağırlığı sıfır oluyor; `1e300` düğüm katsayısıyla
birleşince kurtulabilecek `1e-100` katkı kayboluyordu. Sabit büyük alanların
küçük aralıktaki gradyanı da sıfır yerine `NaN` olabiliyordu.

Ek kontroller temsil edilemeyen normalize hücre kesrinin son subnormal
çıktıyı değiştirebildiğini ve büyük karşıt gradyan terimleri arasında kalan
`0.5` katkının eksen sırasına bağlı olarak kaybolduğunu gösterdi.
Gerçek çıktı aralığı ihlalleri açık alan hatası vermiyordu.

Bağımsız akışkan haritası aynı kaybı motor köprüsüne taşıyordu: sentetik
örnekte beklenen `-2e-100 kg/s` yerine `-2e-300 kg/s` tüketim görülüyordu.
Bu, fiziksel görev önerisi veya ölçülmüş motor sonucu değildir.

Güncel 14 yöntem, arşivlenmiş eski paket kaynaklarıyla yeniden çalıştırıldı:
10 yöntem başarısız, dört kontrol yöntemi başarılı oldu. Paket kökenleri,
kaynak baytları ve test dosyası hash'i kontrol edildi. Alt kontroller ayrı
test yöntemi sayılmadı. Kayıt: `multilinear_historical_replay.json` ve `.xml`.

## Hesap ve çıktı sözleşmesi

Hücre genişliği ve kesri, verilen binary değerler üzerinden rasyonel tutulur.
Ortak özel çekirdek köşe ağırlığını ve türev ağırlıklarını üretir. Türev
ağırlığı, türevi alınan eksenin kesrini çarpıma katmaz; o eksendeki işaret ve
genişliğin tersini kullanır. Bu, düğüm ve iç hücre gradyanlarını ayırır.

Kuvvet, toplam debi ve adlandırılmış akışkan kanalları bağımsız rasyonel
defterlerden toplanır. Ağırlık veya katkı tek başına float'a dönüştürülmez.
Böylece küçük köşe katkıları ve büyük terimlerin iptalinden kalan küçük
gradyanlar korunur. Her çıktı alanı sonunda bir kez dönüştürülür.
Temsil edilemeyen sıfır olmayan gradyan alan adıyla `DomainError` verir;
gerçek sıfır gradyan geçerlidir. Sonlu sonuca rağmen ara taşma gerekçesiyle
bütün girdi reddedilmez.

Mevcut uç nokta toleransı, kabul edilen koşulu hücreye yerleştiren sayısal
bir kuraldır; bu toleransın görev/ölçüm gerekçesini veya deneyle geçerli
çalışma alanını onaylamaz. Birlik ağırlık toplamı ve pozitif düğümlerle kurulan
çok doğrusal değerlerin matematiksel sınırlılığı, fiziksel model doğruluğu değildir.

İki yüzey belgesine de `interpolation_numerics` politikası bağlandı.
Çok doğrusal yüzey hash'leri değişir; okuyucular bilinen politika alanının
silinmesini/değiştirilmesini reddeder. Kübik yol kendi ayrı politikasını
kullanır. Sayısal kimlik denetimi ölçüm kaynağını veya yazarı doğrulamaz.

## Bağımsız kontroller

Referans, üretimdeki köşe ağırlığı döngüsünü kopyalamaz. Hücreyi bağımsız
seçer ve ardışık bir boyutlu Decimal doğrusal harmanlamayla değeri kurar.
Gradyanlar seçilmiş hücrenin polinom uzantısında iki taraflı farkla bulunur;
düğümde sağ hücre türevi bu nedenle sol türevle karıştırılmaz.
1200/1600 basamaklı referansların float sonuçları uyuşmadan karşılaştırma
yapılmaz. Bu örneklenmiş sayısal kanıttır; evrensel hata sınırı değildir.

14 yöntem genişlik/ağırlık/kesir kaybını, sabit alanı, gradyan iptalini ve
ekseni yeniden sıralamayı; sabit tohumlu bir/iki/üç boyutlu bağımsız kanalları;
düğüm türevlerini; gerçek aralık hatalarını; kimlik okuyucularını ve tam binary
afin koordinat dönüşümünü kapsar. Köprü kontrolleri kuvveti, tank/araç
tüketimini ve koşul türevlerini birlikte karşılaştırır. Değişken koşullu
yanma, analitik tüketim integrali ve kütle kapanışıyla ayrıca denetlenir.

Altı seçilmiş ek kusur ağırlık/kesrin erken yuvarlanmasını, iki gradyan
defterinin float'a geri dönmesini, türev ekseni kesrinin yanlışlıkla çarpıma
katılmasını ve kuvvet türevinin akışkan kanalına kopyalanmasını hedefler.
Yalnız assertion başarısızlığı yakalama sayılır; ilişkisiz istisna sayılmaz.
Güncel kaynak bağlı tam sonuçlar `output/assurance/verification_summary.json`
içindedir. Seçilmiş mutasyonlar bütün olası kusurları taramaz.

[Donanım olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz kapı açıktır.

## References

- [SciPy 1.14.1 RegularGridInterpolator](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.interpolate.RegularGridInterpolator.html)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
- [NASA transient-modeling overview](https://ntrs.nasa.gov/citations/20040000363)
