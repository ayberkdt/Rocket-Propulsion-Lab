# Yön ölçeği, kuvvet ve kestirim türevleri

İnceleme, `PropagatorPropulsionBridge`, `ThrottleableEngineModel` ve
`PowerLimitedPropulsionModel` içindeki ortak yön geometrisini kapsar.
Donanım yöneltme doğruluğu, fiziksel performans veya bütün motor modellerinin
doğruluğu için tamamlanmış geçerleme değildir. Test girdileri sentetiktir.

## Yeniden üretilen hatalar

Eski kareler toplamı, `(1e308, 0, 0)` yönünde taşarak 150 N kuvveti sıfır
vektörüne dönüştürüyordu. `(1e-300, 0, 0)` ise sıfır yön diye reddediliyordu.
Her iki durum üç tüketicide de yeniden üretildi. Ayrıca `(1, 1e-300, 0)` ve
`1e300 N` için `1-u_x²` çıkarımı, temsil edilebilir eksenel yön türevini
sıfırlıyordu. Önceki başarısızlıklar `direction_before_fix.xml` içindedir.

## Hesap ve kapsam

Yeni özel yardımcı yönü en büyük bileşene böler ve sınırlı büyüklüklerde
`hypot` kullanır. Norm tahmini ölçekle rasyonel olarak birleştirilir; normun
kendisinin float aralığını aşması tek başına ret sebebi değildir.

Orijinal koordinatlar ve kareler toplamı verilen binary değerler üzerinden
tam rasyonel tutulur. Kuvvet/ivme projeksiyonunda birim vektör ara çıktıya
yuvarlanmaz. Örneğin `(1e308, 1e-308, 0)` yönündeki temsil edilemeyen küçük
birim bileşen, `1e308` katsayısıyla çarpıldığında tekrar temsil edilebilir.
Sonuç alanı hâlâ temsil edilemiyorsa sessiz sıfır yerine alan adı içeren
`DomainError` oluşur. Tam sıfır geçerlidir.

Yön Jacobian'ında `T/(m |d|) (I-ddᵀ/(dᵀd))` kullanılır; köşegen farkları
yuvarlanmış birim bileşenler üzerinden çıkarılmaz. `T` işaretlidir; negatif
basınç kökü deneme uzantısında türevin işareti korunur. Kütle, komut, güç,
parametre, çalışma koşulu ve ortam basıncı ivme türevleri aynı projeksiyonu
kullanır. Manifestlerdeki `direction_numerics` alanı model hash'ini değiştirir.

Bu yöntem tam karekök veya tüm girdilerde doğru son bit iddiası taşımaz:
ölçekli `hypot` ve ölçekli koordinat bölmeleri float tahminleridir. Profil
kuvvet örneklemesi, elektrikli çalışma noktası ve debi eğrileri gibi önceki
hesap katmanları burada bütünüyle yeniden doğrulanmış değildir.

## Bağımsız kontroller

Dokuz yeni yöntem üç tüketicide ölçek değişmezliğini, `2^±900` ölçeklerini,
float üst sınırını aşan normu, eksene yakın yönü, işaretli itkileri, kapalı
durumu, kurtarılan bileşeni ve gerçek çıktı aralığı hatasını sınar.
Referans Jacobian, üretimdeki projektör formülünü kopyalamadan Decimal kuvvet
yolunun iki taraflı farkıyla hesaplanır. 160/240 basamak sonuçları karşılaştırılır;
zor iptal örneklerinde 700/900 veya 800/1000 basamak gerekir. Referansın
float'a dönüşmüş sonuçları iki hassasiyette uyuşmadan karşılaştırma yapılmaz.
Bu, örneklenmiş sayısal kanıttır; evrensel hata sınırı değildir.

İlk tam regresyonda 520/521 yöntem geçti. Kalan entegrasyon kontrolü enine
hızı `sqrt(|v|²-v_parallel²)` ile ölçerken yuvarlama yüzünden negatif kök
girdisi üretti. Kontrol `|v × u|` normuna geçirildi; `1e-9` kabul sınırı
aynıdır. Başarısız koşunun kaynakları ve sonuçları
`history/direction_first_full_run/archive_manifest.json` altında korunur.

Beş seçilmiş kusur enjeksiyonu norm ölçeğinin atlanmasını, birim bileşenin
erken yuvarlanmasını, köşegen çıkarımının geri getirilmesini ve itkinin mutlak
değerinin türevde kullanılmasını ve entegre enine ivme sapmasını hedefler.
Sonuncusu yenilenen entegrasyon kontrolünün fiziksel yön hatasını yakaladığını
ayrıca sınar. Yalnız assertion başarısızlığı kusur
yakalama sayılır; ilişkisiz istisna sayılmaz. Güncel kaynak bağlı sonuçlar
`output/assurance/verification_summary.json` içindedir.

[Donanım olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz kapı açıktır.

## References

- [Python 3.12 hypot](https://docs.python.org/3.12/library/math.html#math.hypot)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
