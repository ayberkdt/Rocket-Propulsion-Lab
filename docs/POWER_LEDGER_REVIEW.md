# Faz enerjisi hesabının sayısal incelemesi

Bu inceleme verilen sabit faz güçlerinin ve sürelerinin sayısal defterini
kapsar. EPS, batarya, PPU veya termal donanım geçerlemesi sağlamaz.

## Yeniden üretilen iki kusur

Çalışma oranı `d=nextafter(1,0)` ve ateşleme süresi 1 s olduğunda bekleme
süresi `(1-d)/d` olur. Eski `t/d-t` hesabı iki yuvarlanmış sayıyı çıkardığı
için bu süreyi iki kat gösterebiliyordu. Bekleme gücü `1e20 W` seçilen yapay
örnekte yaklaşık 11102.23 J yerine 22204.46 J raporlanıyordu. Büyük güç
değeri gerçek tasarım önerisi değildir; küçük fazdaki hatayı görünür yapar.

İkinci kusurda ateşleme ve ısınma gücü en küçük pozitif binary float,
iki fazın süresi de 0.4 s idi. İki çarpım ayrı ayrı sıfıra yuvarlansa da
toplamın doğru yuvarlanmış değeri pozitif en küçük float'tır. Eski uygulama
sıfır enerji ve ortalama güç raporluyordu.

Düzeltme öncesi başarısızlıklar `output/assurance/power_before_fix.xml`;
düzeltme sonrası hedefli koşu `power_after_fix.xml` içinde kayıtlıdır.

## Hesap ve bağımsız kontroller

`power_numerics.solve_power_budget`, sonlu binary girdileri tam rasyonel
sayılar olarak alır. Bekleme süresini `t*(1-d)/d` ile; enerjiyi ayrık
ateşleme, bekleme ve ısınma fazlarının toplamıyla hesaplar. Yalnızca dört
rapor kanalı binary sayıya çevrilir: enerji, toplam süre, en yüksek güç ve
ortalama güç. Pozitif fiziksel çıktının sıfıra veya sonlu olmayan değere
yuvarlanması açık `DomainError` üretir. Sıfır süreli faz tepe güce eklenmez;
pozitif bekleme fazı, yuvarlanmış toplam süre ateşleme süresinden farklı
görünmese de hesaba katılır.

Son çıktı dönüşüm denetimi depolama hesaplarıyla ortak
`rational_numerics.binary_ledger` yardımcısında uygulanır.

Bu aritmetik politikası belirsiz girdileri kesin fiziksel değerler yapmaz.
Faz güçleri zaten tüm bara yüklerini içermeli; ısınma ateşlemeden önce ve
ayrı olmalıdır. Gölge, güneş paneli, dönüştürücü verimi ve geçiş enerjileri
bu modelde çözülmez.

`tests/test_sass_power_boundaries.py`, 115 basamaklı Decimal enerji hesabıyla
dört zaman ölçeğini karşılaştırır; küçük enerji toplamını da doğru binary
yuvarlama sınırıyla denetler. Bekleme süresinde eski çıkarım ve fazları
toplamadan önce eski yuvarlama ayrıca bellekte enjekte edilir. Her iki
kusur assertion başarısızlığı üretmelidir; ilişkisiz istisna başarılı
yakalama sayılmaz. Bu testler deney toleransı veya bütün girdiler için
geçerleme değildir.

## References

- [BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [Donanım ve geçerleme açıkları](PROPULSION_ASSURANCE.md)
