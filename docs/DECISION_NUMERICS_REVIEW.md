# Karar puanlaması ve görev ölçeklemesinin sayısal incelemesi

Puanlar, kullanıcının seçtiği tercih ölçekleriyle yapılan deterministik ön
tasarım hesabıdır. Sayısal eşitliklerin doğru hesaplanması, bir motorun veya
besleme sisteminin fiziksel olarak uygun seçildiğini kanıtlamaz.

## Yeniden üretilen iki kusur

Üç kriterde aynı `float(1/3)` ağırlığı kullanıldığında `(0, 0.5, 1.0)` ve
`(0, 0.625, 0.875)` fayda değerleri aynı ağırlıklı toplamı verir. Eski kod
her ağırlıklı terimi ayrı binary sayıya yuvarladığı için ilk adaya `0.5`,
ikinciye `0.49999999999999994` puanı veriyordu. Böylece matematiksel olarak
eşit iki adaya yapay bir sıra oluşuyordu. Bu örneğin girdileri sıradan
büyüklüktedir. Ayrı bir alt taşma kontrolünde iki pozitif katkının toplamı
temsil edilebilirken terimlerin önce yuvarlanması sıfır puan üretiyordu.

İkinci kusurda pozitif bir görev Δv tahsisi ve pozitif ölçek katsayısının
çarpımı temsil dışına küçüldüğünde, kod bunu sıfır gereksinim olarak kabul
ediyordu. Yeniden üretim `1e-300 * 1e-30` gibi yapay uç ölçek kullanır;
bu bir gerçek görev önerisi değildir. Temsil edilemeyen pozitif değişim,
kullanıcının bilerek seçtiği sıfır tahsisle aynı anlama gelmez.

Başarısız ilk koşu `output/assurance/decision_before_fix.xml`; düzeltme
sonrası hedefli koşu `decision_after_fix.xml` içinde saklanır.

## Hesap ve kanıt zinciri

Sabit en iyi/en kötü çıpalardan elde edilen kırpılmış doğrusal fayda ve
ağırlıklı katkılar tam rasyonel aritmetikte tutulur. Son toplam bir kez
binary puana çevrilir. Görüntülenen fayda değerleri ayrıca yuvarlanır ve
hesabın sonraki girdisi yapılmaz. Puanın kaynak/varsayım soyu asıl metrik,
çıpa ve ağırlıklardan taşınır.

Sıfır veya eksik verinin yorumlanması, sert sayısal taramalar ve adayların
girdi sırası korunur. Ağırlıklar yeniden normalize edilmez; mevcut toplamın
1'e `1e-12` mutlak yakınlık kontrolü sayısal giriş politikasıdır, NASA'nın
belirlediği bir tercih toleransı değildir. Yakın fakat farklı gerçek puanlar
son binary gösterimde eşit görünebilir. Soyut tercih değerleri için son
nicemleme kalır; hiçbir puan farkı fiziksel anlamlılık sınırı olarak sunulmaz.

Görev ölçeklemesi, tam rasyonel çarpımı ortak `binary_ledger` sınırından
geçirir. Pozitif kütle/Δv değişimi sıfıra yuvarlanıyorsa `DomainError` verir.
Açık sıfır Δv katsayısı geçerli bir varsayımsal senaryodur. Başlangıç isteği
değiştirilmez, her senaryo aynı başlangıçtan yeniden hesaplanır.

## Bağımsız kontroller ve kalan karar açığı

Eşitlik testi 115 basamaklı Decimal ağırlıklı toplamıyla karşılaştırılır.
Küçük katkıların korunması, pozitif ölçeklemenin kaybolmaması ve açık sıfır
senaryonun isteği değiştirmemesi ayrıca sınanır. İki mutasyon, katkıların
toplanmadan önce yuvarlanmasını ve ölçekleme sınır denetiminin atlanmasını
yeniden enjekte eder; assertion başarısızlığı üretmeleri gerekir.

NASA SE 6.8, karar belirsizliklerinin değerlendirilmesini ve yakın sıralı
adayların karar vericiye birlikte sunulmasını ele alır. Burada fiziksel
girdi/model belirsizliği nicelenmemiştir; sıralamanın bu belirsizlik altında
dayanıklılığı **değerlendirilmedi**. JSON/HTML bu açığı açık tutar. Henüz
seçilmemiş görev ve donanım için nihai mimari veya TRL kararı verilmez.

## References

- [NASA SE Handbook 6.8: Decision Analysis](https://www.nasa.gov/reference/6-0-crosscutting-technical-management/)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [Donanım ve geçerleme açıkları](PROPULSION_ASSURANCE.md)
