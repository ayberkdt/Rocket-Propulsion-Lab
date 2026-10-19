# Köprüde kütle debisi ve türev defteri

İnceleme, `PropagatorPropulsionBridge.evaluate` içindeki debi örneklemesi,
ölçekleme, tanklara gruplama ve debinin parametre/koşul türevlerini kapsar.
Kuvvet, yön normalizasyonu ve bütün köprü davranışları için tamamlanmış
doğrulama iddiası değildir. Bütün test verileri sentetiktir; motor veya
besleme donanımı için deneysel geçerleme sağlanmamıştır.

## Yeniden üretilen kusurlar ve düzeltme

Önceki taramada `1e-320` ortak katsayısında sıfır olmayan bileşen debi
türevleri oran hesabına ulaşmadan sıfırlanıyordu. Ek testlerde tankların aynı
parametre/koşul türevini farklı çarpım sıralarıyla hesapladığı görüldü:
sonlu sonuca rağmen ara çarpım taşabiliyor veya sıfıra kaybolabiliyordu.
Kanal ve tank toplamları da küçük bir katkıyı büyük karşıt terimler arasında
kaybediyordu. İlk altı yeni test yönteminin eski koddaki başarısızlıkları
`output/assurance/bridge_flow_before_fix.xml` içinde tutulur; alt kontrol
hataları ayrı yöntem veya ayrı fiziksel kusurlar diye sayılmaz.

Köprü, debi ve debi eğimini zaman/debi uç değerlerinden tam rasyonel
aritmetikle örnekler. Böylece tanklar birleştirilmeden önce her kanalın ayrı
float'a yuvarlanması önlenir. Sabit tutma politikasında debi sabittir ve eğimi
sıfırdır. Doğrusal politikada uç noktalar ve süre farkı üzerinden hesaplanır.
Kanal akışları varsa araç toplamı ve türevleri aynı kanal defterinden kurulur;
kanal yoksa toplam debi kanalı kullanılır.

Harita/sürücü çarpımları, tank grupları ve türev toplamları yalnız ihraç
edilen her alanın sonunda float'a çevrilir. Tam sıfır olmayan bir alan float
olarak temsil edilemiyorsa açık `DomainError` oluşur. Bu, bütün küçük
girdilerin reddedilmesi değildir: ayrı ayrı temsil edilemeyen terimler aynı
tankta temsil edilebilir bir toplam oluşturabilir veya tam sıfıra iptal
olabilir. Testler bu durumların başarılı kalmasını da ister.

`numeric_outputs.py` aynı çıktı dönüşümünü çift yakıtlı oran çekirdeğiyle
paylaşır. Tamlık, verilen binary girdilere göredir; ölçümde kaybolmuş bilgiyi
geri getirmez. Manifest `mass_flow_numerics` politikasını içerir; model hash'i
bu nedenle değişir. Kaynak kodu hash'leri ayrıca tutulur.

## Oracle'lar ve kusur enjeksiyonları

Yedi yeni test yöntemi:

- Temsil edilemeyen bileşen türevinde sessiz sıfır yerine açık alan hatası.
- En küçük float'ın yarısı olan iki terimin toplamda kurtarılması; ayrı
  tanklarda temsil edilemeyen çıktının reddi ve tam iptalin kabulü.
- Küçük tank tüketimi ve debi sürücüsü türevinin aynı defterle korunması.
- Karşıt büyük kanallar yanında küçük kanalın analitik olarak bilinen
  `1 kg/s` koşul türevinin altı kanal sıralamasında korunması.
- Sabit toplamlı büyük profil kanallarının yanında `1+t` küçük kanalının
  bilinen ateşleme-zamanı türevinin sıralamadan bağımsız kalması.
- Ara taşma/alt taşmadan kurtulan sonlu sonucun 160/240 basamaklı Decimal
  çarpım hesabıyla karşılaştırılması.
- Örneklemeden önce birleşmesi gereken küçük kanallar ve bitişteki sol limitte
  tam son değer; sağ taraftaki kapanış davranışı değiştirilmez.

Beş yeni seçilmiş mutasyon erken örnek/terim yuvarlamasını, float tank
biriktirmesini ve sıfıra kaybolan çıktının kabulünü hedefler. Ortak yardımcının
mutasyonu bütün yüklü köprü/besleme bağlarına uygulanır; yalnız bir takma adı
değiştirip yanlış bir deney sonucu üretilmez. Bu inceleme aşamasındaki
besleme/köprü/API kümesi 22 seçilmiş kusurdu; güncel sayı aşağıdaki kaynak
bağlı rapordadır. Bu bir eksiksiz kusur kapsamı veya TRL kararı değildir.

Güncel kaynak bağlı sonuçlar `output/assurance/verification_summary.json`
içindedir. [Olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) donanım kapıları
açık kalır.

## References

- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [NASA transient-modeling overview](https://ntrs.nasa.gov/citations/20040000363)
