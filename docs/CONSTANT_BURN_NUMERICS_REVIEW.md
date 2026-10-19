# Sabit yanma hedefi, kesme sınırları ve tüketim defteri

Bu inceleme L0 `simulate_constant_burn` yolunu kapsar. L1 profil çözücüsünün
hedef/zaman terslemesi, eski bağımsız envanter yardımcıları ve donanım modeli
bu incelemeyle doğrulanmış sayılmaz. Bütün saldırı girdileri sentetiktir;
çok küçük hedefler gerçek bir motorun en küçük darbe kabiliyetini göstermez.

## Yeniden üretilen kusurlar

Eski çözücü küçük hedefte `1-exp(-x)` çıkarımı nedeniyle süreyi yanlış
hesaplıyordu: `1e-12 m/s` isteğinde yaklaşık `8.27225e-13 m/s` teslim edip
başarı bildiriyor; `1e-14 m/s` isteğinde `ZeroDivisionError` veriyordu.
Pozitif tüketim, yuvarlanmış ilk/son tank durumlarının çıkarımında ayrıca
sıfırlanıyordu. Sıfır kullanılabilir envanter de Isp hesabında sıfıra bölünüyordu.

Hedef toleransı, maksimum süreyi ve korunan envanter sınırını az aşan bir
isteğin süresini doğrudan uygulayabiliyordu. Yeni karar seçilmiş kesme
süresini ileri taşımaz; eşitlikte hedefin mevcut önceliği korunur.

Ek testler küçük tüketim/kütle oranının float'ta sıfırlanmasını, temsil
edilemeyen ara egzoz hızından doğan sonsuz çıktıyı ve neredeyse tükenmiş
depoda pozitif artığa rağmen logaritmanın tanım alanı dışında kalmasını
yeniden üretti. Gerçek çıktı aralık hataları alan bilgili hata vermiyordu.

## Hesap ve karar sözleşmesi

Debi, akış başına tank gruplaması, kütle/itki hedeflerinin süreleri, bütün
rasyonel kesme süreleri, tüketimler ve artıklar verilen binary girdilerle
rasyonel tutulur. Skaler çıktılar sonunda denetimli float'a dönüştürülür.
Tank tüketimi debi-zaman çarpımıdır; gösterilen iki kütlenin farkından
kurulmaz. Bu nedenle iki kütle aynı float'a yuvarlansa bile tüketim korunur.
Tam defterin kapanışı, yuvarlanmış alanların tam eşit olacağı iddiası değildir.

Δv terslemesindeki üstel ve ileri hesabın logaritması, çağıranın Decimal
bağlamından bağımsız 1100 basamaklı bağlamda hesaplanır. Üstelin rasyonel
sonucu 1'den float'a dönüştürülmeden çıkarılır. Bu hassasiyet, binary64
girdilerin bu denklemlerdeki üs yayılımını kapsamak için seçilmiştir;
evrensel doğru yuvarlama veya bütün çıktılar için kanıtlanmış hata zarfı
iddiası değildir. Üstel argümanı 4096'yı aşarsa kalan oran sıfırla temsil
edilir: en büyük binary64 ilk kütlede dahi gerçek kalan kütle en küçük
pozitif float'ın çok altındadır. Daha erken sonlu kesme sınırı yine çözülebilir.

Isp'ler süreye veya tüketimin sıfırdan farklı olmasına bölünmeden, sabit
performans oranından hesaplanır. Sıfır envanterde süre/tüketim/Δv sıfır,
hedef başarısız ve kesme nedeni tanımlıdır. Gerçek sıfır geçerlidir; pozitif
bir çıktının temsil aralığından sıfıra kaybolması veya taşması alan adıyla
`DomainError` verir. Süre ve kütleler dışa ayrı yuvarlanır; olayın tam süresi
sonuçtaki float zamandan daha hassas olabilir.

İşletim noktasının önceden izin verdiği küçük akış/ilan edilmiş toplam
farkı sessizce kapatılmaz; ayrı araç ve tank defterleri kütle kapanışı
artığında bunu gösterir. İtki kapanışı aynı sabit bağıntının öz tutarlılık
kontrolüdür ve bağımsız fiziksel kanıt değildir. Kayıttaki mevcut tolerans
alanları bütün girdilere ait ispatlanmış hata sınırı olarak kullanılmamalıdır.

Kanıt yöntem adı ve tarayıcıdaki yöntem kaydı
`analytic_constant_binary_ledgers_decimal1100_v2` olur. Sonuç hash'leri
değişir; girdi şeması değişmez. Hash ölçüm kaynağını doğrulamaz.

## Bağımsız kontroller ve kapsam sınırı

On altı test yöntemi küçük/büyük ölçekleri, komşu float kesme sınırlarını,
sıfır envanteri, pozitif artık/logaritmayı, gerçek çıktı aralığı hatalarını,
akış uyuşmazlığını, cant katsayısını, ortam Decimal hassasiyetini ve JSON
dönüşümünü kapsar. Alt kontroller ayrı yöntem sayılmaz. Kuru kütle sınırı
testi, mevcut tanımın kabul ettiği bir-ulp yük/kuru bütçe uyuşmazlığını da
korunan sınırı aşmak için kullanamadığını denetler; fiziksel konfigürasyon
önerisi değildir.

Ters referans doğrudan üretim üstel çıkarımını kopyalamaz:
`(1-exp(-x))/x` için integral serisini kullanır. İleri referans küçük
tüketimlerde logaritma yerine `atanh` tek kuvvet serisini kullanır.
200/300 basamaklı referansların aynı float'a yerleşmesi denetlenir;
karşılaştırma sınırı dört ulp'dir. Sınır ve tüketim kontrolleri ayrıca tam
bilinen binary değerler ister. Altı seçilmiş kusurun hepsi yalnız assertion
başarısızlığıyla yakalanmalıdır; ilişkisiz istisna başarılı deney sayılmaz.

Arşivlenmiş eski kaynakla yeniden oynatma `constant_historical_replay.json`
içindedir. Güncel tam kaynak bağlı sonuçlar
`output/assurance/verification_summary.json` içindedir. Bağımsız modellerin,
L1 profil yolunun veya fiziksel ölçümlerin doğruluğu bu sonuçlardan çıkmaz.
[Olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz donanım kapısı açıktır.

## References

- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
