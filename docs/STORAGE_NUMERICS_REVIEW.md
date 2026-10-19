# Depolama ve basınçlandırıcı hesabının sayısal incelemesi

Bu çalışma verilen yoğunluk ve ideal gaz varsayımı altında envanter hesabını
inceler. Tank dayanımı, gazın gerçek davranışı veya motor/besleme donanımının
ilgili ortamda geçerlemesi yapılmış değildir. Aşağıdaki uç ölçekler gerçek
akışkan özelliği ya da basınç/sıcaklık önerisi değildir.

## Yeniden üretilen kusurlar

| Kimlik | Kusur | Gözlenen etki |
|---|---|---|
| STORAGE-01 | `m*R*T` veya `R*T` ara işleminin binary taşması/alt taşması | Sonlu tank hacmi reddediliyor veya pozitif hacim sıfır raporlanıyor |
| STORAGE-02 | `p/(R*T)` önce yuvarlanıp yer değiştiren hacimle çarpılıyor | Aktarılan gaz ve yüklü/artık envanter sıfır olabiliyor; ters ölçekle sonlu sonuç gereksiz reddediliyor |
| STORAGE-03 | Pozitif son boş hacmin sıfıra yuvarlanması kabul ediliyor | Pozitif ullage, sıfır hacim gibi raporlanıyor |

Örneğin `m=1e-200 kg`, `R=1e-200 J/(kg*K)`, T=1 K,
başlangıç/son basınçlar `2e-200/1e-200 Pa` için ideal denklem pozitif ve
temsil edilebilir `1e-200 m^3` hacim verirken eski hesap sıfır veriyordu.
Bu yalnızca sayısal uygulamanın uç ölçek deneyidir; bu değerlere sahip bir
gerçek gazın var olduğu veya modelin bu koşullarda fiziksel olarak geçerli
olduğu iddia edilmez.

Başarısız ilk koşu `output/assurance/storage_before_fix.xml` içindedir.
Hedefli düzeltme sonrası koşu `storage_after_fix.xml`; güncel paket koşusu
`regression/test_sass_storage_boundaries.xml` içinde saklanır.

## Düzeltme ve sayısal kabul sınırı

Üç depolama fonksiyonu, doğrulanmış sonlu binary girdilerini tam rasyonel
sayılar olarak değerlendirir. Son fiziksel kanallar yalnızca bir kez binary
float'a dönüştürülür. Pozitif kanalın sıfıra veya sonlu olmayan değere
yuvarlanması `DomainError` üretir. Sıfır envanter, sıfır ek hacim/kütle verir;
bu dejenere hesap sonucu gerçek bir sıfır hacimli tank önerisi değildir.

Bu dönüşüm sınırı faz enerjisiyle ortak `rational_numerics.binary_ledger`
yardımcısında tutulur. Sonlu son çıktılar için gerekmeyen `R*T` ara temsil
kısıtı kaldırılmıştır. Basınç sıralaması, mutlak sıcaklık, pozitif R/yoğunluk
ve toplam tank hacmi üzerinden ullage tanımı korunur.

Tam iç aritmetik, ölçüm belirsizliğini gidermez; son binary değerlerde
yuvarlama vardır. Özellikle subnormal çıktılarda bağıl nicemleme büyük
olabilir. Başlangıç/son basınçlar birbirine yaklaşırken gerekli envanter
girdilere hassaslaşır. Bu fonksiyonlar fiziksel belirsizlik veya donanım
kabul aralığı hesaplamaz.

## Bağımsız test ve kusur enjeksiyonu

Oracle 115 basamaklı Decimal aritmetiğiyle başlangıç gaz envanteri ile son
gaz envanterinin farkını değerlendirir. Hacmi kapalı formülden kopyalamak
yerine braket kurup bisection ile çözer. Her uç envanter bu kökten ayrıca
hesaplanır; yalnızca toplamdaki görünür kapanışa bakılmaz.

Kontroller ara çarpım taşması/alt taşması, ara termal enerji taşması,
basınçlandırıcı yoğunluğunun iki yönde taşması, pozitif temsil dışı ullage,
sıfır envanter ve temsil dışı nihai hacim davranışlarını kapsar. Önceki
küçük artık gaz regresyonu da korunur. Bunlar bütün float girdiler için
fiziksel doğruluk ispatı değildir.

İki ek seçilmiş mutasyon, aktarımdan önce yoğunluğu binary yuvarlamayı ve
pozitif çıktının sıfıra yuvarlanmasını kabul etmeyi yeniden enjekte eder.
İkisi de assertion başarısızlığı üretmelidir; ilişkisiz istisna yakalanmış
kusur sayılmaz. Bu iki kusur SASS seçilmiş mutasyon kümesine eklenmiştir.

## Fiziksel model sınırları

İdeal gaz envanteri aynı sabit sıcaklıktaki uç koşulları kullanır. Gerçek
gaz özellikleri, blowdown soğuması, ısı transferi, başlangıç ullage gazı,
regülatör dropout/kayıp payı, çözünme ve hatlarda kalan gaz ayrıca gerekir.
Bir son basıncın regüle tank basıncından büyük olması, gerçek regülatörün
istenen debiyi sağlayabildiğini tek başına kanıtlamaz.

NASA SP-8112 tasarımın sistem gereksinimlerinden başlayıp donanım/entegrasyon
ve sistem değerlendirmesine ilerlediğini belirtir. Buradaki cebirsel hesabın
doğrulanması bu aşamaların veya donanım gösteriminin yerine geçmez.

## References

- [NASA Glenn ideal gas equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/)
- [NASA SP-8112: Pressurization systems for liquid rockets](https://ntrs.nasa.gov/citations/19760015212)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [Açık donanım ve geçerleme kapıları](PROPULSION_ASSURANCE.md)
