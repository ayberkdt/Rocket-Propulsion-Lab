# Görev kütle kapanışının sayısal incelemesi

Görev, motor ve besleme tasarımı henüz seçilmedi. Aşağıdaki uç girdiler bir
donanım önerisi değildir; ideal roket denkleminin sayısal uygulamasını sınar.
Deney geçerlemesi veya TRL kanıtı sağlamaz.

## Yeniden üretilen kusurlar

Önceki hesapta artık oranı `nextafter(1,0)` ve
`delta_v=(1-f)*0.5*g0*200` iken 10 kg kuru kütle için yaklaşık 5 kg yakıt
bulunuyordu. Bağımsız ileri yayılım ve kök araması yaklaşık 10 kg gerektiriyor:
yakıt hatası %50. Buna rağmen eski toplam kütle kapanış artığı yaklaşık
`8.9e-16 kg` idi. Büyük toplam kütlenin kapanması küçük tüketimin doğru
hesaplandığını kanıtlamıyor. Sınırın iki tarafındaki komşu float girdiler de
ayrıca sınanır.

İkinci kusur, son çıktılar sonlu olmasına rağmen `exp(S)` ara kütle oranının
taşmasıydı: `S=1000`, kuru kütle `1e-300 kg`, sıfır rezerv ve artık oranında
yüklenen kütle yaklaşık `1.97e134 kg` ile temsil edilebilir. Bu örneğin fiziksel
bir uzay aracıyla ilişkisi yoktur. Üçüncü kusur, pozitif en küçük delta-v'nin
üs hesabında sıfıra yuvarlanıp sıfır tüketim/zaman olarak kabul edilmesiydi.
Artık temsil edilemeyen pozitif fiziksel çıktı açık `DomainError` üretir.

Düzeltme öncesi başarısızlıklar `output/assurance/mass_before_fix.xml`;
düzeltme sonrası hedefli koşu `mass_after_fix.xml`, son paket koşusu
`regression/test_sass_mass_boundaries.xml` içinde tutulur.

## Hesap yöntemi ve kapsamı

`S=sum(delta_v/(g0*Isp*eta))`, `B=exp(-S)`, kuru kütle D, yüklenen yakıt P,
artık oranı f ve sabit rezerv R için denklem:

`P = [D*(1-B)+R]/(B-f)`.

Bu, NASA Glenn ideal roket denkleminden burada türetilmiş cebirsel düzenlemedir.
Artık oranı, rezerv ve sayısal kabul politikası NASA tarafından verilmemiştir.
`1-exp(-S)` küçük S için yakınsayan Taylor serisiyle hesaplanır. Yalnızca
raporlanan son fiziksel kanallar binary float'a çevrilir; pozitif kanalın
sıfıra veya sonlu olmayan değere yuvarlanması reddedilir.

Hesap, çağıranın Decimal bağlamından ayrılmış 80, 160, 320, 640 basamaklı
bağlamları kullanır. İki ardışık çözülmüş hassasiyette bütün fiziksel binary
kanallar aynı olmalıdır. Her değerlendirmede toplam ve yakıt kütle defterleri
ayrı denetlenir. Payda işareti yuvarlama koruması içinde çözülemiyorsa daha
yüksek hassasiyet denenir; 640 basamakta çözülemeyen durum reddedilir.
On basamaklık koruma payı uygulama politikasıdır; kesin aralık ispatı veya
evrensel hata garantisi değildir.

JSON ve HTML şu ayrımı korur:

- `mass_closure_residual` ve `propellant_closure_residual`, raporlanan binary
  değerlerden `fsum` ile hesaplanır.
- `numerics` içindeki Decimal artıklar iç hesap yuvarlamasını kaydeder.
  Denenen hassasiyetler, payda çıkarımının koşul sayısı ve binary dönüşüm
  farkı ayrıca tutulur.
- Sayısal yakınsama, girdilerdeki ölçüm belirsizliğini gidermez. Özellikle
  uygulanabilirlik sınırı yakınında küçük girdi değişimleri büyük sonuç
  değişimleri yaratabilir. Fiziksel belirsizlik hâlâ hesaplanmamıştır.

## Bağımsız kontroller

`tests/test_sass_mass_boundaries.py` 115 basamaklı ileri roket-denklemi
yayılımı ve braketli kök araması kullanır; üretimdeki kapalı P formülünü
kopyalamaz. Yakıt sınırı, üç artık oranında uygulanabilirliğin iki tarafı,
taşan ara oran, pozitif çıktı alt taşması, küçük üs/büyük kuru kütle,
çağıranın Decimal bağlamı ve raporlanan defter artıkları sınanır. Çok küçük
üs örneğinin oracle'ı, ihmal edilen Taylor terimleri için analitik sınırdır.
Önceki 250 sabit tohumlu senaryonun 75 basamaklı geri yayılım oracle'ı da
korunur. Bunlar bütün float girdiler için doğruluk ispatı değildir.

Seçilmiş anlamsal mutasyonlar, hatalı son kütle raporunu iki hassasiyette
tutarlı biçimde üretse bile bağımsız kuru kütle kontrolüyle yakalanmalıdır.
Mutasyon sonuçları da kaynak/test/araç dosya kümesine bağlanır; değişmiş
girdilerle eski yakalanma kaydı yeniden başarılı gösterilemez.

## References

- [NASA Glenn ideal rocket equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/)
- [NASA Glenn specific impulse](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [Donanım ve geçerleme açıkları](PROPULSION_ASSURANCE.md)
