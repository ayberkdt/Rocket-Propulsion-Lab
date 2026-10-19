# L1 profil maruziyeti, zaman terslemesi ve tüketim incelemesi

Kapsam parça parça doğrusal komutun integralleri, en erken maruziyet-zaman
terslemesi, L1 kesme kararları ve profil sonuçlarının sayısal defterleridir.
İtki ve bütün akışların komutla orantılı ölçeklendiği fiziksel varsayım
henüz deneyle geçerlenmiş değildir. Girdiler sentetiktir; küçük hedefler gerçek
motorun en küçük darbe kabiliyetine kanıt oluşturmaz.

## Yeniden üretilen kusurlar

Eski L1 yolunda pozitif `1e-10 m/s` hedefi, sıfır süre/Δv ve başarılı hedef
kararı üretiyordu. Üstel çıkarımı ve zaman terslemesindeki mutlak tolerans
pozitif maruziyeti kaybediyordu. Uzun bir aralıktaki gerçek eğim de boyutlu
bir eşikle sıfır sayılabiliyordu. Toplam maruziyetin hemen üzerindeki istek
ulaşılabilir kabul ediliyor, bekleme aralığı sonrası en erken zaman yanlış
çözülebiliyordu.

Hedef toleransı maksimum süreyi veya takvim sonunu aşabiliyordu. Tank
tüketimi yuvarlanmış ilk/son kütlelerin farkında kayboluyor; rampa üzerindeki
kesmede tam seçilen envanter korunmuyordu. Küçük bir komut değişimi
`1e-12` toleransıyla tam sabit dışa aktarım olarak etiketlenebiliyordu.
Ara aralık kayıpları sonlu integral ve merkez zamanını da bozabiliyordu.

## Hesap ve karar sözleşmesi

Komut parçalarının zamanları, kesilmiş uç değerleri, maruziyetleri ve birinci
zaman momentleri verilen binary girdilerden rasyonel hesaplanır. Genel
toplam skaler çıktının temsil edilememesi, erken bir sonlu sorguyu veya
yanmayı engellemez. Doğrudan toplam/moment ihracı gerçek aralık ihlalinde
alan bilgili `DomainError` verir. Segment ve aralık kayıtları bozuk, boolean
ve sonlu olmayan sayısal alanları reddeder.

Tersleme mutlak tolerans veya boyutlu eğim eşiği kullanmaz. Tam sınırlar ve
rasyonel karekökler tam kalır; diğer karekökler çağıranın bağlamından bağımsız
1100 basamaklı Decimal bağlamında çözülür. Bekleme aralığında aynı maruziyet
için en erken zaman döner. İstek tam toplam maruziyeti aşarsa sonuç `None`
olur. Pozitif bir zamanın float'a sıfırla ihracı açık hatadır.

Kesme önce zaman sınırındaki tam maruziyetle karşılaştırılır; envanter
sınırı için terslenen zaman, seçilmiş tam maruziyeti yeniden tanımlamaz.
Eşit zamanlarda hedefin önceliği korunur; daha önce ulaşılan envanter sınırı
bekleme aralığı boyunca ertelenmez. Son tank kütleleri ve tüketimler bu
seçilmiş defterden kurulur. Küçük tüketim, iki gösterilen kütle aynı float'a
yuvarlansa da korunur. Profil son noktası da aynı son envanteri taşır.

Δv terslemesi ve ileri logaritma L0'ın ayrı incelenmiş yüksek hassasiyetli
çekirdeğini kullanır. Fiziksel çıktılar sonunda denetimli float'a çevrilir.
Yuvarlanmış zaman ile profil yeniden integre edildiğinde ortaya çıkan itki
farkı da kapanış kaydına katılır. Artıklar dışa doğru yuvarlanan sonlu üst
sınırlar olarak saklanır; sıfırdan farklı pozitif artık en küçük pozitif
float'tan daha küçükse bu float'la sınırlanır. Bu kayıt modelin fiziksel hata
veya belirsizlik zarfı değildir. Kütle defteri artığı ayrı araç/akış toplamları
arasındaki farktır; yuvarlanmış bütün çıktı alanlarının tam eşit olacağı
iddiası değildir.

Tam sabit dışa aktarım yalnız kullanılan komut değerleri gerçekten eşitse
işaretlenir. Bu işaretin anlamı profil biçimidir; bütün zamanların binary64
ile tam temsil edildiği veya aktarım yapılan başka bir çözücünün bütün
uç ölçeklerde doğrulandığı anlamına gelmez. Ayrı iki zamanın tek float'a
çökmesi halka açık aralık listesini dışa aktarırken açık hatadır; doğrudan
integral bu durumda hâlâ temsil edilebilir olabilir. Profil noktalarının
zamanları ayrı yuvarlanır ve anlık komut geçişlerinde yinelenebilir.

Kanıt yöntemi `analytic_profile_binary_ledgers_decimal1100_v2` olarak
değişir. Sonuç hash'leri ve kanonikleşen komut kayıtlarının hash'leri
değişebilir; şemalar aynı kalır. Hash ölçüm kaynağı doğrulaması değildir.
1100 basamak sabit bir sayısal seçimdir; evrensel doğru yuvarlama veya bütün
çıktılara ait ispatlanmış hata zarfı sunulmaz. Mevcut tolerans metaverisi de
böyle bir garanti olarak yorumlanmamalıdır.

## Bağımsız kontroller

25 yöntem pozitif çok küçük hedefleri, artan/azalan eğimleri, uzun aralıkları,
komşu float sınırlarını, beklemeleri, erken kesmeleri, moment/merkez zamanı
aralıklarını, fiziksel çıktı hatalarını, bozuk sorgu/kayıtları, gerçek sabitlik
kontrollerini, envanteri ve JSON dönüşümünü kapsar. Alt kontroller ayrı
yöntem sayılmaz.

200/300 basamaklı bağımsız referans kökleri global ikiye bölmeyle bulur;
üretimin kuadratik tersleme formülünü kopyalamaz. Maruziyet ve birinci moment
doğrudan yerel polinomların ayrı Decimal integrallerinden kurulur. Küçük Δv
hedefi L0 incelemesindeki integral serisiyle ayrıca denetlenir. İki referans
hassasiyeti aynı float'a yerleşmeden dört ulp karşılaştırması yapılmaz.
Seçilmiş sınır/defter örneklerinde tam bilinen binary eşitlikler istenir.

Dokuz ek seçilmiş kusur pozitif maruziyetin silinmesini, erişilemez isteğin
kabulünü, eğimin düzlenmesini, yanlış diskriminantı, eksik global momenti,
yuvarlanmış tank farkını, yaklaşık komutun tam sayılmasını, kesme toleransını
ve gösterilen zamanın hatasının gizlenmesini hedefler. Yalnız assertion
başarısızlığı yakalama sayılır; ilişkisiz istisna geçerli deney değildir.
Arşivlenmiş eski kaynakla sonuç `profile_historical_replay.json` içindedir.
Yeni kayıt sözleşmesi eksiklikleri aritmetik kusurla aynı şey sayılmaz.

Güncel tam kaynak bağlı sonuçlar `output/assurance/verification_summary.json`
içindedir. [Olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz donanım
kapısı açıktır; bu inceleme donanım TRL değerlendirmesi değildir.

## References

- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
