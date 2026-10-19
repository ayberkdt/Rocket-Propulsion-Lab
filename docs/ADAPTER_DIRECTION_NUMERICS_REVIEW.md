# Adaptörlerde yön normalizasyonu ve son vektör izdüşümü

İnceleme Sidera sabit yanma dosyasını, ensemble üretici/okuyucusunu ve
Lunaris dürtü eşdeğerindeki vektör hesabını kapsar. Native kontroller
sabit yanmanın kamuya açık ivme/debi fonksiyonu ile dürtünün kamuya açık
uygulanma sınırıdır. Weighted/ileri provider yayılımının, bütün yörünge
çözücüsünün veya fiziksel yönlendirme donanımının geçerlemesi değildir.

## Yeniden üretilen kusurlar

Sidera sabit ve ensemble yolları doğrudan kareler toplamıyla norm
hesaplıyordu. Sonlu `1e308` ölçekli yön sıfır birim vektöre dönüşebiliyor;
çok küçük, sıfır olmayan yön reddedilebiliyordu. Ensemble çıktısının sıradan
bir eğik yönle yeniden açılması da tekrar normalizasyon yüzünden hash
kapsamındaki koordinatları değiştirebiliyordu.

Lunaris dürtü yolunda ham yönün normu ve `Δv * yön / norm` işlem sırası
ayrıca ara taşma/alt taşma üretebiliyordu. Scalar Δv mevcutken vektör normu
onu karşılamıyor, native dürtü yanlış hız/kütle sonucu verebiliyordu.
Yalnız eksen yönleri veya yalnız norm karşılaştırması bileşenlerin işaret
ve sıra hatalarını yakalamaya yetmez.

## Çıktı ve kayıt sözleşmesi

Üç yol, daha önce incelenmiş geometri çekirdeğinin ölçeklenmiş `hypot`
norm tahminini ve verilen binary koordinatlar üzerindeki rasyonel son
izdüşümünü kullanır. Norm tahmini tam karekök değildir; evrensel doğru
bileşen yuvarlama garantisi sunulmaz. Orijinal koordinatlar son çıktı
oluşana kadar korunur.

Sidera dosyasında çıktı birim vektördür. Sıfır olmayan bir birim bileşenin
binary64 aralığından kaybolması alan adıyla `DomainError` verir; bütün küçük
ham yönler reddedilmez. Lunaris'te çıktı doğrudan Δv vektörüdür: önce birim
vektörü float'a çevirip sonra Δv ile çarpmak yerine son bileşen dönüştürülür.
Böylece tek başına saklanamayan birim bileşen, büyük Δv ile kurtulabiliyorsa
korunur. Son Δv bileşeni gerçekten temsil edilemiyorsa açık alan hatası verir.

Sidera sabit dosyası ve hashed ensemble belgesi
`scaled_hypot_exact_binary_projection_v1` yöntem kaydını taşır. Ensemble
okuyucusu bu alanın değiştirilmesini/silinmesini reddeder. Eski yöntemsiz
belgeler bu okuyucuda yeniden oluşturulmalıdır; sessiz yöntem yükseltmesi
yapılmaz. Hash ölçüm kaynağı veya gönderici kimliği doğrulaması değildir.

Ensemble yeniden açılırken kayıtlı koordinatlar `vector3` sınırından geçer
ve normlarının birden uzaklığı en fazla birin dört ulp'si olarak denetlenir.
Bu, float gösterimi için sayısal kontroldür; fiziksel açı hata payı değildir.
Geçerli koordinatlar yeniden ölçeklenmeden aynen korunur. Yeniden hash
hesaplanmış olsa bile birim olmayan veya boolean içeren yönler reddedilir.
Kaynak ensemble kimliği korunur; dış handoff hash'i yöntem ve yeni
koordinat hesabı nedeniyle değişebilir.

## Bağımsız ve native kontroller

14 saldırgan yöntem en büyük ve subnormal eksenleri, işaretli bileşenleri,
ham norm taşmasını, tam ikinin kuvvetleriyle ölçek değişimini, bütün eksen
sıralamalarını ve sabit tohumlu 32 anizotropik yönü kapsar. 200/300 basamaklı
Decimal norm/son vektör hesaplarının aynı float'a yerleşmesi doğrulanır;
sekiz ulp sınırı her bileşene ayrı uygulanır. Genel uçuş doğruluğu sınırı
bu testlerden çıkarılmaz.

Ek kontroller strict JSON, yinelenen ensemble açılışında aynı koordinat ve
hash, yöntem bütünlüğü, geçersiz/sıfır girdiler, gerçek bileşen aralık kaybı
ve büyük Δv ile kurtulan küçük bileşeni denetler. Alt kontroller ayrı test
veya ayrı fiziksel kusur olarak sayılmaz.

Sabitlenmiş Sidera kaynağından gerçek `FiniteBurn`/`ManeuverPlan` nesneleri
kurulur ve kamuya açık ivme/debi fonksiyonu inertial, RIC ve VNB çerçevelerinde
basit, bağımsız bilinen bazlarla karşılaştırılır. Lunaris'in gerçek dürtü
nesnesi kamuya açık fonksiyonla uygulanır; hız ve kütle değişimi karşılaştırılır.
Native dürtü uygulama kontrolü olağan sonlu Δv büyüklüğündedir; uç Δv
büyüklüğündeki native norm/roket bağıntısının aralık dayanıklılığı onaylanmaz.
Pinned kaynak kopyası, kurulu ikili ortamın veya donanımın kopyası değildir.

Beş ek seçilmiş kusur işaretlerin silinmesini, küçük bileşenlerin kesilmesini,
iki eksenin yer değiştirmesini, vektör büyüklüğünün atılmasını ve kayıtlı
birim yön denetiminin kaldırılmasını hedefler. Modül işlevi mutasyonları
adaptörlerdeki yüklü takma adlara da uygulanır; sınıf fixture'ları unittest
suite yaşam döngüsüyle kurulur. Yalnız assertion başarısızlığı yakalama
sayılır; ilişkisiz istisna geçerli deney değildir.

Arşivlenmiş eski paketle mevcut test dosyasının sonucu
`adapter_direction_historical_replay.json` içindedir. Yeni yöntem/kayıt
sözleşmesi eksiklikleri aritmetik kusurla eş tutulmaz. Güncel kaynak bağlı
tam sonuçlar `output/assurance/verification_summary.json` içindedir.
[Olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz donanım kapısı açıktır.

## References

- [Python 3.12 hypot](https://docs.python.org/3.12/library/math.html#math.hypot)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
