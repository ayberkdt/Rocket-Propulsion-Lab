# Doğrudan tank envanteri yardımcıları

İnceleme `burns/inventory.py` içindeki dört çağrı ve bunların beş iç
yardımcısını kapsar. L0/L1 çözücüleri bu modülün kapsam denetimini çağırır;
süre ve son kütle hesaplarını kendi ayrı defterlerinde yapar. Bu nedenle
doğrudan yardımcıdaki bir kusur, ayrıca çağrı yolu gösterilmeden güncel
çözücü sonucuna atfedilmez.

Önceki uygulamada eksik tank `(inf, None)` sonucuna, NaN zaman NaN envantere,
negatif zaman artan envantere dönüşüyordu. Aynı tankın küçük akış yolları
ardışık float toplamında kaybolabiliyor; mutlak/göreli rezerv toleransı
gerçek bir aşımı gizleyebiliyordu. Örneğin yük `1e16`, rezerv sıfır ve iki
yolun debileri `1e16` ile `1` iken bir saniye güvenli değildir.

## Hesap ve sınır sözleşmesi

Her izlenen yolun kimliği, akışkan anahtarı ve rolü tek bir geçerli tankla
eşleşir. Sıfır debili izlenen yollar da eşleşmek zorundadır. Yinelenen tank
kimlikleri reddedilir. Koleksiyon/kayıt türleri, alanlar, boolean izleme
bayrağı ve sonlu negatif olmayan kütle/debi/zaman değerleri kontrol edilir.
Bilinen rol metinleri enum karşılıklarıyla karşılaştırılır; çağıranın
kayıtları değiştirilmez. İzlenmeyen akış tank envanterini tüketmez.

Debiler, yük, rezerv ve zaman, kabul edilen binary64 değerlerinin tam
rasyonel karşılıklarıyla hesaplanır. Aynı tanka ait yollar son çıktıdan önce
toplanır. `flow_by_tank` tek kez yuvarlanmış debi döndürür. Süre/kütle
yardımcıları bu yuvarlanmış çıktıyı geri beslemez. Kapsam denetimi bir debi
çıktısının float aralığına sığmasına bağlı değildir.

`tank_masses_at` rezerv kararını yuvarlamadan önce verir. Gerçek aşım,
toleransla kabul edilmez. Son kütle yalnız çıktı sınırında yuvarlanır;
pozitif kütlenin sıfıra yuvarlanması veya taşma alanı belirtilmiş
`DomainError` üretir. Sıfır kalan kütle geçerlidir. Bu sözleşme, girdilerde
ölçüm belirsizliği bulunmadığı anlamına gelmez.

`tank_duration_limit`, fiziksel kesme anına en yakın float yerine rezervi
aşmayan **en büyük binary64 zamanını** döndürür. Kesin kesme anı `1/10` ise
`.1` bir miktar büyüktür; dönen sınır `.1` değerinin önceki komşusudur.
Bu yüzden temsil edilebilir güvenli zaman kümesi keyfî olarak daraltılmaz.
Bağlayıcı tank yuvarlamadan önce seçilir; gerçek eşitlikte giriş sırası korunur.

- Pozitif tüketim yoksa sonuç `(inf, None)` olur.
- Pozitif kesme anı en küçük pozitif float değerinden küçükse sınır sıfırdır;
  bu, tankın başlangıçta rezervde olduğu iddiası değildir.
- Kesme anı en büyük sonlu float değerinden büyükse `(inf, tank_id)` döner:
  bütün sonlu zaman girdileri rezerv açısından güvenlidir. Fiziksel kesme
  anının sonsuz olduğu iddia edilmez; sonsuz zaman sorgusu yine reddedilir.
- Sonlu güvenli zaman, bütün son kütlelerin çıktı aralığına sığacağını garanti
  etmez. Pozitif çok küçük kütlenin temsil edilememesi ayrı çıktı hatasıdır.

Bu ayrım doğrudan yardımcıların davranış değişikliğidir. L0/L1 çözümündeki
tam olay defterleri ve yuvarlanmış gösterim zamanı politikası değiştirilmez.
Bu yardımcıların yeni bir seri hale getirme biçimi veya doğrulama sertifikası
ürettikleri iddia edilmez.

## Karşı kanıt ve test kapsamı

`test_burn_legacy_inventory.py` içindeki 16 yöntem gerçek rezerv aşımı,
eksik/uyumsuz/yinelenen eşleme, kötü zaman ve kayıt türleri, akış permütasyonu,
yuvarlama öncesi bağlayıcı tank seçimi, sıfır/izlenmeyen yollar ve gerçek çıktı
aralığı hatalarını sınar. Uç ölçekler yalnız sayısal örneklerdir.
48 sabit tohumlu örnekte kalan kütle ve komşu zaman kararı bağımsız Decimal
defteriyle 800/1200 basamakta karşılaştırılır; uç örnekler ayrıca daha yüksek
hassasiyetle kontrol edilir. Referans, üretim kodunun Fraction yardımcılarını
çağırmaz. Bazı testler motorun genel performans alanlarından bağımsız olarak
akış defterini sınamak için kayıt kopyalarındaki yalnız ilgili alanları değiştirir;
bu örnekler fiziksel olarak tutarlı motor çalışma noktaları sayılmaz.

Sekiz seçilmiş kusur kapsam denetimini atlama, negatif zaman kabulü,
rezerv toleransı, erken toplam yuvarlaması, yukarı yuvarlanmış süre,
yuvarlamadan sonra bağlayıcı seçim, yinelenen tankı ezme ve küçük kütleyi
sıfırlamayı hedefler. Yalnız assertion başarısızlığı yakalama sayılır;
ilişkisiz çalışma zamanı hatası geçerli yakalama değildir.

Güncel testler arşivlenmiş eski kaynak üzerinde yeniden çalıştırılır;
yüklenen paket kökü, kaynak ve test hash'leri kaydedilir. İlk keşif koşusunda
uç örneğin beklenen değeri yanlış yazılmıştı; düzeltme bağımsız Decimal
referansına göre yapıldı. Son test dosyasıyla tarihsel tekrar esas karşı kanıttır.
Tam koşu ve mutasyon sonuçları `output/assurance/verification_summary.json`
içindedir. Bu sınırlı inceleme bütün yazılımın doğruluğunu veya donanımın TRL 5
olgunluğunu kanıtlamaz; dokuz donanım kapısı açık kalır.

## References

- [Python 3.12 Fraction: float değerinin tam rasyonel karşılığı](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 nextafter: komşu float değeri](https://docs.python.org/3.12/library/math.html#math.nextafter)
- [Donanım olgunluğu kapsamı ve açık kapılar](PROPULSION_ASSURANCE.md)
