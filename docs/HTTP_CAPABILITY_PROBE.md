# Sidera yetenek sorgusunun yanıt sınırı

## Yeniden üretilen kusur

`/api/v1/integrations/sidera/capabilities` ilk Sidera yüklemesini HTTP
isteği içinde bekliyordu. Gerçek yerel Sidera ile iki saniyelik soket sınırı
hem paralel koşuda hem diğer modüller başlamadan çalıştırıldığında aşıldı.
Bu nedenle kaynak izolasyonu kök sorunu çözmedi. Önceki başarısız kayıtlar
`output/assurance/http_parallel_failure.*`, `http_isolated_failure.*` ve
`verification_storage_stage.json` olarak korunur. Tek bir ayrı koşunun
geçmesi, kusurun giderildiği kanıtı sayılmaz.

## API davranışı, sürüm 1.1

HTTP sorgusu artık isteğe bağlı paket incelemesini bir arka plan işçisinde
başlatır ve her çağrıda mevcut durumu hemen döndürür. Python entegrasyonundaki
senkron `inspect_sidera_capabilities` fonksiyonunun davranışı korunur.

| Durum | `probe_complete` | Yetenek alanları |
|---|---|---|
| `checking` | false | Henüz bilinmediği için gönderilmez |
| `supported`, `incompatible`, `unavailable` | true | İncelemenin son kaydı gönderilir |

İlk yanıt HTTP 200 ile sorgunun güncel durumunu taşır. Sonuç bilinmeden
`constant_burn_exact=false` gibi gözlenmiş bir destek kararı verilmez.
Arayüz `checking` sırasında bilgi verir ve tekrar sorgular. Tamamlanan
kayıt süreç ömrü boyunca saklanır; kurulum/yol değişiklikleri için sunucu
yeniden başlatılmalıdır. Aynı anda gelen istekler yeni incelemeler başlatmaz.

İşçi başlatılamazsa veya inceleme istisna/`SystemExit` üretirse destek
iddiası taşımayan `incompatible` kaydı yayımlanır. Takılmış bir dış import
`checking` kalır; yeni işçi doğurulmaz. İşçi daemon'dır. Bu, yabancı native
kod için süreç izolasyonu veya zorla iptal mekanizması değildir; böyle bir
kodun sunucunun geri kalanını etkileyemeyeceği iddia edilmez.

## Doğrulama

Gerçek HTTP testi, incelemeyi bir Event üzerinde kilitler. Kilit ancak
HTTP yanıtı alındıktan sonra açılır: senkron bekleyen uygulama iki saniyelik
soket sınırını aşar. Yanıtın `checking` olması, yetenek bayrağı uydurmaması
ve kilit açılınca doğru son kaydın dönmesi denetlenir.

Ayrıca 24 eşzamanlı yoklamanın yalnızca bir inceleme başlatması, inceleme
hatasının önbelleğe alınması, işçi başlatma hatası ve dış importun
`SystemExit` ile bitmesi sınanır. Gerçek isteğe bağlı Sidera kurulumu için
son durum gözlenene kadar her HTTP isteğinin iki saniyelik sınırı korunur;
asenkron işin tamamlanmasını gözlemek için testte ayrı 20 saniyelik pencere
vardır. Bu bir kullanıcı görevinden türetilmiş performans gereksinimi değildir.

Son regresyon HTTP dahil bütün modülleri dört işçili koşuya geri alır.
Kontroller kaldırılmamıştır; önceki başarısız sonuçlar son başarılı koşuyla
örtülmez. Yük altındaki genel gecikme dağılımı ve tarayıcıda görsel doğrulama
bu çalışmada değerlendirilmemiştir.

Bu düzeltme yazılım yanıt davranışını kapsar. Motor/besleme geçerlemesi,
ileri Sidera propulsion-provider tüketimi veya donanım TRL kanıtı değildir.

## References

- [Python 3.12 threading: Thread, Lock, Event and daemon shutdown semantics](https://docs.python.org/3.12/library/threading.html)
