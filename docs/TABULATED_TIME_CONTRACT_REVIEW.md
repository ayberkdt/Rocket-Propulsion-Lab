# Tablolu yanma: zaman koordinatları ve olay sahipliği

Bu inceleme göreli zamanın örneklenmesini, kaydırılmış köprü koordinatını,
timeline olay tarihlerini ve Lunaris adaptörünün profil zaman sınırlarını kapsar.
Donanım saatini, deneysel zaman belirsizliğini veya bütün adaptör hesaplarını
geçerlemez. Kaynak dosyalar ve deney sonuçları `output/assurance` altında bağlanır.

## Yeniden üretilen kusurlar

Önceki varsayılan `max(1e-12, duration*1e-12)` penceresi, `1e-20 s` süren
bir profilin tamamını sınır sayabiliyordu. `7.5e-21 s` iç noktası sıfıra
taşınıyor, sol örnek önceki segmentten geliyor ve köprü kapalı sayılabiliyordu.
Normal ölçekli profilde de sınıra komşu temsil edilebilir noktalar yanlış
tek taraflı limit veya erken kapanış üretebiliyordu.

Lunaris adaptörü ayrıca farklı olay tarihlerini toleransla birleştiriyordu.
Göreli kısa profillerde bütün tarihler tek işarete düşebiliyordu. Büyük epoch
için 16-epsilon penceresi, mutlak saatte gerçekten ayırt edilebilen bir ULP
uzaklıktaki olayları da reddediyordu. Bunlar aynı arıza değildir: saat çözünürlüğü
kontrolü bir tolerans tahmininden türetilmemelidir.

## Göreli koordinat sözleşmesi

- Başlangıç tam sıfır, komşu segmentlerin sınırları tam eşit olmalıdır.
  Boşluk/örtüşme veya başlangıç kayması artık `1e-12` ile saklanmaz. Gerçek
  boşluk, açık bir sıfır itki/debi segmentiyle belirtilir. Asgari süre yoktur.
- Sahiplik `[start,end)` kuralıdır; `side='left'` yalnızca tam iç sınırda
  önceki segmenti seçer. Artifact son noktada son segmenti örnekler; köprüde
  ateşleme öncesi sol limit ve son noktanın sağ limiti kapalıdır.
- `boundary_tolerance_s` varsayılanı sıfırdır. `snap_to_boundary` isteğe bağlı
  pencerede en yakın sınırı seçer. Pencere, en kısa **tam** sürenin yarısından
  küçük olmalıdır. Koruma ve uzaklık karşılaştırmaları yuvarlanmadan yapılır;
  subnormal bir sürenin yarısının sıfıra yuvarlanması geçerli sıfır pencereyi
  reddetmez.
- `t-bias` tam rasyonel konumuyla sahiplik ve itki/debi örneklemesine girer.
  Son görüntülenen `effective_profile_time_s` binary64'tür. Bu alan bir sınıra
  yuvarlansa bile iç noktanın fazı veya türevdeki düzgünlük bilgisi değiştirilmez.

Bu, önceki yaklaşık süreklilik kabulünü bilinçli olarak değiştirir. Üretici
aynı olayın iki ucunda aynı zaman değerini kullanmalıdır; dosya okuyucusu
belirtilmemiş bir fiziksel süreyi sessizce ekleyemez veya silemez.

## Mutlak saat sözleşmesi

Her profil sınırı `b` için kanonik tarih, verilen binary değerlerin
`ignition+bias+b` tam toplamından bir kez binary64'e çevrilir. Ardışık iki sınır
aynı mutlak tarihe yuvarlanıyorsa bu saat profili ayırt edemez; epoch API'si ve
profil segmentasyon üretimi açık hata verir. Bir ULP uzaklıktaki farklı
kanonik tarihler kabul edilir. Göreli API aynı profili sıfır epoch'ta kullanabilir.

İki eski toplama sırasının sonlu sonuçları da yalnızca tek bir sınıra işaret
ediyorsa olay takma adı olarak kabul edilir. Kanonik tarih her zaman önceliklidir;
birden çok sınıra ait takma ad atılır. Kötü bir ara toplama sırasının küçük
profili kaybetmesi, iptalin tam korunmasıyla ayırt edilebilen kanonik saati
geçersiz yapmaz. Olay olmayan diğer zamanlarda `time-ignition-bias` doğrudan tam
koordinat olarak kullanılır. Bu kural keyfi bir çevre penceresindeki bütün
noktaları olaymış gibi sınıflandırmaz.

`PropagatorPropulsionBridge.evaluate_at_epoch` ortak giriş noktasıdır.
Olayın saltation farkı `time_event_jump_at_epoch` ile aynı saatten alınır.
Timeline ve Lunaris profil işaretleri aynı kanonik tarihleri kullanır.
Lunaris yalnızca eşit tarihleri tekilleştirir; RHS ilk çağrıda sağ, başka bir
kanonik segment sonu çağrısında sol limiti kullanır. Bu, sayısal tüketicinin
segment sonunda ve yeniden başlarken uyguladığı çağrı sözleşmesidir; her
entegratörün bu sözleşmeye uyduğu iddia edilmez.

`epoch_time_tolerance` eski, ihtiyatlı bir ölçek tahmini olarak erişilebilir
kalır. Artık timeline olay kimliğini, çözünürlük kararını veya Lunaris olay
birleştirmesini belirlemez. Bunun bir ölçüm belirsizliği bütçesi olduğu söylenemez.

## Kimlik, kanıt ve sınırlar

Artifact, köprü ve timeline için zorunlu hash alanı:
`time_numerics=exact_relative_clock_and_once_rounded_epoch_events_v1`.
Eksik/değiştirilmiş politika reddedilir. Eski belgeler kendi kaynak örneklerinden
yeniden üretilmelidir; eski hash üzerine yeni etiket yazmak yeterli değildir.
Bu kontrol ham JSON/CSV yapısının tamamının denetlendiği anlamına gelmez.

25 yeni zaman yöntemi; sıfır/subnormal/küçük/büyük ölçekler, komşu binary noktalar,
tam kaydırılmış konum, gerçek saat çakışması, bir ULP ile ayrılan olaylar,
1200/1600 basamaklı bağımsız Decimal tarihleri, ham politika ve normal ölçek
kontrollerini kapsar. Ek bir varyasyon yöntemi, serbest uçuşta tam yanmanın
zamanda ötelenmesinin son durumu değiştirmemesini ve zaman kaydırma türevinin
sıfır olmasını, bağımsız trapez kütle defteriyle birlikte sınar. Toplam 26 yeni
yöntem ve 16 seçilmiş semantik kusur eklenmiştir. Kusur deneylerinde
saf saat tablosu önbelleği de temizlenir; hata veren deney yakalama sayılmaz.

Eski segment tarama testinin beklenen sahipliği artık toleranslı taramadan
değil, aralık üyeliği ve tam eşitlikten türetilir; komşu `nextafter` noktaları
eklenir. Eski ihtiyatlı çözünürlük testinde ise gerçekten aynı binary tarihe
çöken iki olay kullanılır. Ayrı yeni kontrol, eski kuralın reddettiği ayırt
edilebilir bir ULP durumunu korur. Bu değişiklikler tolerans gevşeterek geçiş
sağlamak değildir; değişen zaman sözleşmesini matematiksel örneklerle sınarlar.

İlk tam koşuda 84 modülden ikisi başarısız oldu ve 504 kaynak/kanıt dosyası
değiştirilmeden arşivlendi. Rezerv oracle'ındaki eski zaman toleransıyla gerçek
bir kesimi silme kuralı kaldırıldı; artık gerçek ilk geçişin öncesi/sonrası
sınanıyor. Varyasyon testindeki RK tüketicisi de kendi bildiği hedef zamanı
yeniden toplamak yerine doğrudan kullanıyor, olaylarda tam eşitlik uyguluyor
ve hem durum hem saltation hesabında epoch girişini kullanıyor. Mevcut son
hassasiyet/sonlu fark assertion'ları ve zaman adımı değiştirilmedi. Bu tüketici
uyarlaması göreli API'ye gizli pencereyi geri getirmedi.

Lunaris tüketim integrali, rezerv kökü, birden fazla olayın ortak durum
geçişleri ve tüm fiziksel geçerleme bu zaman incelemesiyle kapanmaz. Sayısal
testlerin geçmesi TRL 5 veya sertifikasyon sağlamaz. Güncel tam koşu sonucu
`verification_summary.json`, bağımsız eski kaynak tekrarı
`time_contract_historical_replay.json` içinde tutulur.

## References

- [Python Fraction](https://docs.python.org/3.12/library/fractions.html): float girdinin tam rasyonel değerini koruma.
- [Python math.nextafter](https://docs.python.org/3.12/library/math.html#math.nextafter): komşu temsil edilebilir değerler.
- [Python Decimal](https://docs.python.org/3.12/library/decimal.html): bağımsız hassasiyet kontrollü tarih toplamları.
- [Donanım kanıt kapıları](PROPULSION_ASSURANCE.md): bu yazılım incelemesinin fiziksel sınırları.
