# Test tanılaması ve güncel kanıt raporu

## Gözlenen native çöküş ve kontrollü inceleme

Geniş koşuda saldırgan JSON okuyucu testi sırasında Python 3.12.1 süreci,
20 saniyelik yığın izi dökümü başladıktan sonra erişim ihlaliyle kapandı.
Bu koşu başarılı sayılmadı. Orijinal kayıtlar
`output/assurance/json_reader_native_crash.txt` ve
`regression_native_crash_summary.json` olarak korundu. Aynı JSON modülünün
ayrı koşusu dört test ve 324 alt kontrolle geçti; tek başına bu tekrar,
çöküşün giderildiği kanıtı sayılmaz.

İtki kütüphanesini veya test verilerini içermeyen, yalnızca Python yığın
çerçevelerini değiştirip iz toplayan küçük alt süreçler ayrıca çalıştırıldı:

| Denetim | İki kısa koşunun sonucu |
|---|---|
| C faulthandler zamanlayıcısı | Bir erişim ihlali, bir normal çıkış |
| Python iş parçacığından GIL ile iz toplama | İki normal çıkış |
| İz zamanlayıcısı yok | İki normal çıkış |

Tekrarlanabilir örnek ve byte özetleri `traceback_race_probe.py` ve
`traceback_probe_results.json` içinde tutulur. Kontrol çöküşünün çıkış kodu
`3221225477` idi. Orijinal geniş koşunun eski aracı sayısal çıkış kodunu
kaydetmediği için bu kod orijinal koşuya atanmaz.

Upstream CPython kayıtlarında benzer çerçeve/yığın izi yarışları vardır.
Bu bulgular tanılama yarışını özgün çöküş için güçlü aday yapar; özgün
native C yığını alınmadığından birebir hata kimliği kesinleştirilmiş değildir.
Denetimler uygulamanın bütün native hata olasılıklarını dışlamaz.

## Son test protokolü

C zamanlayıcısıyla canlı çerçevelerin eşzamansız gezilmesi yerine,
`scripts/pytest_diagnostics.py` Python Timer callback'i kullanır. İz
toplamaya girerken GIL tutulur. Aynı 20 saniyelik tanılama gecikmesi
korunur; gerçek fatal hata işleyicisi kapatılmaz. Test seçimi, girdiler ve
assertion'lar değiştirilmez. Ana süreç her modülü eskisi gibi en çok
120 saniye çalıştırır; bu sert sınır tanılama thread'inden bağımsızdır.

GIL'i bırakmayan takılmış native bir çağrı Python callback'ini geciktirebilir.
Bu durumda 20 saniyelik iz en iyi çabadır; 120 saniyelik üst süreç sınırı
işlemeye devam eder. Bu protokol free-threaded Python için geçerli değildir.
Çıkış kodları artık modül kayıtlarında ayrıca tutulur.

Gerçek bir alt pytest sürecinde bekleyen testin yığın izinin yazılması ve
aynı testin sonucu korunarak tamamlanması denetlenir. Bu bir model veya
donanım testi değildir; test altyapısının tanılama kontrolüdür.

## Eksik kanıtın eski yeşil raporu bırakması

İncelemede ayrı bir raporlama kusuru yeniden üretildi: bir süreç JUnit
üretmeden çökerse, özetleyici istisna veriyor fakat önceki başarılı
`verification_summary.json` yerinde kalıyordu. Ayrıca tamamlanmış ama
başarısız bir özetin komut çıkış kodu başarısızlığı yansıtmıyordu.

Özetleyici artık okumaya başlamadan eski başarıyı geçersiz kılar. Eksik,
bozuk veya güncel olmayan kanıt `incomplete_or_invalid` ve
`complete_suite_passed=false` üretir; sayılamayan testler sıfır başarılı
gibi gösterilmez. Başarısız/eksik değerlendirme komutu 1 ile çıkar.
Kaynak ve test kimliği doğrulamasının katılığı korunur. Eksik bir JUnit
dosyası ve önceki yeşil raporla bu davranış gerçek dosyalar üzerinde sınanır.

Güncel rapor yalnızca yazılım doğrulama kanıtını kapsar. Bu düzeltmeler
deney gerçekliğini, fiziksel kalibrasyonu veya donanım TRL seviyesini
belirlemez; ilgili kapılar açık kalır.

## Seçilmiş mutasyonların genel sonuca katılması

Yeni bir raporlama kusuru da geçerli bir sentetik JUnit ve kaynak kaydıyla
yeniden üretildi: mutasyon `survived` veya `invalid` olsa bile eski özetleyici
genel durumu `complete` yapıyor ve 0 ile çıkıyordu. Boş/eksik, yinelenmiş veya
yapılandırılmamış mutasyon satırlarını da reddetmiyordu.
`output/assurance/mutation_summary_before_fix.xml` bu karşı örneği tutar.

`complete_suite_passed` yalnız toplanan pytest yöntemlerinin durumudur.
`curated_mutation_runs_passed`, kaynakları güncel iki seçilmiş mutasyon
koşusunun her birinin bütün tanımlı kusurları assertion ile yakalamasını ister.
`software_checks_passed` ikisinin birlikte sonucudur; genel durum ve komut
çıkışı bu son alanı kullanır. Geçerli bir JUnit başarısı, yakalanamayan veya
geçersiz bir kusuru gizlemez. Geçerli fakat başarısız mutasyon kaydında pytest
sonucu korunur, genel yazılım sonucu başarısız kalır.

Özetleyici, araçların kaynak hash'ine bağlı açık `MUTATIONS` tuple'larından
kusur adlarını AST ile okur; araç kodunu bu amaçla çalıştırmaz. Her tanımlı
kusur için tek sonuç gerekir. Durumun assertion/hata sayılarıyla tutarlı
olması da denetlenir. Testte aynı geçen JUnit korunurken yakalanamayan,
geçersiz ve yakalanmış mutasyonlar ile eksik/yinelenen/yanlış kayıtlar
değiştirilir; geçerli kontrolün başarılı kalması da zorunludur.
Bu tutarlılık denetimi deney gerçekliği veya raporun bağımsız onayı değildir.

Sonraki türev düzeltmesinden önce 497 yöntemlik başarılı koşunun 285 depo
girdisi, JUnit/log kayıtları ve seçilmiş mutasyon raporları
`output/assurance/history/pre_mixture_derivative_fix` altında ayrı kopyalandı.
Kopya girdilerin hash'leri önceki kaynak anlık görüntüsüyle karşılaştırıldı.
Dış bağımlılıklar bu arşivde yalnız hash ve sürüm bilgisiyle kayıtlıdır;
tam çalıştırılabilir ortam veya donanım yeterlilik arşivi değildir.

## Test toplama aşamasının başarısızlığı

Türev incelemesindeki ilk geniş koşu test toplamayı 30 saniyede bitiremedi.
Eski araç bu yolda traceback ile çıkıyor, kısmi koleksiyon çıktısını
saklamıyor ve önceki `regression/summary.json` dosyasını bırakıyordu.
Zaman aşımı ve sıfır olmayan toplama çıkışı için bu eski-yeşil davranış
geçici gerçek dosyalarla yeniden üretildi (`collection_reporting_before_fix.xml`).

Araç artık çalışmanın başında eski başarıyı geçersiz kılar. Toplama
başarısızlığı/zaman aşımı, eldeki stdout/stderr, çıkış kodu ve aşama süresiyle
başarısız bir güncel kayıt üretir. Test modülleri başlatılmaz; beklenen modül
listesi ayrı tutulur. Özetleyici tamamlanmamış toplamayı kabul etmez. Araç
çıkışı, atlanan test veya değişmiş kaynak dahil genel paket sonucunu yansıtır.

Tanılama alt sürecinde aynı test kümesi bir kez 18.30 saniyede tamamlandı;
pytest kendi toplama bölümünü 11.06 saniye bildirdi. Bu gözlem ilk zaman
aşımının nedenini kanıtlamaz; o koşunun kısmi çıktısı eski araç tarafından
korunmamıştır. Tanılama denemesinin 90 saniyelik ayrı üst sınırı, gerçek
regresyon protokolünün sınırı değildir. Regresyonda toplama **30 saniye**,
modül **120 saniye** ve **dört işçi** olarak kalır. Başlangıç gecikmesinin yük
altındaki dağılımı değerlendirilmemiştir. Yeni bir başarılı koşu önceki
başarısız denemeyi geçmişten silmez.

## Canlı dış bağımlılıktan sürüme sabitlenmiş kaynak kopyasına

Köprü debi düzeltmesinin 69 modüllük koşusundan sonra Sidera'nın canlı
`ui/app.py` kaynağı değişti. Test yöntemleri geçmiş olsa da son özet dış
bağımlılık hash farkını yakalayıp güncel başarıyı reddetti. Fark
`bridge_flow_external_dependency_changes.json` içinde saklanır; bu deneme
başarılı son doğrulama olarak kullanılmaz.

Özgün dış dizinler değiştirilmeden, o anda bulunan Sidera/Lunaris Python
paketleri ve kaynak içi kaynak dosyaları çalışma alanına kopyalandı. Kopya
öncesi/sonrası özgün envanter ve kopya içerik hash'leri karşılaştırıldı.
Node yapı bağımlılıkları, bytecode ve günlükler kapsam dışıdır. Proje kökü
dış veri dosyaları ve kurulu ikili ortam klonlanmamıştır. Kopya manifesti ve
kaynak kökü `bridge_flow_pinned_dependency_snapshot.json` ile tanımlanır;
paket dışı ikili/sayısal bağımlılıkların sürümleri ayrıca kaydedilir.

Yeni Python süreçleri hem toplama/ana süreçte hem alt test süreçlerinde bu
kopya kökü kullanır. Sidera ve Lunaris import kökleri kontrol edildi; mevcut
yerel sağlayıcı testlerinin 16 yöntemi bu kopyalarda geçti. Ardından ayrı
API/yetenek incelemesinde yüklenen 79 Python modülünün kopyadan geldiği ve
kopyalanan dosyaların değişmediği kontrol edildi. Bu denetim, mümkün olan
her dinamik import veya ikili bağımlılık için genel bir kanıt değildir.

Son geniş test raporu canlı dış dizinlerin daha sonraki değişiklikleri için
değil, manifestteki kaynak kopyası için geçerlidir. Toplama/modül süreleri,
seçilen testler ve dört işçi protokolü korunur. Bağımlılık kopyası ve bu
yazılım sonuçları donanım olgunluğu veya fiziksel geçerleme kanıtı değildir.

## References

- [Python 3.12 faulthandler](https://docs.python.org/3.12/library/faulthandler.html)
- [CPython 3.12.1 faulthandler implementation](https://raw.githubusercontent.com/python/cpython/v3.12.1/Modules/faulthandler.c)
- [CPython issue 116008](https://github.com/python/cpython/issues/116008)
- [CPython issue 131580](https://github.com/python/cpython/issues/131580)
- [Donanım ve geçerleme açıkları](PROPULSION_ASSURANCE.md)
