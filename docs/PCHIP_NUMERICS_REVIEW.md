# PCHIP eÄŸimleri ve kÃ¼bik harita aritmetiÄŸi

Ä°nceleme ortak PCHIP Ã§ekirdeÄŸini, `MonotoneCubicCurve` deÄŸer/tÃ¼rev hesabÄ±nÄ± ve
`RectilinearPerformanceSurface` iÃ§indeki `cubic_hermite` yolunu kapsar.
Gaz kelebeÄŸi eÄŸrisi ve kÃ¼bik performans haritasÄ± aynÄ± eÄŸim Ã§ekirdeÄŸini kullanÄ±r.
Darbe planlayÄ±cÄ±sÄ± da `MonotoneCubicCurve` Ã¼zerinden bu hesabÄ± tÃ¼ketir.
Ã‡ok doÄŸrusal harita yolu veya bÃ¼tÃ¼n performans modelleri iÃ§in tamamlanmÄ±ÅŸ
inceleme iddiasÄ± deÄŸildir. Veriler sentetiktir; motor/besleme donanÄ±mÄ± geÃ§erlenmez.

## Yeniden Ã¼retilen kusurlar

`1e-200` Ã¶lÃ§ekli pozitif ardÄ±ÅŸÄ±k sekantlarÄ±n Ã§arpÄ±mÄ± sÄ±fÄ±ra yuvarlanÄ±yordu.
Ä°ÅŸaret denetimi bunu sÄ±fÄ±r veya dÃ¶nÃ¼ÅŸ noktasÄ± sanÄ±yor, dÃ¼ÄŸÃ¼m eÄŸimlerini
sÄ±fÄ±rlÄ±yordu. Ã–rnekte doÄŸru Ã¶rnek `1.25e-200` iken `1.15625e-200` Ã§Ä±kÄ±yor;
`1e200` rated debi Ã§arpanÄ± sonrasÄ±nda `1.25` yerine `1.15625 kg/s` oluÅŸuyordu.
Bu Ã¶rnek fiziksel motor Ã§alÄ±ÅŸma noktasÄ± Ã¶nerisi deÄŸildir.

Ek saldÄ±rÄ±lar ÅŸu kusurlarÄ± da yeniden Ã¼retti:

- SonuÃ§ sonlu iken koordinat geniÅŸliÄŸinin veya ordinat farkÄ±nÄ±n taÅŸmasÄ±.
- KÃ¼Ã§Ã¼k aralÄ±kta bÃ¼yÃ¼k sabit deÄŸer iÃ§in sÄ±fÄ±r tÃ¼rev yerine sonlu olmayan sonuÃ§.
- Subnormal koordinat/ordinatlardan hesaplanabilen birim tÃ¼revin kaybÄ±.
- BÃ¼yÃ¼k karÅŸÄ±t terimlerin yanÄ±ndaki kÃ¼Ã§Ã¼k deÄŸerin erken koordinat/terim yuvarlamasÄ±yla bozulmasÄ±.
- Sonlu yÃ¼zey deÄŸerleri ve ilk gradyanlar iÃ§in gereken, float'a sÄ±ÄŸmayan karÄ±ÅŸÄ±k ara tÃ¼revlerin bozulmasÄ±.
- GerÃ§ek Ã§Ä±ktÄ± aralÄ±ÄŸÄ± hatalarÄ± ve bozuk doÄŸrudan Ã§ekirdek girdilerinin aÃ§Ä±k alan hatasÄ± vermemesi.

Bu ilk incelemedeki 19 yÃ¶ntem arÅŸivlenmiÅŸ eski paket kaynaklarÄ±na karÅŸÄ± yeniden Ã§alÄ±ÅŸtÄ±rÄ±ldÄ±:
18 yÃ¶ntem baÅŸarÄ±sÄ±z, anizotropik ayrÄ±labilir kontrol baÅŸarÄ±lÄ± oldu. YÃ¼klÃ¼ bÃ¼tÃ¼n
paket modÃ¼llerinin arÅŸiv kÃ¶keni ve kaynak baytlarÄ± kontrol edildi. BÃ¶ylece
referans hesabÄ±nÄ±n kendi hassasiyet sorunu, Ã¼retim kusuru diye sayÄ±lmadÄ±.
KayÄ±t: `output/assurance/pchip_historical_replay.json` ve `.xml`.

## Hesap ve Ã§Ä±ktÄ± sÃ¶zleÅŸmesi

Sekant iÅŸaretleri Ã§arpÄ±larak sÄ±nanmaz. GeniÅŸlikler, sekantlar, aÄŸÄ±rlÄ±klÄ± harmonik
ortalamalar ve uÃ§ nokta sÄ±nÄ±rlamalarÄ± verilen binary girdiler Ã¼zerinde rasyonel
tutulur. Ä°ki dÃ¼ÄŸÃ¼m doÄŸrusu, dÃ¶nÃ¼ÅŸ noktasÄ±, dÃ¼z bÃ¶lÃ¼m ve sekantÄ±n Ã¼Ã§ katÄ± olan uÃ§ sÄ±nÄ±rÄ±
aynÄ± PCHIP tanÄ±mÄ±nÄ± izler. `pchip_slopes` doÄŸrudan Ã§aÄŸrÄ±da dÃ¼ÄŸÃ¼m eÄŸimlerini
kontrollÃ¼ float olarak verir; gerÃ§ekten temsil edilemeyen eÄŸim aÃ§Ä±k hata verir.

Ä°ki tÃ¼ketici iÃ§teki rasyonel eÄŸimleri kullanÄ±r. Ara eÄŸimin float'a sÄ±ÄŸmamasÄ±
tek baÅŸÄ±na bÃ¼tÃ¼n eÄŸrinin reddi deÄŸildir: baÅŸka bir noktadaki dÃ¼z bÃ¶lÃ¼mÃ¼n
deÄŸeri/tÃ¼revi temsil edilebilir. Bir boyutlu polinom bir kez ortak paydalÄ±
tamsayÄ± kuvvet katsayÄ±larÄ±na derlenir. Her sorguda binary koordinat farkÄ± ve
polinom toplamÄ± tam tamsayÄ± aritmetiÄŸiyle hesaplanÄ±r. Son deÄŸer ve analitik tÃ¼rev ayrÄ± Ã§Ä±ktÄ±
dÃ¶nÃ¼ÅŸÃ¼mleriyle verilir; sÄ±fÄ±r olmayan temsil edilemez alan sessiz sÄ±fÄ±r olmaz.

KÃ¼bik haritada dÃ¼ÄŸÃ¼m deÄŸerleri, karÄ±ÅŸÄ±k tÃ¼rev tablolarÄ±, hÃ¼cre koordinatlarÄ±,
baz fonksiyonlarÄ± ve tÃ¼m tensor toplamlarÄ± aynÄ± dÃ¶nÃ¼ÅŸÃ¼m kuralÄ±nÄ± kullanÄ±r.
Ã–rneÄŸin `1e-180` eksen aralÄ±klarÄ±nda Ã¼Ã§Ã¼ncÃ¼ karÄ±ÅŸÄ±k tÃ¼rev float'a sÄ±ÄŸmayabilir;
son deÄŸer ve ilk gradyanlar yine sonlu olabilir. Bu ara alanlar erken yuvarlanmaz.
Her ihraÃ§ edilen kuvvet/debi katsayÄ±sÄ± ve gradyan kendi alan adÄ±yla denetlenir.

PCHIP'in bir boyutlu ÅŸekil korumasÄ±, N boyutlu iÃ§ bÃ¶lgede pozitiflik veya taÅŸmama
kanÄ±tÄ± deÄŸildir. Mevcut pozitif olmayan kÃ¼bik katsayÄ± ret davranÄ±ÅŸÄ± sÃ¼rer.
Rasyonel ara hesap, Ã¶lÃ§Ã¼m doÄŸruluÄŸu, fiziksel belirsizlik veya bÃ¼tÃ¼n kabul
edilen girdiler iÃ§in bu testlerden tÃ¼retilmiÅŸ evrensel hata sÄ±nÄ±rÄ± saÄŸlamaz.

Gaz kelebeÄŸi ve darbe manifestlerine `curve_numerics`, kÃ¼bik yÃ¼zey belgesine
`interpolation_numerics` eklendi; ilgili model/yÃ¼zey hash'leri deÄŸiÅŸir.
Okuyucular bilinen politikalarÄ± karÅŸÄ±laÅŸtÄ±rÄ±r. Ã–nceki kimlikler gÃ¼ncel hesap
anlamÄ±yla sessizce eÅŸ tutulmaz. Ã‡ok doÄŸrusal haritalarÄ±n ayrÄ± sayÄ±sal kimliÄŸi
ve hesap yolu [Ã§ok doÄŸrusal incelemede](MULTILINEAR_NUMERICS_REVIEW.md) aÃ§Ä±klanÄ±r.

## Oracle ve kusur enjeksiyonu

BaÄŸÄ±msÄ±z referans PCHIP tanÄ±mÄ±ndan dÃ¼ÄŸÃ¼m eÄŸimlerini kurar; deÄŸer/tÃ¼revi Ã¼retimdeki
kuvvet polinomu yerine Bezier kontrol noktalarÄ± ve de Casteljau yoluyla bulur.
Ä°lk 200/300 basamak denemesi bazÄ± uÃ§ eÄŸimlerin float yuvarlama orta noktasÄ±nda
aynÄ± sonuca yerleÅŸmedi. KayÄ±t `pchip_reference_precision_mismatch.xml` iÃ§inde
korunur. KarÅŸÄ±laÅŸtÄ±rma toleransÄ± deÄŸiÅŸtirilmeden referans 1200/1600 basamaÄŸa
Ã§Ä±karÄ±ldÄ±; iki referansÄ±n float sonuÃ§larÄ± eÅŸleÅŸmeden uygulama karÅŸÄ±laÅŸtÄ±rÄ±lmaz.
Bu Ã¶rneklenmiÅŸ referans tutarlÄ±lÄ±ÄŸÄ±dÄ±r; her olasÄ± veri iÃ§in ispat deÄŸildir.

19 yÃ¶ntem uÃ§ farklarÄ±, iÅŸaret/sÄ±nÄ±rlama kurallarÄ±nÄ±, 12 sabit tohumlu eÄŸriyi
dÃ¶rt `2^Â±500` Ã¶lÃ§ek birleÅŸiminde, sÄ±fÄ±r ve gerÃ§ek aralÄ±k hatalarÄ±nÄ±, gaz
kelebeÄŸi tÃ¼keticisini, bir/iki/Ã¼Ã§ boyutlu haritalarÄ± ve politika kaydÄ±nÄ± kapsar.
Anizotropik ayrÄ±labilir harita, baÄŸÄ±msÄ±z bir boyutlu Bezier yollarÄ±nÄ±n deÄŸer
ve tÃ¼rev Ã§arpÄ±mlarÄ±yla karÅŸÄ±laÅŸtÄ±rÄ±lÄ±r; Ã¶lÃ§Ã¼m verisi veya fiziksel model deÄŸildir.

Yedi seÃ§ilmiÅŸ ek kusur iÅŸaret Ã§arpÄ±mÄ±nÄ±, uÃ§ sÄ±nÄ±rlamasÄ±nÄ±n atlanmasÄ±nÄ±,
doÄŸrusal katsayÄ± ve polinom teriminin erken yuvarlanmasÄ±nÄ±, tensor gradyan
iÅŸaretini, karÄ±ÅŸÄ±k tÃ¼revlerin atlanmasÄ±nÄ± ve kÃ¼bik politika denetiminin
atlanmasÄ±nÄ± hedefler. YalnÄ±z assertion baÅŸarÄ±sÄ±zlÄ±ÄŸÄ± yakalama sayÄ±lÄ±r.
GÃ¼ncel kaynak baÄŸlÄ± tam sonuÃ§lar `output/assurance/verification_summary.json`
iÃ§indedir; seÃ§ilmiÅŸ mutasyonlar eksiksiz kusur kapsamÄ± deÄŸildir.

Ä°lk tam koÅŸuda darbe saldÄ±rÄ± modÃ¼lÃ¼ 120 saniyelik sÄ±nÄ±rÄ± aÅŸtÄ±; kaynaklar ve
eksik sonuÃ§ `history/pchip_first_full_run/archive_manifest.json` iÃ§inde korunur.
TanÄ± yÄ±ÄŸÄ±nÄ± eÄŸri katsayÄ±larÄ±nÄ±n her sorguda yeniden rasyonel kurulmasÄ±nÄ±
gÃ¶sterdi. Ã–n derlenen tamsayÄ± katsayÄ±larÄ± aynÄ± polinomu korur; iÃ§ kÃ¶k/bÃ¼tÃ§e
aramalarÄ± kullanÄ±lmayan tÃ¼revi hesaplamaz, son `PulseBit` iki tÃ¼revi de verir.
1000 seÃ§ilmiÅŸ sorguda Ã¶nceki rasyonel hesapla deÄŸer/tÃ¼rev float'larÄ± aynÄ±ydÄ±;
profil genel gecikme garantisi deÄŸildir. Test girdileri, doÄŸruluk toleranslarÄ±,
200 adÄ±mlÄ±k arama Ã¼st sÄ±nÄ±rlarÄ±, dÃ¶rt iÅŸÃ§i ve 120 saniye sÄ±nÄ±rÄ± deÄŸiÅŸtirilmedi.

[DonanÄ±m olgunluk dosyasÄ±ndaki](PROPULSION_ASSURANCE.md) dokuz kapÄ± aÃ§Ä±ktÄ±r.

## Derlenmiş tamsayı oranlarının dışa aktarımı

Native Lunaris düzeltmesinin ilk tam koşusunda darbe saldırı modülü 120 saniyelik
sınırı aştı. Başarısız koşunun girdi baytları ve mevcut kayıtları
`history/native_lunaris_first_full_timeout/archive_manifest.json` içinde tutulur.
20 saniyelik tanı örneği, derlenmiş eğrinin son `Fraction` kurulumunu gösterdi;
bu örnek tek başına bütün sürenin nedensel dağılımını kanıtlamaz.

Derlenmiş pay/payda yeniden EBOB ile sadeleştirilmeden aynı tamsayı bölmesiyle
float'a çevrilir. Pozitif payda önkoşulu, taşma ve sıfır olmayan sonucun sıfıra
kaybolması denetimleri korunur. Bu, mevcut rasyonel polinomu veya şema/politika
anlamını değiştirmez. Eski çekirdekle 3022 darbe benzeri/işaretli sorguda değer,
türev ve yalnız-değer çıktıları bit düzeyinde aynı kaldı. Beş iç içe zaman
ölçümünün medyanında çekirdek yaklaşık 3.38 kat hızlandı; genel gecikme veya
başka Python sürümleri için garanti verilmez.

Yirminci yöntem, çok büyük ortak çarpanla sadeleştirilmemiş katsayıların
bağımsız Decimal Bezier referansı, işaret ve uç ölçekler karşısında aynı
çıktıları vermesini denetler. Ek seçilmiş kusur, yeni çıktı yolundaki gerçek
alt taşmanın sessiz sıfıra çevrilmesini hedefler. Önceki terim-yuvarlama kusuru
yeni ifadeye taşındı; yeni kusur diye tekrar sayılmadı. İş yükü, dört işçi ve
30/120 saniyelik sınırlar aynı kaldı; tamamlanmış yeni tam koşu ayrıca gerekir.

## References

- [SciPy 1.14.1 PCHIP tanÄ±mÄ± ve algoritma kaynaklarÄ±](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.interpolate.PchipInterpolator.html)
- [SciPy 1.14.1 kÃ¼bik interpolasyon kaynak kodu](https://raw.githubusercontent.com/scipy/scipy/v1.14.1/scipy/interpolate/_cubic.py)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
