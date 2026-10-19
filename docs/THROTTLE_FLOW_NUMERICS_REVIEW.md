# Gaz kelebeği debi ve kütle kapanışı incelemesi

Bu inceleme `ThrottleableEngineModel` içindeki debi defterini, tank paylarını,
komut/sürücü debi türevlerini, vakum özgül itkisini ve `mass_closure` artıklarını
kapsar. Eğri örnekleyicisinin verdiği sonlu değerler ve eğimler giriş kabul edilir.
Eğri çekirdeği, bütün motor davranışları ve fiziksel donanım doğruluğu için
tamamlanmış doğrulama veya geçerleme değildir. Bütün örnekler sentetiktir.

## Yeniden üretilen kusurlar

İlk 13 yeni yöntemin dokuzu eski uygulamada başarısız oldu; alt kontroller
ayrı yöntem sayılmaz. `output/assurance/throttle_flow_before_fix.xml` kaydı:

- Pozitif toplam debi sıfıra yuvarlanınca özgül itki bölümünde `ZeroDivisionError`.
- Araç tüketimi sıfırdan farklıyken iki pozitif tank tüketiminin de sıfır olması.
- Komut türevinde son sonuç temsil edilebilirken ara çarpımın sıfır veya sonsuz olması.
- İzin verilen birlik toplamı toleransının tank/araç defterinde kalıcı açık bırakması.
- Özgül itkinin, daha önce yuvarlanan subnormal debi nedeniyle yaklaşık %0,2 sapması.
- Sonlu/temsil edilebilir artıkların, yuvarlanmış envanter veya kütle farkından kaybolması.
- Temsil edilemeyen debi türevleri ve kapanış alanlarının açık alan hatası vermemesi.

## Defter ve temsil kuralları

Verilen binary debi, eğri katsayısı/eğimi ve parametre sürücüsü rasyonel tutulur.
Araç ve tank tüketimleri ile komut/sürücü türevleri bu ortak defterden hesaplanır.
Özgül itki, ihraç edilen debi float'ından geri hesaplanmaz; tam ara debi kullanılır.
Her çıktı alanı sonunda bir kez float'a çevrilir. Gerçek sıfır kabul edilir;
temsil edilemeyen sıfır olmayan alan, alan adı içeren `DomainError` verir.
Bu kural, kurtarılabilir ara taşmaları gerekçesiz reddetmek için kullanılmaz.

Kullanıcı paylarının toplamı tam binary girdiler üzerinden `1e-12` yakınlıkla
denetlenir. Kabul edilen paylar bu tam toplama bölünür; göreli oranları korunur.
Rated toplam debi değişmez. İç defter tam kapanır; ayrı float çıktılarının
toplanmasında son yuvarlama farkı kalabilir. Bütün girdiler ve durumlar için
bit düzeyinde float toplam kapanışı iddia edilmez. Manifestte
`mass_flow_numerics` bu politikayı belirtir ve model hash'i değişir.

Kütle kapanışı, izlenen tankların tam toplamını kullanır. Envanter, artık ve
göreli artık ayrı ayrı ancak son noktada dönüştürülür. Böylece gösterilen
envanter toplamına sığmayan küçük bir katkı, artık hesabında korunabilir.
`mass_closure_numerics` bu anlamı model kimliğine bağlar. Negatif matematiksel
deneme durumlarına fiziksel kabul verilmez; mevcut durum kabul alanı bu
aritmetik düzeltme ile değiştirilmemiştir.

## Kayıttaki sayısal politika alanları

Ek kontrolde üç modelin okuyucusunun bilinen `*_numerics` alanları değiştirilse
veya silinse bile eski model hash'iyle dosyayı kabul ettiği yeniden üretildi.
`numerical_policy_before_fix.json` ve `.xml` kayıtları bu hatayı korur.
Köprü, gaz kelebeği ve elektrikli model okuyucuları artık tanıdıkları sayısal
politikaları yeniden kurulan modelle karşılaştırır ve uyuşmazlığı alan adıyla
reddeder. Bu, bütün bilinmeyen JSON alanları, bütün meta veriler veya yazar
kimliği için tamlık/özgünlük denetimi iddiası değildir.

## Bağımsız kontroller ve seçilmiş kusurlar

Toplam 14 yeni yöntem vardır. Debi yolları ayrı Decimal koduyla kurulur;
komut ve sürücü türevleri iki taraflı farkla hesaplanır. 160/240 basamaklı
referansların float sonuçları eşleşmeden uygulamayla karşılaştırılmaz.
Doğrusal eğri örneklerinde referans uç noktalar arasında bağımsız doğrusal
interpolasyon kullanır; üretimin kübik değerlendiricisini çağırmaz.

Diğer kontroller birlik toplamına yakın payları, sıfır payı, gerçek çıktı
aralığı hatalarını, kapanışta sıra değişimini ve küçük artık iptalini kapsar.
Zamanla değişen komut, RK4 ile entegre edilir; toplam/tank tüketimi analitik
komut integraliyle ve hız roket denklemiyle karşılaştırılır. Ayrıca üç modelde
altı bilinen sayısal politika için değiştirme/silme saldırıları sınanır.

Yedi ek seçilmiş kusur erken komut çarpımı yuvarlamasını, pay normalizasyonunun
atlanmasını, özgül itki öncesi debi yuvarlamasını, kapanışın iki farklı erken
yuvarlamasını, komut türevinden sürücünün düşmesini ve politika denetiminin
atlanmasını hedefler. Yalnız assertion başarısızlığı yakalama sayılır; ilişkisiz
istisna sayılmaz. Bu küme eksiksiz mutasyon kapsamı değildir.

Güncel kaynak bağlı sonuçlar `output/assurance/verification_summary.json`
içindedir. [Donanım olgunluk dosyasındaki](PROPULSION_ASSURANCE.md) dokuz kanıt
kapısı açıktır; bu testler donanım TRL 5 sonucunu değiştirmez.

## References

- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
- [NASA transient-modeling overview](https://ntrs.nasa.gov/citations/20040000363)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
