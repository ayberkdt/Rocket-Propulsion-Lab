# Serbest uzay hesabı ve tablolu yanma defterleri

Bu inceleme, `burn_metrics.py` sayısal çekirdeğini; tablolu segmentlerin
integral/moment, örnekleme, ölçekleme ve L1 üretim yollarını; köprünün bunları
tükettiği eksenel kanal yolunu kapsar. Epoch/sınır bulma algoritmalarının,
tablolama okuyucularının bütün yapısal kurallarının veya fizik modellerinin
eksiksiz doğrulandığı iddia edilmez. Dokuz donanım kanıtı kapısı açıktır.

## Somut kusurlar ve sonuç

Kütle `1e200 kg`, sabit kuvvet `1e-200 N`, debi `1e-150 kg/s` ve süre
`1e150 s` olan yapay sayısal örnek eski hesapta Δv, merkez zamanı ve etkin
egzoz hızını sıfır; iki kuadratür farkını da sıfır veriyordu. Bağımsız
Decimal500/700 kapalı form referansı sırasıyla `1e-250 m/s`, `5e149 s` ve
`1e-50 m/s` verir. Bu örnek gerçek donanım önerisi değildir.

Başka örneklerde sonlu moment, süre karesinin erken taşmasıyla sonsuza;
birkaç alt-normal segmentin temsil edilebilir toplamı sıfıra dönüşüyordu.
Ölçek oranı taşması/kaybı ve L1 akışlarının erken yuvarlanması da üretici
yollarını etkiliyordu. Artık verilen binary girdilerin rasyonel defterleri
segment/kanal ölçeklemesi, toplama ve oran kararları bitene kadar korunur.
Bir tanka ait küçük akışlar komutla ölçeklenmeden önce birleştirilir.

Segmentin doğrudan istenen skaler momenti temsil edilemiyorsa alanı belirtilmiş
hata doğabilir; aynı momentin iç defterde tutulup toplam impulsla bölünmesiyle
elde edilen sonlu merkez zamanı yine hesaplanabilir. Ölçekleme de temel
geçmişin sonlu ham örneklerini hash'e bağlar; önce temsil edilemeyen temel
özet çıktısını dışa aktarmayı gerektirmez.

`HOLD` segmenti gerçek bir eğimi yaklaşık eşitlik diye saklayamaz: başlangıç
ve son değerleri eşit olmalıdır. CSV'deki bilinen `"hold"` yazımı önce enum'a
çevrilir; böylece enum ile aynı kesikli profil üretilir. L1 dönüştürücüsü
kullanacağı sonuç kaydının hash'ini ve girdilerinden yeniden üretimini de denetler.

## İntegral ve sayısal sınır

Her segmentin zaman/debi/kütle hesabı rasyoneldir. Gerekirse zaman aralığı,
bir yaprağın tüketimi başlangıç kütlesinin yarısını aşmayıncaya kadar ikiye
bölünür. Bu, fiziksel süreye keyfî bir alt sınır koymaz. İstenen `n` değeri
1–100 aralığındadır; normal yaprak için normalize edilmiş `2n` noktalı
Gauss–Legendre kuralı kullanılır. Düğüm katkısında önce `F/m` float'a
çevrilmez; süre ve ağırlık da hesaba girdikten sonra izole Decimal96 toplamına
geçilir. `quadrature_points` yapılandırılan rafine kuralın yaprak başına
derecesidir; toplam iş miktarı veya bütün yaprakların örneklendiği iddiası değildir.
`integration_subintervals` yaprak sayısını verir.

Kütle sabitken polinom ivme integralleri ve `F` ile kütlenin tam orantılı
olduğu sabit ivme integralleri doğrudan rasyonel hesaplanır. `F=c q` kimliği
de etkin egzoz hızının tam değerini verir. Bunlar özellikle tam binary çıktı
sınırlarını, yaklaşık kuadratürün son bitine göre yanlış sınıflamayı önler.

Bağımsız sayısal aralık şu pozitif Bernstein açılımına dayanır. Bir yaprakta
yerel zaman `x∈[0,1]`, genişlik `h`, başlangıç kütlesi `m` ve uç debiler
`q0,q1` olsun. Tüketilen kütle oranı `u(x)` için derece iki Bernstein
katsayıları `0, h q0/(2m), r` olur; `r=h(q0+q1)/(2m)≤1/2` ve bütün katsayılar
`[0,r]` içindedir. İtkinin katsayıları da negatif değildir.

`F/(1-u)=Σ F u^k` özdeşliğinde ürünün Bernstein katsayıları pozitiftir.
Derecesi `d` olan polinomun integrali katsayı ortalamasıdır; `x` ile çarpılmış
integralde i'nci katsayının ağırlığı `(i+1)/((d+1)(d+2))` olur. Sonlu toplam
aşağı/yukarı yuvarlanarak hesaplanır. Atlanan kuyruk için üst sınır
`integral(F) r^(K+1)/(1-r)`; zaman momenti kuyruğu için bunun yaprağın son
global zamanı ile çarpımı kullanılır. Kuyruk kuvveti pozitif tekrar çarpımıyla
yukarı yuvarlanır. Üst kuyruk paydası aşağı yuvarlanır.

Böylece Δv ve moment için `lo/hi` aralıkları elde edilir. Pozitif Δv için
merkez zamanı `[moment_lo/dv_hi, moment_hi/dv_lo]` içindedir. Kütle oranının
logaritması da yönlü rasyonel dönüşüm, doğru yuvarlanan Decimal `ln` ve komşu
decimal değerleriyle çevrelenir. Etkin egzoz hızı sınırı bu log aralığını
kullanır. Bu argüman, Python Fraction/Decimal işlemlerinin belirtilen
semantiğine bağlıdır; çalıştırılan C kitaplığının bağımsız biçimsel doğrulaması değildir.

Her yayımlanan skalerin sayısal hata sınırı, kendi gerçek referans aralığına
olan en büyük uzaklığın **dışa** binary64 yuvarlanmasıdır. Son çıktı yuvarlaması
da dahil edilir; pozitif çok küçük bir sınır sıfıra silinmez. Alanlar:

- `delta_v_error_estimate_m_s`: Δv sayısal sınırı; eski iki-kural farkı değildir.
- `centroid_error_bound_s`: ağırlıklı zamanın sayısal sınırı.
- `effective_exhaust_velocity_error_bound_m_s`: tüketim varsa egzoz hızının sınırı.
- `error_estimate_kind`: yönlü aralık ve çıktı yuvarlama sözleşmesini belirtir.

Düşük kural derecesi geniş bir sınır verebilir. Bu sınırlar sensör hatası,
kalibrasyon, fiziksel model yetersizliği, yörünge indirgeme kaybı veya görev
kabul eşiği değildir. Native Lunaris kaydı şu an Δv sınırını taşır; diğer
sınırların bütün dış formatlara yayıldığı iddia edilmez.

## Temsil sınırları ve kanıt kapsamı

Aralık binary64 sıfır/sonlu/taşma sınırını keserse referans hassasiyeti
96'dan 192, 384, 768 ve 1536 basamağa yükseltilir. Hâlâ çözülemeyen bir aralık
`ConvergenceError` üretir; bu fiziksel girdinin geçersiz olduğu iddiası değildir.
Sonlu doğru sonucun kurtarıldığı yarım-subnormal örnek ve tam bağ/tie örnekleri
ayrı test edilir. Temsil edilemediği belirlenen gerekli skaler açık hata verir.
Evrensel doğru-yuvarlama veya bütün kabul edilebilir girdiler için genel süre
garantisi verilmez.

`integral_numerics=exact_binary_segment_aggregate_ledger_v1` artifact kimliğine;
`tabulated_state_numerics=exact_binary_axial_and_grouped_flow_channels_v1` köprü
kimliğine dahildir. Eksik veya değiştirilmiş tanınan politika reddedilir.
Eski belgeler kaynak örneklerinden yeniden üretilmelidir; yeni etiketi eski
sonucun üzerine elle koymak doğrulama sayılmaz. Bu, bütün Tabulated/CSV
okuyucularının alan şeması veya ham metin saldırıları için tamamlanmış denetim değildir.

38 yeni yöntem; bağımsız polinom/log kapalı formları, 100–1600 basamakta
referans yakınsamasını, artan/azalan debiyi, düşük/yüksek kuadratür derecesini,
doğrudan aralık/tail eşitsizliklerini, sahte sonuçları, ölçekleme ve gerçek çıktı
aralığı hatalarını sınar. 21 seçilmiş kusur bu yöntemlerle ilişkilidir. Yalnız
assertion başarısızlığı yakalama sayılır. Önbellekli işlevin kaynak bağlamı
`inspect.unwrap` ile alınır; her varyantta saf hesap önbellekleri temizlenir.

Köprü regresyonu, kamuya açık bütün örnek alanlarını dışa aktarmanın burada
gereksiz bir alt-normal akışı erken reddettiğini gösterdi. Köprü artık yalnız
gereken eksenel alanı iç defterden alır; toplu debiler kendi tam yolunda kalır.
Eski kontrolün assertion'ları gevşetilmedi. Bir kuadratür kusurunu aralık
çözücüsü düzeltebildiğinde, test ayrıca temel kuadratür sonucunu bağımsız
kapalı forma karşı denetler; aksi halde bu sağ kalan kusur yakalanmış sayılmaz.

Tamamlanmamış ilk birleşik entegrasyon koşusu süreç kimliği/CPU ölçümüyle
durduruldu ve ayrı kaydedildi. Tekrarlanan değişmez toplamın önbellek deneyi
beş iç içe ölçümde aynı rasyonel sonuçla yaklaşık 21 kat mikroçekirdek oranı
verdi; genel gecikme garantisi değildir. Son doğrulamada dört işçi,
30 saniye koleksiyon ve 120 saniye/modül sınırları korunur.

Güncel sonuç `output/assurance/verification_summary.json` içindedir.
`tabulated_boundary_input_followup_screening.json`, kısa profillerde zaman
sınırı sahipliği/snap ve ham alan kapsamı için ayrıca açık kalan gözlemleri
tutar. Sonraki [zaman incelemesi](TABULATED_TIME_CONTRACT_REVIEW.md) göreli
sahipliği ve kanonik epoch olaylarını ele alır; ham belge sözleşmeleri ile
kapsamlı fonksiyon/katalog/iddia incelemesi tamamlanmış değildir.

## References

- [NIST DLMF 3.5(v), Gauss quadrature](https://dlmf.nist.gov/3.5#v)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html): izole bağlam, yuvarlama ve doğru yuvarlanan `ln`.
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html): verilen binary girdilerin rasyonel defteri.
- [Donanım kanıtı ve olgunluk sınırları](PROPULSION_ASSURANCE.md)
