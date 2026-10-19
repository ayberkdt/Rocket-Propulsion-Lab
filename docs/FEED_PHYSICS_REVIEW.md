# Besleme fiziği: bağımsız sayısal doğrulama

Bu inceleme regülatör, gaz enerji denklemleri ve iki durumlu sıvı hattının
uygulamasını hedefler. Kullanıcı tasarımın ve görev gereksinimlerinin henüz
belirlenmediğini bildirdi. Fiziksel doğrulama alanı ve donanım test nesnesi açık
kalır. Sentetik kontroller hot-fire veya soğuk akış deneyinin yerine geçmez.

## Bulgular ve düzeltmeler

| Kimlik | Bulgu | Düzeltme / kontrol |
|---|---|---|
| FEED-01 | Yakın basınçlarda kuvvet farkı debiyi kaybettiriyor; gamma→1 sınırında kritik oran yanlış yuvarlanıyor | Basınç farkından log1p, kuvvet farkından expm1; bağımsız 75 basamaklı oracle |
| FEED-02 | Sonlu R/gamma girdilerinden sonsuz cp/cv üretilebiliyor | Türetilmiş ısı sığaları için pozitiflik ve sonluluk denetimi |
| FEED-03 | L*C çarpımı taşınca/sıfırlanınca sonlu doğal frekans kayboluyor | Ayrı kareköklerle frekans ve sönüm; iki uç ölçekli sınama |
| FEED-04 | Hat hesabı sonlu olmayan ODE ve tanı çıktısı döndürebiliyor | Türev ve sonuç kayıtlarında sonluluk denetimi |
| FEED-05 | q² ara işlemi denetimsiz OverflowError üretebiliyor | Direnç*q*q sırası; temsil dışı talebin DomainError ile reddi |
| FEED-06 | Kapasite karşılaştırmasını geçen kütle için ullage sıfıra yuvarlanabiliyor | Basınç bölmesinden önce temsil edilebilir pozitif hacim denetimi |
| FEED-07 | cp-cv ve cp*Tb-cv*Tt çıkarımları gamma→1 / eşit T sınırında enerji denkliğini bozabiliyor | cp-cv=R özdeşliği ve ayrılmış taşıma/iş terimleri; Decimal birinci yasa oracle'ı |
| FEED-08 | Sonlu ODE ile birlikte sonsuz Jacobian dönebiliyor | Hat Jacobian sütunlarında sonluluk ve iki bileşenli şekil denetimi |

İlk beş hata sınıfının düzeltme öncesi koşusu `output/assurance/feed_before_fix.xml`
içindedir. Ullage ve Jacobian kusurları ayrıca başarısız regresyonlarla yeniden
üretildi. Isıl çıkarım, eski cp-cv işlemini yeniden enjekte eden mutasyonla da
sınanır. Güncel kanıt `feed_after_fix.xml` ve tam regresyonun modül XML'leridir.
Testleri geçirebilmek için sayısal toleranslar gevşetilmedi.

## Oracle ve sınırlar

Regülatör oracle'ı NASA Glenn ideal gaz/izentropik kütle akışı bağıntılarından
türetilen basınç formunu 75 basamaklı Decimal aritmetiğinde hesaplar; üretimdeki
logaritmik düzenlemeyi kopyalamaz. Otuz koşul: gamma `nextafter(1,+inf)`,
1.000001, 1.1, 1.4, 1.667; basınç oranı `nextafter(1,0)`, 0.999999999, 0.999,
0.8, 0.6, 0.3. Yukarı basınç 2 MPa, T=300 K; R, alan, açıklık ve Cd sentetiktir.
Gamma'nın 1'e çok yakın değerleri gerçek gaz önerisi değildir.

Ek kontroller:

- Kapalı, adyabatik ullage için p*V^gamma sabiti; bağımsız RK4 ile 20/40/80
  adımda hata azalması. Adımlar testin valf zaman ölçeğine uygun seçilir.
- Şişe boşalmasında T/m^(gamma-1) sabiti ve toplam basınçlandırıcı kütlesi;
  sıfır olmayan akışla sınanır.
- Açık sistem birinci yasasıyla eşit sıcaklıkta yüksek hassasiyetli enerji
  hesabı; cp=cv+R yalnızca ideal gaz modelinde kullanılır.
- Sabit yukarı basınç ve sıfır motor talebinde S=(L*q²+C*(p-p_up)²)/2 için
  dS/dt=-R*abs(q)^3; ters akış dahil iki işaret kontrol edilir. S doğrudan joule
  diye etiketlenmez; sabit yoğunlukla fiziksel enerjiye ölçeklenebilir.

Bu hata sınırları donanım kabul toleransı değildir. Test edilmiş örnekler
dışındaki bütün pozitif float girdiler için doğruluk iddia edilmez. Temsil
edilemeyen sonuçlar açık hata üretir. Gerçek gaz, akışkan özelliklerinin
sıcaklık bağımlılığı, valf karakterizasyonu, kavitasyon, su darbesi, çözünmüş
gaz ve tabakalaşma için bu incelemede deney kanıtı yoktur.

## Kusur enjeksiyonu ve kanıt güncelliği

`scripts/run_feed_mutations.py` ilk incelemede kaynakları değiştirmeden sekiz kusur enjekte
eder: debi kazancı, boğulma rejimi, ısı sığası denetimi, modal frekans, şişe
entalpi işareti, ullage sınır işi, ters akış sürtünmesi ve cp-cv çıkarımı.
Sekizi de assertion başarısızlığıyla yakalandı; ilişkisiz istisna başarı sayılmaz.
Bu seçilmiş küme bütün kusurları kapsamaz.

Sonraki [bağlı model incelemesi](COUPLED_FEED_REVIEW.md), negatif olay deneme
durumlarındaki tank/hat kütle açığını ve yanlış sıfır hata raporunu ele alır.
Güncel küme ilk sekiz fizik kusuruna ek olarak üç olay-uzantısı ve iki API
sınırı, dört karışım-türevi/çıktı ve beş köprü defteri kusuru içerir: toplam 22 seçilmiş kusur. Bu sayı bütün giriş
uzayının veya bütün yazılım kusurlarının kapsamı değildir.

Özet rapor artık sabit test sayılarıyla yazılmaz. Regresyon aracı toplanan test
kimliklerini ve JUnit sonuçlarını; kaynak, test, araç, belge ve veri girdilerini
SHA-256 ile bağlar. Yerel Sidera/Lunaris kaynakları da kaydedilir. Özetleyici
değişmiş dosyayı, kayıp testi veya başarısız alt testi başarılı raporlayamaz.
Bu davranış gerçek geçici dosyaların değiştirilmesi ve JUnit kayıtlarının
eksiltilmesiyle ayrıca sınanır. Hash deneyin gerçekliğini/onayını kanıtlamaz.

PowerShell, worktree kökünde:

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) 'src') + ';' + (Get-Location)
python -m pytest tests/test_feed_physics_adversarial.py tests/test_burn_feed_system.py tests/test_burn_feed_line.py tests/test_burn_bipropellant.py -q
python scripts/run_feed_mutations.py
python scripts/run_bounded_regression.py --include-deferred --sidera-path D:\sidera\src
python scripts/summarize_verification.py
```

Yerel Sidera yolu yalnızca bu bilgisayarın açık test yapılandırmasıdır; üretim
kodu bu dizini aramaz. Güncel özet `output/assurance/verification_summary.json`
içindedir. Donanım kapıları: [ana dosya](PROPULSION_ASSURANCE.md).

## References

- [NASA Mass Flow Choking](https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html)
- [NASA Specific Heats](https://www.grc.nasa.gov/www/k-12/airplane/specheat.html)
- [NASA SP-8112](https://ntrs.nasa.gov/citations/19760015212)
- [NASA-STD-7009B](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
- [NASA feedline dynamics](https://ntrs.nasa.gov/citations/19740028545)
