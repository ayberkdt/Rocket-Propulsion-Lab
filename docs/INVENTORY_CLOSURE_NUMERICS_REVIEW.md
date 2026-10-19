# Tank envanteri ve kütle kapanışı kararının incelenmesi

Kapsam köprü, elektrikli model ve throttle modelinin `mass_closure`
yardımcıları ile ortak `PropellantMassClosure` kaydıdır. Bu, verilen binary
girdilerin sayısal hesabıdır; donanım envanteri ölçümü veya fiziksel geçerleme
değildir. İmzalı deneme durumları olay kökü aramalarındaki muhasebeyi sınar;
negatif fiziksel tank kütlesine izin verildiği anlamına gelmez.

## Yeniden üretilen kusurlar

Araç kütlesi `1e16+2`, diğer kütle `1` ve tank kütlesi `1e16` olduğunda
tam artık `1 kg` olur. Eski köprü ve elektrikli yardımcılar sıfır raporluyor,
`0.1 kg` toleransında yanlış kapanış veriyordu. Throttle yardımcısı bu örneği
zaten doğru hesaplıyordu; onun ara toplamın erken yuvarlanmasına karşı
kontrolü de ortak çekirdeğe taşındı.

Üç modelin kapanış kararı ise dışa verilen, yuvarlanmış artığa dayanıyordu.
Tam artık toleransı az aşarken gösterilen float artık toleransa eşit
olabiliyor. Testler aynı görünen beş skaler alanın farklı doğru kararlara
sahip iki durumu ayırt etmesini ister. Sınırdaki ve sınırın altındaki geçerli
kapanışlar da başarılı kalmalıdır.

## Hesap ve kayıt sözleşmesi

Seçilmiş tankların toplamı, artık ve göreli artık rasyonel defterde hesaplanır;
her skaler alan sonunda bir kez denetimli float'a dönüştürülür. Köprüde bir
tanka bağlı birden çok akış tankı yalnız bir kez sayar. Bağlanmamış ek
durumlar hesaba katılmaz. Gerçek aralık ihlalleri alan adıyla `DomainError`
verir; sonlu sonucun kurtulduğu ara taşmalar girdi reddine dönüşmez.

Kayıtta yeni isteğe bağlı altıncı alan `residual_ratio` bulunur. Pay ve payda
tam sayı metinleri olarak saklanır; strict JSON dönüşümü sayısal bilgiyi
korur. `is_closed` verilen toleransın tam binary değeriyle bu oranı karşılaştırır.
Kayıt eşitliği oranı da kapsar. Oranı oluşturmak için sonlu bir dış float
sınırı gerekmez: tam artık en büyük float'ı az aşsa bile gösterilen alan sonlu
kalabilir ve bütün sonlu toleranslara karşı doğru karar `False` olur.

Beş alanla elle oluşturulan eski kullanım desteklenir; eksik oran gösterilen
artığın nokta değerinden türetilir. Bu, kaybolmuş asıl envanteri yeniden
oluşturmaz. Bütün skaler alanlar sonlu gerçek değerler olarak denetlenir;
araç kütlesi pozitif, diğer kütle negatif olmayan değerdir. Oranın gösterilen
artığa ve göreli artığa uygunluğu denetlenir, değişebilir liste kopyalanır.
Bu kontroller ölçüm kaynağını veya özgün defteri doğrulamaz.

Üç manifest `exact_binary_inventory_residual_and_ratio_decision_v2`
politikasını bağlar; model kimlikleri değişir. Bilinen politika alanının
silinmesi veya değiştirilmesi okuyucuda reddedilir.

## Bağımsız sayısal kontroller

On üç yöntem, 1600/2000 basamaklı bağımsız Decimal toplam/çıkarım ve bölme
sonuçlarının aynı float'a yerleşmesini doğrular. Artık oranı Decimal'ın
tam sayı oranıyla ayrıca karşılaştırılır. Testler iptal, toplam sırası,
paylaşılan tank, dışa verilen alanların aralık sınırları, yanlış kapalı
kararı, gerçek sınır kabulü, bozuk tolerans/kayıt, strict JSON ve politika
kimliğini kapsar. Alt kontroller ayrı yöntem sayılmaz.

Beş ek seçilmiş kusur yuvarlanmış artıkla kararı, oranın erken yuvarlanmasını,
boolean tolerans kabulünü, kayıt doğrulamasının atlanmasını ve ortak tankın
akış başına tekrar sayılmasını hedefler. Daha önceki iki envanter kusuru
ortak çekirdeğe taşındı; yeni kusur diye tekrar sayılmadı. Yalnız assertion
başarısızlığı yakalama sayılır; ilişkisiz istisna başarılı deney değildir.

Arşivlenmiş eski paketle yeniden oynatma `inventory_historical_replay.json`
içindedir. Yeni sözleşme alanını eski pakette bulamayan testler aritmetik
kusurun bağımsız yeniden üretimi diye sunulmaz. Güncel kaynaklara bağlı tam
sonuçlar `output/assurance/verification_summary.json` içindedir; seçilmiş
testler evrensel doğruluk ispatı değildir. [Olgunluk dosyasındaki](PROPULSION_ASSURANCE.md)
dokuz donanım kapısı açık kalır.

## References

- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
