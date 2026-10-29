# Sentinel serileri — hesap türetimleri ve doğrulama (V4)

**29 Ekim 2026.** Kaynaklı büyüklükler, kataloglar ve varsayımlar ayrıdır. Hesaplar
yeni donanımın yeterliliğini, tam uçuş bütçesini veya TRL'sini kanıtlamaz. Önceki V3'ün
"bütün sayılar hatasız" genel hükmü yerine kapsamı belirli kontroller kullanılır.
Kaynak kodları: [KAYNAKLAR.md](KAYNAKLAR.md).

## 1. Girdiler ve yeniden üretim

Python 3.11+ standart kütüphanesi yeterlidir. Depo kökünden:

```powershell
python scripts/sentinel_series_budget.py --output docs/missions/sentinel/SERI_SONUCLARI.json
python -m unittest discover -s tests -p test_sentinel_series_budget.py
```

Veri: [SERI_VERILERI.json](SERI_VERILERI.json); sonuç: [SERI_SONUCLARI.json](SERI_SONUCLARI.json).
Kod: [sentinel_series_budget.py](../../../scripts/sentinel_series_budget.py);
kontroller: [test_sentinel_series_budget.py](../../../tests/test_sentinel_series_budget.py).
SI birimleri: kg, m, s, N, m/s. 1 cm/s = 0,01 m/s; JSON şeması
`sentinel-series-v4-SI`. Farklı şema reddedilir; yanlış birimde yazılmış etiketsiz
bir sayıyı otomatik tanımak mümkün değildir.

| Sabit/girdi | Değer | Kapsam |
|---|---:|---|
| g0 | 9,80665 m/s² | Isp standart yerçekimi |
| μ | 3,986004418×10¹⁴ m³/s² | Dünya iki cisim modeli |
| RE | 6378137 m | Ekvator yarıçapı; küresel yaklaşım |
| ωE | 7,2921159×10⁻⁵ rad/s | Yer izi lineerleştirmesi |
| S1C başlangıç / hidrazin | 2186 / 154 kg | VV25 kiti / ESA dolumu |
| S2C başlangıç / hidrazin | 1143 / 133 kg | VV24 kiti / ESA dolumu |

A/B/D için eşleşmiş kabul/dolum çifti bu hesap kümesine atanmadı. Kütle ve yakıt araç
kimlikleri uyuşmazsa hesap durur; bu kontrol belge gerçekliğini kendi başına denetlemez.
2032/1010 kg farkları sertifikalı kuru kütle değildir; basınçlandırıcı gazı içerebilir.

## 2. Tsiolkovsky türetimi ve kapasite

Vakum, sabit yerel kuvvet yönü ve dış kuvvetlerden doğan yanma kaybının ihmaliyle
`m dv = −ce dm`, `ce = g0 Isp`. İntegrasyon:

```text
ΔV = ce ln(mi/mf)
mf = mi exp(−ΔV/ce)
mp = mi [1 − exp(−ΔV/ce)]
```

Kod çok küçük ΔV'de iptal hatasını önlemek için `−expm1(−ΔV/ce)`; kapasitede
`−log1p(−mp/mi)` kullanır. Kuvvet ideal toplam ΔV'yi doğrudan belirlemez; süreyi ve
sonlu yanma kayıplarını etkiler.

| Varsayılan sabit Isp | S1C bütün 154 kg ideal kullanılırsa | S2C bütün 133 kg ideal kullanılırsa |
|---:|---:|---:|
| 200 s | — | 242,63 m/s |
| 205 s | 146,86 m/s | — |
| 214 s | 153,31 m/s | — |
| 220 s | 157,61 m/s | 266,89 m/s |
| 223 s | — | 270,53 m/s |

Bu **ideal kapasite**, tahsis değildir. 220 s S1 için duyarlılık varsayımıdır;
katalog >205/>214 alt sınırlarından gerçek C Isp'si çıkarılamaz. Basınç haritası,
erişilemeyen yakıt ve hat kaybı bu tabloda yoktur.

## 3. Eğik itici ve skaler debi

Yakıt skaler, kuvvet vektörel toplanır:

```text
Feksen = Σ di Fi cosθi
ṁtoplam = Σ di Fi / (g0 Ispi)
ce,etkin = Feksen / ṁtoplam
ηeksen = Feksen / ΣdiFi     (eş Isp: ce,etkin = g0 Isp ηeksen)
```

di açık kalma oranıdır. Karşı tork darbeleri, yönler, plume ve yedek dal seçimi
"motor adedi × nominal N" hesabını geçersiz kılar. Roket denklemindeki ce, net eksensel
ΔV'yi **bütün** tüketilen yakıtla eşleştirmelidir.

1 N ve 220 s eşleşmiş ideal noktada `ṁ=0,000463507 kg/s=0,463507 g/s`.
ArianeGroup'un ayrı nominal 0,44 g/s ve 220 s değerleri **0,949284 N** verir;
bu nominal tanımların aynı ölçüm noktası olduğu varsayılmaz.

## 4. Kademeli kütle defteri: S2C yayımlanmış manevraları

SpaceOps #617 s.15–16 yaklaşık edinim değerleri ve Tablo 3 kalibre bileşenleri.
İlk OOP kaleminde makale planlama dili kullanır: toplam, eksiksiz telemetri defteri
değil, **yayımlanmış özetten rekonstrüksiyon**dur.

| Kalem | ΔV | Sayım |
|---|---:|---|
| Dört büyük IP edinim | yaklaşık 6,33 m/s | Toplu |
| 29 Ekim ince edinim | 0,012 m/s | Bir kez |
| İki ilk OOP | yaklaşık 0,57 m/s | Toplu büyüklük |
| Tablo 3: 31 Ekim–18 Aralık | 1,0564418 m/s | 8 olay; 31 Ekim ayrıca eklenmez |
| Toplam | **yaklaşık 7,97 m/s** | Aritmetik 7,9684418; kaynak kesinliği daha düşük |

6 Kasım tek yanma için `sqrt(0,9344²+0,0537²)=0,9359418 m/s`. Farklı zamanlardaki
23 Kasım −0,0207 ve +0,0206 ise **0,0413 m/s** maliyettir; net −0,0001 değildir.
Toplam yakıt maliyetinde `Σ||Δv_i||` kullanılır, `||ΣΔv_i||` kullanılmaz.

Her olayda `mi+1=mi−mpi`. 220 s/η=1 ideal defter: son kütle **1138,7862 kg**, kalan
hidrazin **128,7862 kg**, tüketim **4,2138 kg**. Bunlar **ölçülmüş yakıt değildir**;
gösterilen basamaklar yalnız hesap tekrar üretimi içindir.

| Isp [varsayım] | ηeksen [varsayım] | Hesap tüketimi |
|---:|---:|---:|
| 200 s | 0,8 | 5,7900 kg |
| 200 s | 1,0 | 4,6343 kg |
| 220 s | 0,8 | 5,2648 kg |
| 220 s | 1,0 | 4,2138 kg |

Olay kimliği kontrolü aynı kimliğin tekrarını yakalar. Farklı kimlikle aynı olayı
yeniden eklemeyi önlemek ayrıca kaynak eşleştirmesi ister; şema fiziksel bağımsızlığı
tek başına kanıtlamaz.

## 5. Kuvvet, süre, debi

Sabit ce ve debi yaklaşımıyla `t=mp/ṁ`, `Feksen=mp ce,etkin/t`.
1 m/s, 220 s ve tek net 1 N için:

| Araç | Yakıt | Süre | Başlangıç ivmesi |
|---|---:|---:|---:|
| S1C 2186 kg | 1,0130 kg | 36,42 dk | 0,4575 mm/s² |
| S2C 1143 kg | 0,5297 kg | 19,05 dk | 0,8749 mm/s² |

S2C 6 Kasım modelinde mi=1139,336 kg; mp=0,494155 kg; t=314,1 s;
Feksen≈**3,39421 N**. `mi ΔV/t` küçük kütle değişimini ihmal eder. Gerçek ce/F
zamanla değişebilir; sonuç çıkarılmış ortalamadır. Tek sayıdan cant, basınç ve duty
ayrı ayrı bulunamaz. Dört motor tam 1 N/duty=1 olsaydı cant yaklaşık 32° olurdu;
motor başına yaklaşık 0,85 N ile sıfıra yakın cant da aynı net kuvveti verebilir.

## 6. Yörünge geçişi: plan ve gerçek son durum

İki cisim, aynı düzlem, vis-viva: `v²=μ(2/r−1/a)`. r0 daireselden alçak rp×ra'ya:
önce r0×rp transferine impuls; rp'de apojeyi ra'ya indiren ikinci impuls. İki büyüklük
toplanır. Model J2, atmosfer, sonlu yanma veya gerçek faz programı içermez.

| Geçiş | İdeal ΔV | Anlamı |
|---|---:|---|
| S1 693→530 km dairesel | 88,0576 m/s | B'nin eski hedefi; gerçekleşen son yörünge değil |
| S1 693→580×610 km | 52,5822 m/s | B son yörüngesiyle geometrik kıyas; gerçek burn toplamı değil |
| S2 786→653 km dairesel | 70,2157 m/s | 2014 aile EOL geometrisi |

S2'nin 73 m/s planı idealden **%3,9654** fazladır; kaynak kayıpları ayırmadığı için
tamamı "sonlu yanma marjı" denmez. 116×600=69600 s=19,33 saat planlı ateşleme.
73 m/s ve 47 kg ile ters çözüm, EOL başlangıç kütlesi 1057–1139 kg varsayıldığında
etkin Isp≈164–177 s verir. Cant, blow-down ve marjın payı ayrıca bilinmelidir.

1B'nin yaklaşık 70 kg yakıtla sonlandırılan 2024 bertarafı ve 24 yıllık balistik
öngörüsü [#226 s.20], rezervin tek başına görev ömrü olmadığını gösterir. Kısmi
pasivasyon ve manevra belirsizliğinin CAM taramasına etkisi de önemlidir.

## 7. CONOPS senaryosu ve koşullu ek ömür

S1C örnek **varsayım**: edinim 11,01; rutin 7,25 yıl×1; CAM 3;
B son yörünge geometrisi EOL 52,5822 m/s. Toplam **73,8422 m/s**, 205 s/η=1 ile
**78,8367 kg**. 11,01; ISSFD 2017 planlama kalemlerinin aritmetik toplamıdır,
kovaryanslı birleşik 1σ bütçesi değildir. Bu C'nin kabul edilmiş tahsisi/EOL planı değildir.

S2C'de 133−4,2138−47=81,7862 kg aritmetik kalan örneği; ideal rekonstrüksiyon ile
başka sürüm EOL planını birleştirir. Safe-mode, CAM, gerçek OOP, erişilemeyen yakıt ve
basınç zarfı yoktur. "Yakıt görevi sınırlamaz" sonucu çıkarılamaz.

**Tutarlı koşullu ek yıl:** Bugünkü toplam kütle m; rutin için ayrılabilecek yakıt
`p=p_mevcut−p_EOL−p_emniyet−p_erişilemeyen`. Pozitif p, sabit ce ve sabit yıllık
toplam ΔVyıl varsayımıyla:

```text
N_ideal = ce ln[m/(m−p)] / ΔVyıl
```

EOL yakıtı rutin boyunca taşındığı için m'nin içinde kalır; başta kütleden ayrıca
düşülmez. m=1000 kg, p=0,5 kg, ce=220g0, ΔVyıl=1 m/s için yaklaşık **1,079 yıl**.
Bu belirtilen varsayımların cevabıdır; motor, valf, basınç, batarya ve yük ömürleri
ayrıca sınanır. Değişen yıllık ihtiyaçta yıl yıl yakıt/impuls/çevrim defteri gerekir.

## 8. Basınçlı besleme modeli ve sınırı

İdeal polytropik blow-down `P Vgazⁿ=sabit`; n=1 izotermal, n≈γ adyabatik sınır.
`Vgaz=Vtank−mp/ρp` ancak diyafram/ölü hacim ve yoğunluk doğruysa uygulanır.
`Pgiriş=Ptank−ΔPhat−ΔPvalf−ΔPfiltre`; kayıplar debi, sıcaklık ve viskoziteye bağlıdır.

Örnek **varsayım**: 200 L tank, 150 L sıvı, 20 bar **mutlak**, izotermal. Sıvı bittiğinde
gaz 50→200 L; basınç 5 bar olur. Minimum motor girişi 5,5 bar ise tüm sıvı yararlı
değildir. Bu Rafael/C uçuş hesabı değildir; katalog tank maksimumunu motor minimumuyla
doğrudan karşılaştırmanın neden yeterli olmadığını gösterir.

Gaz kütlesi ilk yaklaşımda `mg=P Vgaz/(Rgaz T)`; gerçek gaz Z, sızıntı/çözünme,
diyafram, ısıl geçmiş ve gerçek kabul basıncı ayrıca gerekir. S2C'nin 21 bar dolumu
ile katalogdaki 5,5–22 bar **motor giriş** aralığı farklı sınır noktalarıdır.

## 9. Tandem kontrolü: yaklaşık fizik karşılaştırması

Küçük teğetsel hata: `δa/a≈2δv/v`; `δn/n≈−3δv/v`.
`|δt_faz|≈3|δv|Δt/v`; dönen Dünya yer izi için `|δXT|≈ωE RE |δt_faz|`.
J2, atmosfer, eş manevra, enlem ve gerçek propagasyon yerine geçmeyen lineerleştirmedir.

| Hata / süre | Model | SpaceOps #617 s.5 yuvarlak örneği |
|---|---:|---:|
| 1 cm/s / 1 gün: yer izi | 161,62 m | 160 m |
| 1 cm/s / 1 gün: zaman | 0,3475 s | 0,342 s |
| 3,85 cm/s / 1 gün: yer izi | 622,24 m | 614 m |

Fark %2'nin altındadır; katsayılar bu sayılara sonradan uydurulmadı. Bu makale tahminiyle
model tutarlılığıdır, bağımsız uçuş yakıt ölçümü veya hot-fire doğrulaması değildir.
cm/s→m/s faktör 100 hatası kontrolü geçemez. Eş ve dönüş manevraları kontrolsüz
drift hesabını değiştirir; gerçek 50 günlük performansa tek hata yayılmaz.

## 10. Güç, veri, link ve ısıl örnekler

`T=2πsqrt(r³/μ)`; β=0 için `t_tutulma=T asin(RE/r)/π`.
S1: 98,6265/35,2982 dk; S2: 100,5786/35,1518 dk (tur/geometrik tutulma).
Gerçek tutulma β ve mevsime bağlıdır. `Ebat≈Vbat Ah DoD η`;
`Eüretim≥Eyüklere+Ekayıp` en kötü mevsim ve mod profiliyle sağlanmalıdır.

Veri: 1410 Gbit=176,25 GB; 2,4 Tbit=300 GB; 520/560 Mbit/s×600 s=39/42 GB;
300 GB/560 Mbit/s=71,43 dk. Bit/byte, GB/GiB, ham/kodlu hız ve kanal adedi ayrıdır.

Link **varsayımı**: 8,2 GHz, 1800 km, 10 W RF, 12 dBi anten, G/T=25 dB/K,
diğer kayıp 3 dB. `Lfs=20log10(4πRf/c)=175,83 dB`;
`C/N0=EIRP+G/T−Lfs−Ldiğer+228,6≈96,77 dBHz`.
`Eb/N0=C/N0−10log10(Rb)`; tek taşıyıcı 260/560 Mbit/s için 12,62/9,29 dB.
Eşik ve uygulama kaybı olmadığından **marj** verilemez; uçuş link bütçesi değildir.

İdeal derin uzay radyatörü `A=Q/(εσT⁴)`: 300 W, ε=0,8, 300 K → **0,81646 m²**.
Güneş/albedo/Dünya IR ve iletim yok sayılmıştır. Gerçek termal kapanış değildir.
Rafael valf 9,2 W@28 V; her 257 Ω ısıtıcı 3,05058 W. Açık/kapalı durum ve duty
bilinmeden 14 ile çarpılıp gerçek modül gücü bulunamaz.

## 11. Kontroller ve doğrulama sınırı

25 test yöntemi: yüksek hassasiyetli Decimal oraklı, çok küçük ΔV, kütle korunumu,
sıralı/tek toplam eşdeğerliği, CAM işaretleri, vektör normu, duplikasyon, yanlış araç
eşleştirme, yakıt yetersizliği, cant, bağımsız açısal momentumla geçiş, yayımlanmış
drift örneği, veri birimleri ve Python `-I -S` ile bağımlılıksız tekrar üretim.
16 alt test yanlış tür/NaN/sonsuz/negatif vb. girdileri sınar. Test sayısı fiziksel
güven veya TRL düzeyi değildir.

Tam sonlu yanma propagasyonu, gerçek motor haritaları, motor başına yaşam defteri,
uçuş EPS/link/termal bütçeleri ve tank telemetrisi model dışında kalır. Yeni
motor+besleme sistemi için temsilî ortam deneyleri halen gereklidir.
