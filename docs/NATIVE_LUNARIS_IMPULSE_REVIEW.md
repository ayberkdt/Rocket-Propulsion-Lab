# Native Lunaris dÃ¼rtÃ¼ normu ve kÃ¼tle gÃ¼ncellemesi

Kapsam `lunaris.core.propagation.plans.apply_impulsive_maneuver` ve norm/kÃ¼tle
yardÄ±mcÄ±sÄ±dÄ±r. Yama, RocketProp'a yeni bir fizik modeli eklemez; mevcut
`m_after = m_before exp(-|dv|/(g0 Isp))` baÄŸÄ±ntÄ±sÄ±nÄ±n sayÄ±sal uygulanÄ±ÅŸÄ±nÄ± dÃ¼zeltir.
Bu inceleme bÃ¼tÃ¼n Lunaris kodunun veya ilgili yÃ¶rÃ¼nge/propulsion modellerinin
geÃ§erlemesi deÄŸildir. UÃ§ Ã¶lÃ§ekler sentetik sayÄ±sal stres girdileridir; gÃ¶rev
Ã¶nerisi veya donanÄ±m test koÅŸulu deÄŸildir.

## Yeniden Ã¼retilen kusurlar

Sonlu yaklaÅŸÄ±k `1e297 m/s` komutunda doÄŸrudan NumPy normu taÅŸÄ±yor; envanter
eÅŸdeÄŸerinin Ã¶ngÃ¶rdÃ¼ÄŸÃ¼ `999.9 kg` yerine native uygulama sÄ±fÄ±r kÃ¼tle veriyordu.
Ã‡ok kÃ¼Ã§Ã¼k vektÃ¶rde de norm sÄ±fÄ±rlanÄ±yor ve sonlu bir kÃ¼tle oranÄ± kayboluyordu.
`g0 Isp` ara Ã§arpÄ±mÄ± ayrÄ±ca taÅŸabiliyor; Ã¼stel tek baÅŸÄ±na float'ta sÄ±fÄ±rlansa
bile baÅŸlangÄ±Ã§ kÃ¼tlesiyle Ã§arpÄ±mda kurtulabilecek pozitif son kÃ¼tle kayboluyordu.
GerÃ§ek Ã§Ä±ktÄ± aralÄ±k ihlalleri de sessiz sÄ±fÄ±r/sonsuz kayÄ±tlara dÃ¶nÃ¼ÅŸebiliyordu.

## Somut yama ve korunan davranÄ±ÅŸ

[Yama kaydÄ±](assurance/lunaris_impulse_patch.json), eski/yeni dosyanÄ±n SHA-256
deÄŸerlerini ve [birleÅŸik farkÄ±](assurance/lunaris_impulse.patch) baÄŸlar. Tam
aday kaynak `scripts/native_lunaris_candidate/plans.py` iÃ§indedir. Uygulama
aracÄ± canlÄ± dosyanÄ±n hash'i inceleme temelinden farklÄ±ysa durur; baÅŸka bir
dosyayÄ± veya kullanÄ±cÄ± deÄŸiÅŸikliÄŸini Ã¼zerine yazmaz. Eski kaynak ve kanÄ±t
kopyasÄ± korunur. GÃ¼ncel uygulama durumu JSON kaydÄ±nda ve
`output/assurance/native_lunaris_apply_result.json` iÃ§indedir.

Pozitif kareler toplamÄ±, karekÃ¶k, `g0 Isp` oranÄ±, Ã¼stel ve baÅŸlangÄ±Ã§ kÃ¼tlesi
Ã§arpÄ±mÄ± Ã§aÄŸÄ±ranÄ±n baÄŸlamÄ±ndan baÄŸÄ±msÄ±z 96 basamaklÄ± Decimal baÄŸlamÄ±nda
hesaplanÄ±r. Norm ve son kÃ¼tle en sonda float'a Ã§evrilir. Bu sabit hassasiyet
seÃ§imidir; evrensel doÄŸru yuvarlama veya ispatlanmÄ±ÅŸ bÃ¼tÃ¼n-Ã§Ä±ktÄ± hata zarfÄ±
iddiasÄ± deÄŸildir. KÃ¼tle oranÄ± Ã¼ssÃ¼ 4096'yÄ± aÅŸarsa en bÃ¼yÃ¼k binary64 ilk
kÃ¼tlede bile son kÃ¼tle temsil edilemez; devasa Ã¼stel hesap yerine aÃ§Ä±k hata
verilir. Daha Ã¶nceki ara taÅŸma/alt taÅŸmalar gerekÃ§esiyle temsil edilebilir
bir sonucun reddi yapÄ±lmaz.

Normun skaler kayÄ±t alanÄ±, pozitif son kÃ¼tle veya gÃ¼ncellenmiÅŸ hÄ±z gerÃ§ekten
binary64 aralÄ±ÄŸÄ±na sÄ±ÄŸmÄ±yorsa `ValueError` oluÅŸur. Ã‡aÄŸÄ±ranÄ±n durum dizisi
baÅŸarÄ±da ve hatada deÄŸiÅŸtirilmez. SÄ±fÄ±r dÃ¼rtÃ¼, Isp olmadan kinematik kullanÄ±m,
RIC dÃ¶nÃ¼ÅŸÃ¼mÃ¼ ve kayÄ±t alanlarÄ±nÄ±n anlamÄ± korunur. DeÄŸiÅŸiklik kÃ¼tle kaydÄ±nÄ±n
gerÃ§ek son durumla aynÄ± olmasÄ±nÄ± da denetler.

## BaÄŸÄ±msÄ±z kanÄ±tlar

15 yeni saldÄ±rgan yÃ¶ntem bÃ¼yÃ¼k/subnormal normlarÄ±, bÃ¼yÃ¼k/kÃ¼Ã§Ã¼k Isp faktÃ¶rlerini,
Ã¼stel alt taÅŸmadan kurtulan normal ve subnormal son kÃ¼tleleri, gerÃ§ek Ã§Ä±ktÄ±
aralÄ±ÄŸÄ± hatalarÄ±nÄ±, RIC dÃ¶nÃ¼ÅŸÃ¼mÃ¼nÃ¼, sÄ±fÄ±r ve geÃ§ersiz durumlarÄ±, sabit tohumlu
24 vektÃ¶r/kÃ¼tle Ã¶lÃ§eÄŸini, strict JSON kaydÄ±nÄ± ve RocketProp'tan native uygulamaya
ilk kusurun aktarÄ±mÄ±nÄ± kapsar. Alt kontroller ayrÄ± yÃ¶ntem sayÄ±lmaz.

Referans, kareler toplamÄ±nÄ± Fraction Ã¼zerinden kurar; 200/300 basamaklÄ±
karekÃ¶kler aynÄ± float'a yerleÅŸmeden karÅŸÄ±laÅŸtÄ±rma yapÄ±lmaz. Ãœstel referansÄ±
Ã¼retimin Decimal `exp` Ã§aÄŸrÄ±sÄ±nÄ± kopyalamaz: kÃ¼Ã§Ã¼k argÃ¼mana indirgeme, Taylor
serisi ve yinelenen kare alma kullanÄ±r. BileÅŸen/kÃ¼tle karÅŸÄ±laÅŸtÄ±rmasÄ± dÃ¶rt
ulp'dir; RocketProp eÅŸdeÄŸerinin mevcut yaklaÅŸÄ±k integralinden gelen son
kÃ¼tle kontrolÃ¼ sekiz ulp'dir. Bu sayÄ±sal kanÄ±t fiziksel belirsizlik deÄŸildir.

943 dosyalÄ± ayrÄ± Lunaris kopyasÄ±nda aday Ã¶nce 40 ilgili test ve kaynaÄŸÄ± hash'le
baÄŸlanan beÅŸ upstream manevra testiyle denetlendi. CanlÄ± kaynak bu hazÄ±rlÄ±kta
deÄŸiÅŸtirilmedi. Son uygulama ve tam regresyon kanÄ±tlarÄ± ayrÄ±ca kaydedilir;
yalnÄ±z adayÄ±n baÅŸarÄ±lÄ± olmasÄ±, canlÄ± baÄŸÄ±mlÄ±lÄ±ÄŸÄ±n dÃ¼zelmiÅŸ sayÄ±lmasÄ± iÃ§in
yeterli deÄŸildir.

AltÄ± seÃ§ilmiÅŸ native kusur norm kazancÄ±nÄ±, erken `g0 Isp` float Ã§arpÄ±mÄ±nÄ±,
Ã¼Ã§ Ã§Ä±ktÄ± aralÄ±ÄŸÄ± denetiminin kaldÄ±rÄ±lmasÄ±nÄ± ve kÃ¼tle deÄŸiÅŸiminin kayÄ±tta
gizlenmesini hedefler. YalnÄ±z assertion baÅŸarÄ±sÄ±zlÄ±ÄŸÄ± yakalama sayÄ±lÄ±r.
Native mutasyon kaydÄ± hem RocketProp girdilerine hem yÃ¼klenen Lunaris Python
kaynaklarÄ±na baÄŸlanÄ±r. BaÄŸÄ±mlÄ±lÄ±k deÄŸiÅŸtiÄŸinde eski yakalamalar kullanÄ±lamaz;
genel doÄŸrulama Ã¼Ã§ ayrÄ± mutasyon kaydÄ±nÄ±n tamamÄ±nÄ± ister. Native mutasyonun
kaynak kimlikleri tam regresyondaki baÄŸÄ±mlÄ±lÄ±k kimlikleriyle de eÅŸleÅŸmelidir;
baÅŸka bir kopyanÄ±n geÃ§erli kanÄ±tÄ± bu koÅŸuyu onaylayamaz.

CanlÄ± kaynak taramasÄ± da eski snapshot ile aynÄ± `.git`, `__pycache__` ve
`node_modules` dÄ±ÅŸlamalarÄ±nÄ± kullanÄ±r. BÃ¶ylece UI derleme baÄŸÄ±mlÄ±lÄ±klarÄ±nÄ±
tarama/izleme native Python kanÄ±tÄ± gibi sunulmaz; ilgili Python kaynaÄŸÄ±ndaki
ekleme ve iÃ§erik deÄŸiÅŸiklikleri hÃ¢lÃ¢ yakalanÄ±r. Ä°lk canlÄ± denetim giriÅŸimi bu
kapsam hatasÄ±nÄ± dÃ¼zeltmek iÃ§in durduruldu ve ayrÄ± kaydedildi; tamamlanmÄ±ÅŸ
doÄŸrulama gibi sayÄ±lmadÄ±.

GÃ¼ncel kaynak baÄŸlÄ± tam sonuÃ§ `output/assurance/verification_summary.json`
iÃ§indedir. Eski sabitlenmiÅŸ kaynak Ã¼zerine yazÄ±lmaz; yeni doÄŸrulama iÃ§in
ayrÄ± kaynak/kaynak-verisi kopyasÄ± Ã¼retilir. Kurulu ikili ortam ve fiziksel
donanÄ±m kopyalanmÄ±ÅŸ veya nitelendirilmiÅŸ sayÄ±lmaz.
[Olgunluk dosyasÄ±ndaki](PROPULSION_ASSURANCE.md) dokuz donanÄ±m kapÄ±sÄ± aÃ§Ä±ktÄ±r.

## References

- [Python 3.12 Decimal](https://docs.python.org/3.12/library/decimal.html)
- [Python 3.12 Fraction](https://docs.python.org/3.12/library/fractions.html)
