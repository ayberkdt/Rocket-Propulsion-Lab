# Sentinel-1C ve Sentinel-2C — hesaplar (Revizyon V3)

Tüm sayılar bağımsız olarak yeniden hesaplandı. V2 paketindeki her sayı (ideal ΔV tabloları,
Hohmann, perigee düşürme, yakıt, kütle debisi, yanma süresi, tutulma, link, veri) birebir
doğrulandı. V3, keyfî senaryo değerlerinin yerine kaynaklı uçuş değerlerini koyar.
Kaynak kodları: [KAYNAKLAR.md](KAYNAKLAR.md).

## 1. Sabitler ve girdiler

| Değişken | Değer | Not |
|---|---|---|
| $g_0$ | 9,80665 m/s² | Isp tanımındaki standart yerçekimi |
| $\mu$ | 398600,4418 km³/s² | Dünya kütle çekim parametresi |
| $R_E$ | 6378,137 km | WGS84 ekvator yarıçapı; J₂ ihmal |
| S1C $m_0$ / itici | 2186 / 154 kg | [D] VV25 kiti; ESA |
| S2C $m_0$ / itici | 1143 / 133 kg | [D] VV24 kiti; ESA |
| S1 Isp aralığı | > 205 s (5,5 bar) … > 214 s (22 bar) | [M] Rafael 1 N kataloğu |
| S2 Isp | nominal 220 s, aralık 200–223 s | [M] ArianeGroup 1 N veri sayfası |

## 2. Tsiolkovsky

$v_e = I_{sp}g_0$, $\Delta V = v_e\ln(m_i/m_f)$, $m_p = m_i\left(1-e^{-\Delta V/v_e}\right)$.

Bütün itici ideal kullanılırsa (yalnızca kapasite, tahsis değil):

| Isp | S1C (2186 → 2032 kg) | S2C (1143 → 1010 kg) |
|---:|---:|---:|
| 200 s | — | 242,6 m/s |
| 205 s | 146,9 m/s | 248,7 m/s |
| 214 s | 153,3 m/s | 259,6 m/s |
| 220 s | 157,6 m/s | 266,9 m/s |
| 223 s | — | 270,5 m/s |

`m_0 − itici` farkı (2032 / 1010 kg) doğrulanmış kuru kütle değildir (basınçlandırıcı gaz ve
artık sıvılar dahil olabilir).

## 3. Hohmann ve ömür sonu geçişleri

$v_c=\sqrt{\mu/r}$, $a_t=(r_1+r_2)/2$,
$\Delta V=|v_t(r_1)-v_c(r_1)|+|v_c(r_2)-v_t(r_2)|$.

| Geçiş | ΔV (m/s) | Not |
|---|---:|---|
| S1C 700 → 693 km | 3,71 | yalnızca karşılaştırma; gerçek enjeksiyon elemanları bilinmiyor |
| S2C 775 → 786 km | 5,73 | gerçekte 2C nominalin 10,6 km altına bırakıldı ve sürüklenme ile faz aldı; bu bölüm için ~6,33 m/s düzlem içi ateşleme yapıldı [D] |
| **S1C 693 → 530 km dairesel** (1B bertaraf hedefi) | **88,06** | hedef [M] CGMS-52 s.25 |
| **S2 786 → 653 km dairesel** (Airbus planı) | **70,22** | plan 73 m/s → %4 sonlu yanma payı [M] |
| S1C 693 km → 400 km yerberi (V2 senaryosu) | 79,85 | tek teğetsel impuls; artık kullanılmıyor [S] |
| S2C 786 km → 400 km yerberi (V2 senaryosu) | 103,98 | aynı [S] |

## 4. Gerçekleşen S2C manevraları (SpaceOps 2025 #617)

| Kalem | m/s |
|---|---:|
| 4 büyük düzlem içi ilerleme manevrası | ~6,33 |
| İnce edinim, 29 Ekim 2024 | 0,012 |
| LEOP sonrası 2 düzlem dışı manevra | 0,57 |
| Tablo 3 (31 Ekim – 18 Aralık 2024), vektör büyüklükleri; 6 Kasım: √(0,9344² + 0,0537²) = 0,936 | 1,056 |
| **Toplam** | **7,97** |

220 s ile 1143 kg'dan **4,21 kg** (200 s: 4,63 kg) [H]. 6 Kasım manevrası 314,1 s sürdü:
ortalama ivme 0,936/314,1 = **2,98 mm/s²**, etkin itki ≈ 1139 kg × 2,98 mm/s² ≈ **3,4 N**
(dört eğik 1 N iticiden). Aynı manevra 9 yıllık 2A'da 422,5 s sürdü (~%27 daha düşük ivme);
blow-down basınç düşüşüyle tutarlıdır [H].

**Rutin bakım:** 2A'da 2 Temmuz 2015 – Mart 2024 arasında 64 manevra × ortalama 1,56 cm/s =
0,998 m/s → **0,115 m/s/yıl** düzlem içi (220 s ve 1100 kg'da yılda ~58 g) [H]. Bu sayım
en büyüğü 3,85 cm/s olan manevraları kapsar; eğim (düzlem dışı) manevraları dahil değildir.
2A'da 24 çift çarpışmadan kaçınma manevrası yapıldı, en büyüğü 19,0 cm/s; toplam için üst
sınır 24 × 2 × 0,19 = 9,1 m/s'dir (gerçek toplam yayımlanmadı) [H].

## 5. CONOPS'a göre ΔV ve yakıt (V3 senaryoları)

### 5.1 Sentinel-1C

Kalemler: LEOP = S1B planlama bütçesinin 1σ değerleri (eğim 3,93 + dışmerkezlik 1,78 + yuva
edinimi 5,3 = **11,01 m/s**) [M]; rutin bakım $d$ m/s/yıl (**1C için yayımlanmadı**, duyarlılık
olarak verilir) [S]; çarpışmadan kaçınma 3 m/s [S]; ömür sonu 693 → 530 km = 88,06 m/s [H].

| d (m/s/yıl) | Süre | Toplam ΔV | Yakıt @220 s | Yakıt @205 s | 154 kg'a göre fark @220 s |
|---:|---:|---:|---:|---:|---:|
| 0,5 | 7,25 yıl | 105,7 m/s | 104,5 kg | 112,0 kg | 49,5 kg |
| 1,0 | 7,25 yıl | 109,3 m/s | 108,0 kg | 115,7 kg | 46,0 kg |
| 1,0 | 12 yıl | 114,1 m/s | 112,6 kg | 120,6 kg | 41,4 kg |
| 2,0 | 12 yıl | 126,1 m/s | 124,1 kg | 132,9 kg | 29,9 kg |

### 5.2 Sentinel-2C

Kalemler: gerçekleşen 7,97 m/s (4,21 kg) [D/H]; düzlem içi bakım 0,115 m/s/yıl [H];
düzlem dışı manevra yılda $n$ adet × 0,936 m/s ($n$ yayımlanmadı) [S]; ömür sonu planı
**47 kg** [M].

| Süre | n (yıl⁻¹) | Rutin ΔV | Rutin yakıt @220 s | 47 kg EOL sonrası kalan |
|---:|---:|---:|---:|---:|
| 7 yıl | 0 | 0,81 m/s | 0,42 kg | 81,4 kg |
| 7 yıl | 1 | 7,36 m/s | 3,88 kg | 77,9 kg |
| 12 yıl | 0,5 | 7,00 m/s | 3,69 kg | 78,1 kg |
| 12 yıl | 1 | 12,61 m/s | 6,64 kg | 75,1 kg |

Kalan tutar güvenli mod tüketimini, çarpışmadan kaçınmayı ve kullanılamayan kalıntıyı
kapsamaz. Sonuç: S2C için itici, 12 yıllık ufukta bile **sınırlayıcı kaynak görünmüyor**;
ömrü batarya, MSI ve diğer bileşenler belirleyecektir. Bu bir yıl sayısı tahmini değildir.

## 6. Kuvvet, kütle debisi, yanma süresi

$\dot m = F/(I_{sp}g_0)$, $t = (m_i-m_f)/\dot m$. 1 N, 220 s: $\dot m$ = 0,46351 g/s.

| | S1C | S2C |
|---|---:|---:|
| 1 m/s için yakıt | 1,0130 kg | 0,5297 kg |
| Tek net 1 N ile süre | 36,42 dk | 19,05 dk |
| Başlangıç ivmesi (1 N) | 0,457 mm/s² | 0,875 mm/s² |
| EOL yakıtı @220 s | 87,43 kg (693 → 530 km) | 47 kg (plan) |
| EOL toplam impulsu | 188 618 N·s | 101 401 N·s |
| EOL için tek 1 N eşdeğeri ateşleme | 52,4 sa | 28,2 sa |

**İtici impuls kapasitesiyle karşılaştırma:** S1C'nin EOL impulsu tek bir Rafael Tip 1
(100 000 N·s) veya Tip 2 (60 000 N·s) iticinin katalog değerini aşar; yük paylaşımı ideal
olsa bile en az 2 Tip 1 ya da 4 Tip 2 itici gerekir [H]. 1B bu nedenle bertarafta birden çok
itici kullandı. S2'de dört itici aynı anda ateşlendiğinden iticide düşen pay (~25 000 N·s)
ArianeGroup kalifikasyonunun (135 000 N·s) altındadır [H].

**Airbus planı tutarlılık kontrolü:** 73 m/s için 47 kg, EOL başlangıç kütlesi 1057–1139 kg
aralığında **etkin Isp ≈ 164–177 s** eder. 220 s nominalin altında; kaynak nedenini açıklamaz
(düşük blow-down basıncında performans kaybı ve/veya marj olası) [H].

## 7. Güç ve tutulma

Dairesel, küresel Dünya, $\beta=0$ için periyot $T=2\pi\sqrt{r^3/\mu}$ ve en uzun tutulma
$T\cdot\arcsin(R_E/r)/\pi$:

| | h (km) | T (dk) | β = 0 tutulma (dk) |
|---|---:|---:|---:|
| S1C | 693 | 98,63 | 35,30 |
| S2C | 786 | 100,58 | 35,15 |

Bunlar **üst sınırdır**. Sentinel-1'in 18:00 LTAN şafak–alacakaranlık yörüngesinde güneş
yörünge düzlemine neredeyse diktir (β büyük); tutulma yılın yalnızca bir döneminde görülür
(SentiWiki SAR görev çevrimini "tutulma mevsimi" için ayrı tanımlar). S2C'nin 10:30 LTDN
yörüngesinde her turda tutulma vardır. S2C kaba güç kontrolü: 1700 W EOL üretim − 1,4 kW
nominal tüketim ≈ 300 W fark [H]; batarya enerjisi için gerilim ve izin verilen deşarj
derinliği [Y].

## 8. Link bütçesi örneği (tamamen varsayımsal)

$L_{fs}=20\log_{10}(4\pi Rf/c)$, $C/N_0=EIRP+G/T-L_{fs}-L_{diğer}+228{,}6$. Varsayım [S]:
f = 8,2 GHz, R = 1800 km, 10 W (10 dBW), 12 dBi, G/T = 25 dB/K, 3 dB diğer kayıp →
$L_{fs}$ = 175,83 dB, $C/N_0$ = 96,77 dBHz; 260 Mbit/s için Eb/N0 = 12,62 dB, 560 Mbit/s
için 9,29 dB. Uçuş donanımının EIRP/G/T değerleri yayımlanmadığından **gerçek marj değildir**.

## 9. Veri birimleri

S1: 1410 Gbit = 176,25 GB; 2 × 260 Mbit/s × 600 s = 39,0 GB. S2: 2,4 Tbit = 300 GB;
560 Mbit/s × 600 s = 42,0 GB. (1 GB = 8 Gbit; ondalık birimler.)

## 10. Rafael itici ısıtıcısı ve valf

Valf 9,2 W @ 28 Vdc; ısıtıcı eleman başına 257 Ω → 28 V'ta $V^2/R$ = 3,05 W/eleman.
Isıtıcı görev çevrimi ve eşzamanlı çalışan eleman sayısı bilinmediğinden sistem gücü
toplanamaz [Y].

## 11. TRL

Uçuş mirası olan bileşen TRL 9 düzeyinde uçmuş tasarımdır; uçuş öncesi seviye ancak belgeli
bir olgunluk değerlendirmesiyle yazılabilir. Hesapla TRL türetilmez. [TRL-ESA]

## 12. "Kaç yıl daha?" sorusu

Görev başına sorulan "hesaplanan miktar ile taşınan miktar arasındaki fark kaç yıl
eder?" sorusunun dürüst yanıtı:

- **S2C:** Gerçekleşen tüketim (4,2 kg), ömür sonu planı (47 kg) ve düzlem içi bakım
  (~58 g/yıl) bilindiğinde ~75–81 kg kalır; düzlem içi bakımla bölmek binlerce yıl verir.
  Bu sonuç anlamsızdır ve asıl sınırlayıcının yakıt **olmadığını** gösterir: eğim
  manevralarının sıklığı, güvenli modlar, kullanılamayan kalıntı ile batarya ve MSI ömrü
  belirleyicidir. Ailenin tasarımı itici ve bataryayı 12 yıl için boyutlandırır. [M]
- **S1C:** Rutin bakım tüketimi yayımlanmadığından yıl hesabı yapılamaz; §5.1 tablosu
  yalnızca senaryodur. Somut veri: S1B'de 2016–2022 işletmesi ve bertaraf manevralarının
  bir kısmından sonra Nisan 2024'te hâlâ ~88 kg yakıt vardı. [M]

## Kaynaklar

[KAYNAKLAR.md](KAYNAKLAR.md); başlıcaları SpaceOps 2025 #617
<https://publications.spaceops.org/2025/download_by_id.php?id=0617>, ISSFD 2017
<https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf>, Rafael katalog
<https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf>,
ArianeGroup 1 N <https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf>.
