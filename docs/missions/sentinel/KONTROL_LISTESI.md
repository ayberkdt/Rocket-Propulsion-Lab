# Sentinel-1 / Sentinel-2 — 27 maddelik ödev kontrolü (V4)

**29 Ekim 2026.** Her satır iki seri için ayrı değerlendirilir; A/B/C/D ayrı araçlardır.
**Kaynaklı:** raporda birincil kanıt var. **Kısmi:** konu işlendi ama kamuya açık tam
uçuş bütçesi/araç kaydı yok. **Hesap/V:** türetim ve varsayım; uçuş kabulü değil.
"İşlendi" ile "bütün mühendislik bütçesi kapandı" eş anlamlı değildir.

Ana raporlar: [Sentinel-1 serisi](SENTINEL-1_SERISI.md) ve
[Sentinel-2 serisi](SENTINEL-2_SERISI.md). Sayısal ek: [HESAPLAR.md](HESAPLAR.md).

| Kimlik | Hocanın maddesi | Sentinel-1 kanıt / durum | Sentinel-2 kanıt / durum |
|---|---|---|---|
| A1 | Tanım, amaç, operatör, üretici | A1; kaynaklı ESA/kit | A1; kaynaklı ESA |
| A2 | Görev fotoğrafı | A2; gerçek C fotoğrafı, kredi/lisans | A2; gerçek C dolum fotoğrafı, kredi/lisans |
| A3 | Kütle, boyut, güç, yörünge, ömür | A3; 4 araç kitleri ayrı, yakıt kaydı kısmi | A3; eski tasarım/C dolumu ayrı; B/D kabul kütlesi açık |
| A4 | Fırlatıcı, tarih, rampa, enjeksiyon | A4; hedef/gerçek ayrımı, kesin enjeksiyon vektörü açık | A4; A/B/C Vega; D fırlatılmadı; gerçek vektörler kısmi |
| A5 | Dated operasyonel durum | A3/A5; A/B sonu, C/D ve dinamik 1D olayı | A3/A5; B/C, A uzatması, D depolama |
| A6 | Fırlatıcı/üst kademe dahil sistem | A4/C; Soyuz-Fregat, Vega-C-AVUM+, Ariane6-Vinci | A4/C; Vega-AVUM, adaptör ve uydu sınırı |
| A7 | Genel tartışma | A5/B; SAR ve hassas yer izi/itki bağı | A5/C; optik yük ve tandem/geometri bağı |
| B1 | Kütle bütçesi | B1; alt kümeler, toplam kabul bütçesi açık | B1; hidrazin dışı fark, toplam kabul bütçesi açık |
| B2 | Güç bütçesi | B2; RF/DC/ortalama ve valf-ısıtıcı ayrımı; kısmi | B2; EOL/BOL ve tur enerjisi; kısmi |
| B3 | Yönelim bütçesi | B3; yörünge bandı pointing değil; kovaryans açık | B3; gereksinim/güven düzeyi ve varsayımlı RSS; kısmi |
| B4 | Link bütçesi | B4; kanal ayrımı ve denklem, gerçek EIRP/G/T açık | B4; RF/EDRS ve denklem, gerçek link kapanışı açık |
| B5 | ΔV bütçesi | B7/C; evreler ve B gerçek bertarafı, C/D tahsisi açık | B7/C; C yaklaşık rekonstrüksiyon, A istatistiği; yıllık toplam açık |
| B6 | Veri bütçesi | B5; kapasite, indirme, taşma denklemi; kısmi | B5; bit/byte, net hız, kapasite/indirme; kısmi |
| B7 | Isıl bütçe | B6; denge, katalog sıcaklıkları, uçuş limitleri açık | B6; denge ve ideal radyatör; gerçek model açık |
| B8 | Payload / bütün alt sistem detayları | B8/A3; SAR/AIS/EPS/AOCS/veri/termal arayüzler | B8/A3; MSI/EPS/AOCS/OCP/modül arayüzleri |
| C1 | Fırlatma→EOL CONOPS | C evre tablosu ve B bertaraf deneyimi | C evre tablosu, tandem ve görev devri |
| C2 | CONOPS'tan ΔV çıkarımı | C; açık varsayımlı C senaryosu ve hesap | C1; norm/olay bazlı C defteri; C2/C3 A bakım/drift |
| D1 | İtki tipi | D1; monoprop hidrazin/blow-down | D1; monoprop hidrazin/He blow-down |
| D2 | Kütle, güç, boyut, yerleşim | D1/D2/B1/B2; tank katalog + işlevsel şema; toplam modül açık | D1/D2/B1/B2; halka içi aile modülü + şema; hacim/güç açık |
| D3 | İtkinin işlevleri | D3/C; edinim, bakım, CAM, güvenli mod, EOL | D3/C; faz/tandem, bakım, CAM, güvenli mod, EOL |
| D4 | Tsiolkovsky, ΔV/kütle/kuvvet | D4/hesap §2–5; kuvvet süreye bağlı | D4/hesap §2–5; net/skaler kuvvet ayrımı |
| D5 | Süre ve ṁ | D4–6/hesap §5; ideal eşleşmiş nokta | D4–6; gerçek süre, çıkarılmış kuvvet, ideal ṁ |
| D6 | Propellant hesabı | D4–6/C; ideal kapasite ve senaryo | D4–6; kademeli defter ve 4 duyarlılık durumu |
| D7 | Taşınan/hesap rezervi ve ek ömür | D7; B 70 kg karşı örneği; C kalan ölçümü açık | D7; 47 kg eski plan/ideal tüketim ayrımı; ek yıl koşullu |
| D8 | Impulsive / continuous | D8; kimyasal ayrık; 700 s sonlu yanma oranı | D8; saniyelik/314/600 s ve yörünge periyodu |
| D9 | Güvenlik/regülasyon | D9; tehlike, FDIR, 25/5 yıl baseline; uyum kısmi | D9; toksisite/basınç/termal, 25/5 yıl baseline; uyum kısmi |
| D10 | Uçuş öncesi/sonrası TRL | D10; miras/tekil konfigürasyon ayrımı; resmi eski puan açık | D10; A/B uçuşu, C kayıt sınırı, D FAR uçuş değil |

Sonraki kendi motor+besleme tasarımımızda gereksinim ve temsilî ortam deneyi eksikleri
ayrı kapatılmalıdır. Bu tabloda hiçbir "kaynaklı" satır, yeni sistemin TRL 5 olduğu
anlamına gelmez.
