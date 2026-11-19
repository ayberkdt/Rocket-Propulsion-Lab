"""Sentinel-1/2 itki sunumunun hesaplarını bağımsız yöntemlerle ve sunum metniyle sınar.

Katmanlar:
  A. Roket denklemi: 60 basamaklı Decimal, ters denklem, sıralı/toplu eşdeğerlik, mevcut seri modeli,
     kütüphanenin SASS ters kütle kapanışı.
  B. Yörünge: periyot, Hohmann (SASS), elips indirme (seri modeli), sonlu yanma cezası (sayısal integral).
  C. Sentinel-2C bütçesi: toplamlar, 2A kaydından öngörü, 6 Kasım etkin kuvveti, sonucun dayanıklılığı.
  D. Kaynak tutarlılığı: makalelerdeki sayıların bağımsız fizikle uyumu.
  E. Donanım sınırları: tank ve itici kataloğuna göre blowdown, artık yakıt, itici ömrü.
  F. Kanıt defteri: her girdi ve tablo hücresi birincil alıntıya bağlı.
  G. Sunum metni: v32'deki sayılar hesapla aynı.
Girdilerin dayanağı: docs/missions/sentinel/sunum/uretim/kanitlar.json.
"""
import importlib.util
import json
import math
import re
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUNUM = ROOT / "docs/missions/sentinel/sunum"
URETIM = SUNUM / "uretim"
DECK = SUNUM / "surumler/Sentinel_1_2_Propulsion_ITU_v32.pptx"


def _yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    modul = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modul)
    return modul


H = _yukle("sunum_hesap", URETIM / "hesap.py")
SERI = _yukle("sentinel_budget_for_sunum", ROOT / "scripts/sentinel_series_budget.py")
KANIT = json.loads((URETIM / "kanitlar.json").read_text(encoding="utf-8"))
D = H.sunum_degerleri()

try:
    from pptx import Presentation
except ImportError:  # python-pptx yoksa yalnız G katmanı atlanır
    Presentation = None


def dec_yakit(m0, dv, isp):
    with localcontext() as c:
        c.prec = 60
        m0, dv, isp = (Decimal(str(x)) for x in (m0, dv, isp))
        return m0 * (1 - (-dv / (isp * Decimal("9.80665"))).exp())


# ====================================================================== A
class A_RoketDenklemi(unittest.TestCase):
    def test_butce_kalemleri_60_basamakli_decimal_ile_ayni(self):
        b = D["butce"]
        with localcontext() as c:
            c.prec = 60
            m = Decimal(1143)
            kalemler = []
            for dv in (b["dv_edinim"], H.isletme_dv_yillik()[0] * 12, H.isletme_dv_yillik()[1] * 12):
                mp = dec_yakit(m, dv, 220)
                kalemler.append(mp)
                m -= mp
            eol = dec_yakit(m, 73, 220)
        self.assertTrue(math.isclose(b["kg_edinim"], float(kalemler[0]), rel_tol=1e-13))
        self.assertTrue(math.isclose(b["kg_isletme"], float(kalemler[1] + kalemler[2]), rel_tol=1e-13))
        self.assertTrue(math.isclose(b["kg_eol_ideal"], float(eol), rel_tol=1e-13))

    def test_ters_denklem_ileri_denklemi_geri_verir(self):
        for m0, dv, isp in [(1143, 7.398, 220), (2157, 6.3, 214), (1136.97, 73, 220), (1000, 1e-9, 200)]:
            with self.subTest(m0=m0, dv=dv):
                self.assertTrue(math.isclose(H.kapasite(m0, H.yakit(m0, dv, isp), isp), dv, rel_tol=1e-12))

    def test_sirali_defter_tek_ustele_esittir(self):
        dvler = [D["butce"]["dv_edinim"], 1.375, 2.645, 73.0]
        satirlar = H.sirali(1143, dvler, 220)
        self.assertTrue(math.isclose(sum(s[1] for s in satirlar), H.yakit(1143, sum(dvler), 220), rel_tol=1e-13))
        self.assertTrue(math.isclose(satirlar[-1][2] + sum(s[1] for s in satirlar), 1143, rel_tol=1e-15))

    def test_kucuk_dv_iptal_etmez(self):
        self.assertTrue(math.isclose(H.yakit(1000, 1e-12, 220), float(dec_yakit(1000, 1e-12, 220)), rel_tol=1e-12))

    def test_mevcut_seri_modeliyle_ayni(self):
        for m0, dv, isp in [(1143, 7.4, 220), (2186, 100, 205), (2157, 4.1, 214)]:
            self.assertEqual(H.yakit(m0, dv, isp), SERI.propellant(m0, dv, isp))
        satirlar = SERI.ledger(1143, 133, json.loads(
            (ROOT / "docs/missions/sentinel/SERI_VERILERI.json").read_text(encoding="utf-8"))["s2c_events"], 220)
        model_dv = sum(r["dv_m_s"] for r in satirlar)
        # Seri modeli planlanan 0,57 m/s'yi de içerir; sunum bütçesi içermez.
        self.assertAlmostEqual(model_dv - D["butce"]["dv_edinim"], 0.57, places=9)

    def test_kutuphanenin_ters_kutle_kapanisi_1143_kg_verir(self):
        from rocket_propulsion.studies.sass import Maneuver, assumption, size_mission
        q = lambda x, u="1": assumption(x, u, "Sunum testi: hesap.py değeri")
        b = D["butce"]
        kuru = 1143 - (b["kg_edinim"] + b["kg_isletme"] + b["kg_eol_ideal"])
        burns = tuple(Maneuver(n, q(dv, "m/s"), q(1, "N"), q(220, "s"), q(1), q(1)) for n, dv in
                      (("edinim", b["dv_edinim"]), ("bakim", H.isletme_dv_yillik()[0] * 12),
                       ("kacinma", H.isletme_dv_yillik()[1] * 12), ("bertaraf", 73.0)))
        sonuc = size_mission(q(kuru, "kg"), burns, residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        islak = sonuc.summary.to_dict()["values"]["wet_mass"]["value"]
        self.assertTrue(math.isclose(islak, 1143, rel_tol=1e-9))

    def test_gecersiz_girdiler_reddedilir(self):
        for args in [(1000, -1, 220), (1000, float("nan"), 220), (1000, 1, 0), (True, 1, 220), (1000, 1, 220, 1.1)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                H.yakit(*args)
        with self.assertRaises(ValueError):
            H.kapasite(1000, 1000, 220)


# ====================================================================== B
class B_Yorunge(unittest.TestCase):
    def test_periyot_seri_modeliyle_ayni(self):
        for h in (693e3, 786e3):
            self.assertTrue(math.isclose(H.periyot(h), SERI.orbit(h)["period_s"], rel_tol=1e-15))
        self.assertAlmostEqual(D["yanma"]["periyot_dk"], 100.5786, places=3)

    def test_hohmann_kutuphanenin_sass_fonksiyonuyla_ayni(self):
        from rocket_propulsion.studies.sass import assumption, hohmann_transfer
        q = lambda x, u: assumption(x, u, "Sunum testi")
        sonuc = hohmann_transfer(q(H.R_E + 786e3, "m"), q(H.R_E + 653e3, "m"), q(H.MU, "m^3/s^2"))
        self.assertTrue(math.isclose(D["bertaraf"]["hohmann"], sonuc.get("total_delta_v").value, rel_tol=1e-12))
        self.assertTrue(math.isclose(D["bertaraf"]["hohmann"], SERI.lowering_dv(786e3, 653e3, 653e3), rel_tol=1e-12))

    def test_1b_elips_indirme_seri_modeliyle_ayni(self):
        self.assertTrue(math.isclose(D["bertaraf"]["s1b_ideal"], SERI.lowering_dv(693e3, 580e3, 610e3), rel_tol=1e-12))
        self.assertTrue(52.5 < D["bertaraf"]["s1b_ideal"] < 52.7)

    def test_sonlu_yanma_cezasi_sayisal_integralle_ayni(self):
        # Düğüme ortalanmış sabit itkide etkin eğim değişimi ∝ ortalama cos(u), u ∈ [−θ/2, θ/2].
        theta = math.radians(D["yanma"]["yay_derece"])
        n = 20000
        h = theta / n
        simpson = sum((1 if i in (0, n) else 4 if i % 2 else 2) * math.cos(-theta / 2 + i * h) for i in range(n + 1))
        ortalama = simpson * h / 3 / theta
        self.assertAlmostEqual(D["yanma"]["ceza"], 1 - ortalama, places=12)
        self.assertAlmostEqual(D["yanma"]["ceza"], theta ** 2 / 24, delta=theta ** 4 / 1500)
        self.assertTrue(0.004 < D["yanma"]["ceza"] < 0.005)

    def test_sonlu_yanma_cezasi_sinirlari(self):
        self.assertEqual(H.sonlu_yanma_cezasi(0), 0.0)
        self.assertAlmostEqual(H.sonlu_yanma_cezasi(2 * math.pi), 1.0, places=12)
        self.assertLess(H.sonlu_yanma_cezasi(0.1), H.sonlu_yanma_cezasi(0.2))


# ====================================================================== C
class C_Sentinel2CButcesi(unittest.TestCase):
    def test_toplamlar_tutarli(self):
        b = D["butce"]
        self.assertAlmostEqual(b["dv_toplam"], b["dv_edinim"] + b["dv_isletme"] + 73, places=12)
        self.assertAlmostEqual(b["kg_toplam"], b["kg_edinim"] + b["kg_isletme"] + 47, places=12)
        self.assertAlmostEqual(b["ayrilmamis_kg"], 133 - b["kg_toplam"], places=12)
        self.assertAlmostEqual(b["bertaraf_payi"], 73 / b["dv_toplam"], places=15)

    def test_edinim_dv_bilesenleri(self):
        t3 = 0.0146 + math.hypot(0.9344, 0.0537) + 0.0025 + 0.0207 + 0.0206 + 0.0309 + 0.0029 + 0.0283
        self.assertAlmostEqual(t3, 1.0564418, places=7)
        self.assertAlmostEqual(D["butce"]["dv_edinim"], 6.33 + 0.012 + t3, places=12)
        # Ters yönlü kaçınma çifti yakıtta birbirini götürmez: büyüklükler toplanır.
        self.assertGreater(0.0207 + 0.0206, abs(-0.0207 + 0.0206) * 100)

    def test_isletme_ongorusu_2a_kaydindan(self):
        yil = 2024.21 - 2015.50
        bakim, kacinma = H.isletme_dv_yillik()
        self.assertAlmostEqual(bakim, 64 * 0.0156 / yil, places=14)
        self.assertAlmostEqual(kacinma, 24 * 2 * 0.04 / yil, places=14)
        self.assertAlmostEqual(D["butce"]["dv_isletme"], 4.02, delta=0.005)
        self.assertAlmostEqual(D["butce"]["isletme_yillik"], 0.335, delta=0.001)

    def test_6_kasim_etkin_kuvveti_birden_cok_itici_gerektirir(self):
        f = D["yanma"]["etkin_kuvvet_N"]
        self.assertTrue(3.38 < f < 3.41)
        self.assertGreater(f, 1.1, "tek bir itici kataloğun 1,1 N üst sınırını aşamaz")
        satirlar = SERI.ledger(1143, 133, json.loads(
            (ROOT / "docs/missions/sentinel/SERI_VERILERI.json").read_text(encoding="utf-8"))["s2c_events"], 220)
        self.assertAlmostEqual(f, satirlar[4]["inferred_average_axial_thrust_N"], delta=0.005)

    def test_bertaraf_en_buyuk_kalem_kalir(self):
        b = D["butce"]
        # İşletme Δv'si üç katına çıksa ve Isp 200 s olsa bile bertaraf toplamın yarısından fazladır.
        diger = b["dv_edinim"] + 3 * b["dv_isletme"]
        self.assertGreater(73 / (73 + diger), 0.7)
        # Sıralamayı değiştirmek için yayımlanmamış eğim kontrolü 12 yılda en az bu kadar Δv istemeli:
        esik = 73 - (b["dv_edinim"] + b["dv_isletme"])
        self.assertGreater(esik, 60)
        self.assertGreater(esik / D["yanma"]["dv_6kasim"], 60, "6 Kasım büyüklüğünde 60'tan fazla eğim manevrası")

    def test_isp_duyarliligi(self):
        alt, ust = H.yakit(1143, D["butce"]["dv_edinim"], 223), H.yakit(1143, D["butce"]["dv_edinim"], 200)
        self.assertTrue(3.85 < alt < 3.87 and 4.30 < ust < 4.32)
        self.assertEqual(f"{D['butce']['kg_edinim']:.1f}", "3.9")


# ====================================================================== D
class D_KaynakTutarliligi(unittest.TestCase):
    k = D["kaynak"]

    def test_1a_enjeksiyon_duzeltmesi_gauss_ile_uyumlu(self):
        self.assertLess(abs(4.1 / self.k["dv_sma_min"] - 1), 0.02)
        self.assertGreaterEqual(2.2, self.k["dv_inc_min"])
        self.assertLess(2.2 / self.k["dv_inc_min"] - 1, 0.06)
        self.assertGreaterEqual(6.3, self.k["dv_sma_min"] + self.k["dv_inc_min"])

    def test_1a_yer_izi_kaymasi_gunde_65_km(self):
        self.assertLess(abs(self.k["kayma_km_gun"] / 65 - 1), 0.03)

    def test_tork_yaklasik_on_kat(self):
        self.assertTrue(8 < self.k["tork_orani"] < 11)

    def test_30_ve_300_saniyelik_yanma_buyuklukleri(self):
        self.assertTrue(0.009 < self.k["dv_30s"] < 0.011, "≈ 1 cm/s")
        self.assertTrue(0.135 < self.k["dv_300s"] < 0.145, "neredeyse 14 cm/s")

    def test_2c_dort_buyuk_manevra_alt_sinirin_ustunde(self):
        self.assertGreaterEqual(6.33, self.k["dv_2c_min"])

    def test_1a_ucus_tuketimi_ideal_ve_basit_tahminin_ustunde(self):
        s = D["s1a"]
        self.assertTrue(s["ideal_kg"] < s["kayipli_alt_kg"] < s["kayipli_ust_kg"] < s["ucusta_kg"])
        self.assertEqual((f"{s['ideal_kg']:.1f}", f"{s['kayipli_alt_kg']:.1f}", f"{s['kayipli_ust_kg']:.1f}"),
                         ("6.5", "7.6", "8.3"))

    def test_bertaraf_plani_ideal_ustunde_ve_etkin_isp_katalogun_altinda(self):
        b = D["bertaraf"]
        self.assertAlmostEqual(b["plan_fazlasi"], 0.0397, delta=0.0005)
        self.assertGreater(47, b["kg_eol_ideal"])
        lo, hi = b["isp_etkin"]
        self.assertTrue(163 < lo < 165 and 176 < hi < 178)
        self.assertLess(hi, 200, "47 kg, katalog Isp alt sınırından düşük etkin Isp'ye karşılık gelir")


# ====================================================================== E
class E_DonanimSinirlari(unittest.TestCase):
    def test_blowdown_orani_ve_izotermal_son_basinc(self):
        t = D["tank"]
        self.assertAlmostEqual(t["gaz_bas_L"], 50.61, delta=0.01)
        self.assertAlmostEqual(204 / t["gaz_bas_L"], 4.0, delta=0.05, msg="4:1 blowdown oranı")
        p_son = 27.4 * t["gaz_bas_L"] / 204
        self.assertGreater(p_son, 5.5, "izotermal genişlemede son basınç katalog alt sınırının üstünde")

    def test_artik_yakit_ust_siniri(self):
        self.assertAlmostEqual(D["tank"]["artik_ust_kg"], 1.54, places=9)

    def test_ideal_kapasiteler(self):
        t = D["tank"]
        self.assertEqual([round(x) for x in t["kap_1c"]], [147, 158])
        self.assertEqual([round(x) for x in t["kap_2c"]], [243, 271])
        self.assertGreater(t["kap_2c"][0] / D["butce"]["dv_toplam"], 2.8)

    def test_kutle_duyarliligi(self):
        t = D["tank"]
        self.assertAlmostEqual(t["kg_per_m_s_1c"], 1.013, delta=0.001)
        self.assertAlmostEqual(t["kg_per_m_s_2c"], 0.530, delta=0.001)
        self.assertAlmostEqual(D["yanma"]["dk_per_m_s_1N"], 19.05, delta=0.01)

    def test_arianegroup_itici_omru_bertarafa_yeter(self):
        o = D["omur"]
        self.assertAlmostEqual(o["eol_yanma_h"], 19.333, places=3)
        self.assertLess(o["eol_yanma_payi"], 1)
        self.assertLess(o["tek_yanma_payi"], 1)
        self.assertLess(o["yakit_gecis_payi"], 1)
        self.assertLess(o["impuls_2c_payi"], 1)

    def test_rafael_tip2_impulsu_en_az_alti_itici_ister(self):
        kat = D["omur"]["impuls_1c_tip2_kati"]
        self.assertTrue(5 < kat < 6)
        self.assertLess(kat, 7, "dalda yedi itici var")

    def test_katalog_nominalleri_tek_calisma_noktasi_degil(self):
        self.assertAlmostEqual(0.44e-3 * 220 * H.G0, 0.949, delta=0.001)
        self.assertAlmostEqual(D["yanma"]["debi_g_s"], 0.4635, delta=0.0001)


# ====================================================================== F
KAYIT = {k["id"]: k for k in KANIT["kanitlar"]}


class F_KanitDefteri(unittest.TestCase):
    kayit = KAYIT

    def test_her_girdi_defterde(self):
        for anahtar, (_, kaynak) in H.KAYNAK.items():
            with self.subTest(girdi=anahtar):
                self.assertIn(kaynak, self.kayit)

    def test_birincil_kayitlarin_alintisi_ve_belgesi_var(self):
        for k in KANIT["kanitlar"]:
            with self.subTest(kayit=k["id"]):
                self.assertIn(k["belge"], KANIT["belgeler"])
                self.assertIn(k["durum"], {"birincil", "depo", "varsayim"})
                if k["durum"] == "birincil":
                    self.assertTrue(k["alintilar"])
                    self.assertTrue(KANIT["belgeler"][k["belge"]]["url"].startswith("https://"))

    def test_girdi_degeri_alintida_gecer(self):
        beklenen = {  # hesap girdisi → birincil alıntıda geçmesi gereken metin
            "tank_hacim_L": "204 liters", "tank_kapasite": "154 kg", "tank_p_bas": "27.4", "tank_p_son": "5.5 bar",
            "tank_bosaltma": "> 99%", "rafael_isp": ">214 @ 22 bar >205", "rafael_toplam_impuls": "100,000 60,000",
            "ag_isp": "200 … 223", "ag_isp_nom": "220 s", "ag_debi_nom": "0.44 g/s", "ag_toplam_impuls": "135,000",
            "ag_cevrim": "59,000", "ag_yakit_gecis": "67 kg", "ag_tek_yanma_h": "12 h", "ag_toplam_yanma_h": "50 h",
            "yakit_1c": "154 kg", "yakit_2c": "133 kg", "dv_1a_sma": "4.1 m/s", "dv_1a_inc": "2.2 m/s",
            "da_1a": "-7.8 km", "di_1a_mdeg": "-16 mdeg", "kayma_1a_km_gun": "65 k", "yakit_1a_edinim": "8.4 kg",
            "kayip_1a_lo": "-21 %", "kayip_1a_hi": "30%", "yanma_once_s": "300", "yanma_sonra_s": "30 seconds",
            "dv_yanma_once": "14 cm/s", "dv_yanma_sonra": "1 cm/s", "tork_olculen": "0.22, 0.15",
            "tork_ongorulen": "0.0283", "dv_2c_dort_buyuk": "6.33", "dv_2c_29ekim": "1.2 cm/s",
            "sure_6kasim_s": "314.1", "da_2c": "10.6 km", "occ_sayi": "64 orbit control", "occ_ort": "1.56 cm/s",
            "cam_cift": "24 pairs", "cam_ornek": "4 cm/s", "gozlem_bas": "02 July 2015", "gozlem_son": "March 2024",
            "h_s2": "786 km", "h_1b_apogee": "610 km", "ufuk_yil": "12 years",
        }
        for anahtar, metin in beklenen.items():
            kaynak = H.KAYNAK[anahtar][1]
            k = self.kayit[kaynak]
            with self.subTest(girdi=anahtar, kaynak=kaynak):
                self.assertEqual(k["durum"], "birincil")
                self.assertIn(metin, " ".join(k["alintilar"]))

    def test_tablo_hucreleri_kanita_bagli(self):
        for tablo, hucreler in KANIT["tablolar"].items():
            for hucre, kimlikler in hucreler.items():
                with self.subTest(tablo=tablo, hucre=hucre):
                    self.assertTrue(kimlikler)
                    for i in kimlikler:
                        self.assertIn(i, self.kayit)

    def test_depo_kaynakli_kayitlar_acikca_isaretli(self):
        depo = sorted(k["id"] for k in KANIT["kanitlar"] if k["durum"] == "depo")
        self.assertIn("S2-AOCS", depo)
        self.assertTrue(all(not self.kayit[i]["alintilar"] for i in depo))

    def test_alinti_dogrulama_raporu_tamam(self):
        rapor = (URETIM / "kanit_dogrulama_2026-11-19.txt").read_text(encoding="utf-8")
        birincil_alinti = sum(len(k["alintilar"]) for k in KANIT["kanitlar"] if k["durum"] == "birincil")
        self.assertIn("sonuç: TAMAM", rapor)
        self.assertEqual(rapor.count("BULUNDU"), birincil_alinti)
        self.assertNotIn(" YOK ", rapor)


# ====================================================================== G
@unittest.skipIf(Presentation is None, "python-pptx kurulu değil")
class G_SunumMetni(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.slayt = []
        for s in Presentation(DECK).slides:
            parca = []
            for sh in s.shapes:
                if sh.has_text_frame:
                    parca.append(sh.text_frame.text)
                if getattr(sh, "has_table", False) and sh.has_table:
                    parca += [" | ".join(c.text for c in r.cells) for r in sh.table.rows]
            cls.slayt.append(re.sub(r"\s+", " ", " / ".join(parca)))

    def icinde(self, no, *metinler):
        for m in metinler:
            with self.subTest(slayt=no, metin=m):
                self.assertIn(m, self.slayt[no - 1])

    def test_slayt8_butce_tablosu(self):
        b = D["butce"]
        self.icinde(8, f"Acquisition, 2024 | {b['dv_edinim']:.1f} | {b['kg_edinim']:.1f} | Flown",
                    f"Operations, 12 years | {b['dv_isletme']:.1f} | {b['kg_isletme']:.1f} |",
                    f"| 73.0 | 47.0 | A/B design value ({b['kg_eol_ideal']:.1f} kg ideal)",
                    f"Total | ≈ {b['dv_toplam']:.1f} | ≈ {b['kg_toplam']:.1f} |",
                    f"Disposal is {b['bertaraf_payi'] * 100:.0f} % of the illustrative",
                    "m₀ = 1,143 kg, Isp = 220 s")

    def test_slayt6_dv_satiri(self):
        b = D["butce"]
        self.icinde(6, f"{b['dv_edinim']:.1f} m/s", f"{b['dv_isletme']:.1f} m/s", "73 m/s", "3 days", "~500 burns")

    def test_slayt7_s1a_ve_s1b(self):
        s = D["s1a"]
        self.icinde(7, f"8.4 kg ({s['ideal_kg']:.1f} kg ideal)", "21–30 % low", "limited to 30 s, ≈ 1 cm/s",
                    "m₀ 2,157 kg, Isp 214 s", "RCT-3B", "~580 × 610 km", "electrically passivated only")

    def test_slayt9_yakit_ve_yanma(self):
        y = D["yanma"]
        self.icinde(9, f"{y['debi_g_s']:.2f} g/s", f"{y['oran'] * 100:.1f} % of an orbit",
                    f"≈ {y['ceza'] * 100:.2f} %", "1,143 kg", "220 s")

    def test_slayt10_sonuclar(self):
        b = D["butce"]
        self.icinde(10, f"{b['dv_edinim']:.1f} m/s for 2C's acquisition", f"~{b['isletme_yillik']:.1f} m/s per year",
                    f"{b['bertaraf_payi'] * 100:.0f} % of the illustrative 12-year")

    def test_slayt3_mimari_tablosu_katalogla_ayni(self):
        self.icinde(3, "204 L", "27.4 → 5.5 bar", "154 kg (1C)", "133 kg (2C)", "21 bar at fuelling",
                    "2 × 7 thrusters of 1 N (3 orbit, 4 attitude)", "2 × 4 canted thrusters of 1 N",
                    "Module inside the launcher interface ring", "orbit thrusters on ±X, −Y")

    def test_slayt5_itici_tablosu_katalogla_ayni(self):
        self.icinde(5, "5.5–24.5 bar", "5.5–22 bar", "1.0–1.4 N", "0.32–1.1 N", "205–214 s", "200–223 s",
                    "0.008 N·s", "0.010–0.043 N·s", "Sentinel-1 since 2014", "Sentinel-2 since 2015")

    def test_slayt1_filo_durumu(self):
        self.icinde(1, "2A extended to end-2027", "693 km", "786 km")

    def test_sunumdaki_her_tablo_hucresinin_kaniti_var(self):
        sutunlar = {"Table 2": ("Sentinel-1", "Sentinel-2"), "Table 3": ("Rafael", "ArianeGroup")}
        for tablo, slayt_no in (("Table 2", 3), ("Table 3", 5)):
            tbl = next(sh for sh in Presentation(DECK).slides[slayt_no - 1].shapes
                       if getattr(sh, "has_table", False) and sh.has_table).table
            for r in list(tbl.rows)[1:]:
                satir = r.cells[0].text
                for sutun in sutunlar[tablo]:
                    with self.subTest(tablo=tablo, hucre=f"{satir} / {sutun}"):
                        self.assertIn(f"{satir} / {sutun}", KANIT["tablolar"][tablo])

    def test_tablo_hucre_anahtarlari_sunumdaki_satirlarla_eslesir(self):
        for tablo, slayt_no in (("Table 2", 3), ("Table 3", 5)):
            for hucre in KANIT["tablolar"][tablo]:
                satir = hucre.split(" / ")[0]
                with self.subTest(tablo=tablo, satir=satir):
                    self.assertIn(satir + " |", self.slayt[slayt_no - 1])


if __name__ == "__main__":
    unittest.main()
