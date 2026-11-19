"""Sentinel-1/2 itki sunumunun bütün sayıları: tek kaynak, formül ve dayanakla.

Her girdi `KAYNAK` sözlüğündeki bir kanıt kimliğine bağlıdır (bkz. kanitlar.json). Fonksiyonlar saftır,
dışarıya dosya yazmaz ve standart kütüphane dışında bağımlılık kullanmaz. `sunum_degerleri()` sunumda
görünen bütün değerleri tek sözlükte döndürür; tests/test_sentinel_sunum.py bu değerleri bağımsız
yöntemlerle (Decimal, sayısal integral, kütüphanenin SASS fonksiyonları, mevcut seri modeli) ve
v32 sunumunun metniyle karşılaştırır.

Varsayımlar [H]: ideal (sabit Isp, η = 1) roket denklemi; sıralı kütle güncellemesi; iki cisim yörüngesi.
"""
from __future__ import annotations

import math

# ---------------------------------------------------------------- sabitler
G0 = 9.80665                      # m/s², standart yerçekimi (BIPM)
MU = 3.986004418e14               # m³/s², Yer'in kütleçekim parametresi (WGS-84)
R_E = 6378137.0                   # m, Yer ekvator yarıçapı (WGS-84)
OMEGA_E = 7.2921159e-5            # rad/s, Yer'in dönme hızı

# ---------------------------------------------------------------- yayımlanmış girdiler (kanıt kimlikleriyle)
KAYNAK = {
    # kütle ve yakıt
    "m0_1c": (2186.0, "VV25-KIT"), "m0_2c": (1143.0, "VV24-KIT"), "m0_1a": (2157.0, "VS07-KIT"),
    "yakit_1c": (154.0, "S1-FUEL"), "yakit_2c": (133.0, "S2-FUEL"),
    # yörüngeler
    "h_s1": (693e3, "S1-ESA"), "h_s2": (786e3, "S2-SP1322"), "h_s2_eol": (653e3, "S2-AOCS"),
    "h_1b_perigee": (580e3, "SO226-P20"), "h_1b_apogee": (610e3, "SO226-P20"),
    # Sentinel-1A edinimi (ISSFD 2015 #134 s.8 Tablo 2, s.14; #137 s.11-14)
    "dv_1a_sma": (4.1, "ISSFD134-T2"), "dv_1a_inc": (2.2, "ISSFD134-T2"),
    "da_1a": (-7.8e3, "ISSFD134-T2"), "di_1a_mdeg": (-16.0, "ISSFD134-T2"),
    "kayma_1a_km_gun": (65.0, "ISSFD134-P8"), "yakit_1a_edinim": (8.4, "ISSFD134-P14"),
    "kayip_1a_lo": (0.21, "ISSFD134-P8"), "kayip_1a_hi": (0.30, "ISSFD134-P8"),
    "yanma_once_s": (300.0, "ISSFD137-P14"), "yanma_sonra_s": (30.0, "ISSFD137-P14"),
    "dv_yanma_once": (0.14, "ISSFD134-P8"), "dv_yanma_sonra": (0.01, "ISSFD134-P8"),
    "tork_olculen": ((0.00, -0.22, 0.15), "ISSFD137-P13"),
    "tork_ongorulen": ((0.0003, 0.0283, -0.0054), "ISSFD137-P13"),
    # Sentinel-2C edinim ve tandem (SpaceOps 2025 #617 s.15-16, Tablo 3)
    "dv_2c_dort_buyuk": (6.33, "SO617-P15"), "dv_2c_29ekim": (0.012, "SO617-P15"),
    "dv_2c_tablo3": ((0.0146, (0.9344, 0.0537), 0.0025, 0.0207, 0.0206, 0.0309, 0.0029, 0.0283), "SO617-T3"),
    "sure_6kasim_s": (314.1, "SO617-T3"), "da_2c": (-10.6e3, "SO617-P15"),
    # Sentinel-2A işletme istatistiği (SpaceOps 2025 #617 s.5)
    "occ_sayi": (64, "SO617-P5"), "occ_ort": (0.0156, "SO617-P5"),
    "cam_cift": (24, "SO617-P5"), "cam_ornek": (0.04, "SO617-P5"),
    "gozlem_bas": (2015.50, "SO617-P5"), "gozlem_son": (2024.21, "SO617-P5"),
    # bertaraf tasarımı (Airbus GNC 2014 s.24, depo raporu üzerinden)
    "dv_eol": (73.0, "S2-AOCS"), "yakit_eol": (47.0, "S2-AOCS"), "eol_yanma_sayi": (116, "S2-AOCS"),
    "eol_yanma_s": (600.0, "S2-AOCS"),
    # tank ve itici katalogları
    "tank_hacim_L": (204.0, "RAFAEL-P44"), "tank_kapasite": (154.0, "RAFAEL-P44"),
    "tank_p_bas": (27.4, "RAFAEL-P44"), "tank_p_son": (5.5, "RAFAEL-P44"), "tank_bosaltma": (0.99, "RAFAEL-P44"),
    "rho_n2h4": (1004.0, "DERIN-P330"),
    "rafael_isp": ((205.0, 214.0), "RAFAEL-P49"), "rafael_toplam_impuls": ((100000.0, 60000.0), "RAFAEL-P49"),
    "ag_isp": ((200.0, 223.0), "AG-P6"), "ag_isp_nom": (220.0, "AG-P6"), "ag_debi_nom": (0.44e-3, "AG-P6"),
    "ag_toplam_impuls": (135000.0, "AG-P6"), "ag_cevrim": (59000, "AG-P6"), "ag_yakit_gecis": (67.0, "AG-P6"),
    "ag_tek_yanma_h": (12.0, "AG-P6"), "ag_toplam_yanma_h": (50.0, "AG-P6"),
    # sunum varsayımları [H]
    "isp_2c": (220.0, "VARSAYIM"), "isp_1a": (214.0, "VARSAYIM"), "ufuk_yil": (12.0, "SP1322-P23"),
}


def v(key):
    """Kayıtlı girdinin değeri."""
    return KAYNAK[key][0]


# ---------------------------------------------------------------- temel fizik
def yakit(m0: float, dv: float, isp: float, verim: float = 1.0) -> float:
    """İdeal roket denklemi: m_p = m0 (1 − exp(−Δv / (η Isp g0))). expm1 küçük Δv'de iptali önler."""
    for ad, x in (("m0", m0), ("isp", isp), ("verim", verim)):
        if isinstance(x, bool) or not math.isfinite(x) or x <= 0:
            raise ValueError(f"{ad} pozitif ve sonlu olmalı")
    if isinstance(dv, bool) or not math.isfinite(dv) or dv < 0:
        raise ValueError("dv negatif olmayan sonlu bir sayı olmalı")
    if verim > 1:
        raise ValueError("verim 1'den büyük olamaz")
    return m0 * -math.expm1(-dv / (verim * isp * G0))


def kapasite(m0: float, mp: float, isp: float) -> float:
    """Ters roket denklemi: Δv = Isp g0 ln(m0 / (m0 − m_p))."""
    if not 0 <= mp < m0:
        raise ValueError("yakıt kütlesi [0, m0) aralığında olmalı")
    return -isp * G0 * math.log1p(-mp / m0)


def sirali(m0: float, dvler, isp: float):
    """Sıralı defter: her kalemden sonra kütle güncellenir. [(dv, mp, m_sonra), ...] döndürür."""
    m, satirlar = m0, []
    for dv in dvler:
        mp = yakit(m, dv, isp)
        m -= mp
        satirlar.append((dv, mp, m))
    return satirlar


def periyot(h: float) -> float:
    a = R_E + h
    return 2 * math.pi * math.sqrt(a ** 3 / MU)


def hiz(h: float) -> float:
    return math.sqrt(MU / (R_E + h))


def hohmann(h1: float, h2: float) -> float:
    """Eş düzlemli iki darbeli dairesel geçiş, büyüklüklerin toplamı."""
    r1, r2 = R_E + h1, R_E + h2
    a = (r1 + r2) / 2
    return abs(math.sqrt(MU / r1) * (math.sqrt(r2 / a) - 1)) + abs(math.sqrt(MU / r2) * (1 - math.sqrt(r1 / a)))


def elipse_indirme(h0: float, hp: float, ha: float) -> float:
    """Daireselden hp × ha elipsine: önce rp'ye inen geçiş, rp'de apojeyi ha'ya indiren ikinci darbe."""
    r0, rp, ra = R_E + h0, R_E + hp, R_E + ha
    v0 = math.sqrt(MU / r0)
    vt_apo = math.sqrt(MU * (2 / r0 - 2 / (r0 + rp)))
    vt_per = math.sqrt(MU * (2 / rp - 2 / (r0 + rp)))
    v_hedef = math.sqrt(MU * (2 / rp - 2 / (rp + ra)))
    return abs(v0 - vt_apo) + abs(vt_per - v_hedef)


def sonlu_yanma_cezasi(yay_rad: float) -> float:
    """Düğüme ortalanmış, sabit itkili düzlem dışı yanmanın etkinlik kaybı: 1 − sin(θ/2)/(θ/2)."""
    x = yay_rad / 2
    return 1 - math.sin(x) / x if x else 0.0


# ---------------------------------------------------------------- sunumun hesapları
def s1a_edinim() -> dict:
    m0, isp = v("m0_1a"), v("isp_1a")
    ip = yakit(m0, v("dv_1a_sma"), isp)
    ideal = ip + yakit(m0 - ip, v("dv_1a_inc"), isp)

    def kayipli(k):  # düzlem içi itki k oranında düşük; düzlem dışı itici etkilenmedi (#134 s.8)
        ip_k = yakit(m0, v("dv_1a_sma"), isp, 1 - k)
        return ip_k + yakit(m0 - ip_k, v("dv_1a_inc"), isp)

    return {"ideal_kg": ideal, "kayipli_alt_kg": kayipli(v("kayip_1a_lo")),
            "kayipli_ust_kg": kayipli(v("kayip_1a_hi")), "ucusta_kg": v("yakit_1a_edinim")}


def s2c_edinim_dv() -> float:
    t3 = sum(math.hypot(*x) if isinstance(x, tuple) else x for x in v("dv_2c_tablo3"))
    return v("dv_2c_dort_buyuk") + v("dv_2c_29ekim") + t3


def isletme_dv_yillik() -> tuple[float, float]:
    """2A oranları: (bakım, kaçınma) m/s/yıl. Kaçınma için makalenin 4 cm/s örnek boyu kullanılır [H]."""
    yil = v("gozlem_son") - v("gozlem_bas")
    return v("occ_sayi") * v("occ_ort") / yil, v("cam_cift") * 2 * v("cam_ornek") / yil


def s2c_butce() -> dict:
    bakim, kacinma = isletme_dv_yillik()
    ufuk = v("ufuk_yil")
    dv = {"edinim": s2c_edinim_dv(), "bakim": bakim * ufuk, "kacinma": kacinma * ufuk}
    defter = sirali(v("m0_2c"), [dv["edinim"], dv["bakim"], dv["kacinma"]], v("isp_2c"))
    kg_edinim, kg_isletme = defter[0][1], defter[1][1] + defter[2][1]
    m_eol = defter[-1][2]
    eol_ideal = yakit(m_eol, v("dv_eol"), v("isp_2c"))
    dv_isletme = dv["bakim"] + dv["kacinma"]
    dv_toplam = dv["edinim"] + dv_isletme + v("dv_eol")
    return {
        "dv_edinim": dv["edinim"], "dv_isletme": dv_isletme, "dv_eol": v("dv_eol"), "dv_toplam": dv_toplam,
        "kg_edinim": kg_edinim, "kg_isletme": kg_isletme, "kg_eol_tasarim": v("yakit_eol"),
        "kg_eol_ideal": eol_ideal, "kg_toplam": kg_edinim + kg_isletme + v("yakit_eol"),
        "bertaraf_payi": v("dv_eol") / dv_toplam,
        "ayrilmamis_kg": v("yakit_2c") - kg_edinim - kg_isletme - v("yakit_eol"),
        "isletme_yillik": dv_isletme / ufuk, "m_eol_bas": m_eol,
    }


def s2c_yanma() -> dict:
    T = periyot(v("h_s2"))
    oran = v("sure_6kasim_s") / T
    yay = 2 * math.pi * oran
    m_once = v("m0_2c") - yakit(v("m0_2c"), s2c_edinim_dv() - _tablo3_sonrasi(), v("isp_2c"))
    dv_6k = math.hypot(*v("dv_2c_tablo3")[1])
    return {"periyot_dk": T / 60, "oran": oran, "yay_derece": math.degrees(yay),
            "ceza": sonlu_yanma_cezasi(yay), "debi_g_s": 1000 / (v("isp_2c") * G0),
            "dv_6kasim": dv_6k, "etkin_kuvvet_N": m_once * dv_6k / v("sure_6kasim_s"),
            "dk_per_m_s_1N": yakit(v("m0_2c"), 1.0, v("isp_2c")) / (1 / (v("isp_2c") * G0)) / 60}


def _tablo3_sonrasi() -> float:
    """6 Kasım ve sonrasındaki Tablo 3 kalemlerinin büyüklük toplamı (6 Kasım'dan önceki kütleyi bulmak için)."""
    t3 = v("dv_2c_tablo3")
    return sum(math.hypot(*x) if isinstance(x, tuple) else x for x in t3[1:])


def tank_ve_kapasite() -> dict:
    gaz_L = v("tank_hacim_L") - v("yakit_1c") / v("rho_n2h4") * 1000
    return {
        "gaz_bas_L": gaz_L,
        "kap_1c": (kapasite(v("m0_1c"), v("yakit_1c"), 205.0), kapasite(v("m0_1c"), v("yakit_1c"), 220.0)),
        "kap_2c": (kapasite(v("m0_2c"), v("yakit_2c"), 200.0), kapasite(v("m0_2c"), v("yakit_2c"), 223.0)),
        "artik_ust_kg": v("tank_kapasite") * (1 - v("tank_bosaltma")),
        "kg_per_m_s_1c": yakit(v("m0_1c"), 1.0, 220.0), "kg_per_m_s_2c": yakit(v("m0_2c"), 1.0, 220.0),
        "yakit_orani_1c": v("yakit_1c") / v("m0_1c"), "yakit_orani_2c": v("yakit_2c") / v("m0_2c"),
    }


def bertaraf() -> dict:
    hoh = hohmann(v("h_s2"), v("h_s2_eol"))
    b = s2c_butce()
    # 47 kg'dan geriye etkin Isp: EOL başlangıç kütlesi 1057–1139 kg aralığında (depo raporu, HESAPLAR §6)
    isp_etkin = [v("dv_eol") / (G0 * -math.log1p(-v("yakit_eol") / m)) for m in (1057.0, 1139.0)]
    return {"hohmann": hoh, "plan_fazlasi": v("dv_eol") / hoh - 1, "isp_etkin": tuple(isp_etkin),
            "s1b_ideal": elipse_indirme(v("h_s1"), v("h_1b_perigee"), v("h_1b_apogee")),
            "kg_eol_ideal": b["kg_eol_ideal"]}


def kaynak_tutarliligi() -> dict:
    """Makalelerdeki sayıların kendi içinde tutarlılığı (bağımsız fizikle)."""
    a1 = R_E + v("h_s1")
    vv = hiz(v("h_s1"))
    # Gauss: dairesel yörüngede küçük Δa için Δv = v Δa / (2a); eğim için Δv = v Δi
    dv_sma_min = vv * abs(v("da_1a")) / (2 * a1)
    dv_inc_min = vv * math.radians(abs(v("di_1a_mdeg")) / 1000)
    # yer izi kayması: ΔP = 1,5 P Δa/a; her turda ekvatorda ω_E R_E ΔP kayar
    P = periyot(v("h_s1"))
    kayma = 1.5 * P * abs(v("da_1a")) / a1 * OMEGA_E * R_E * (86400 / P) / 1000
    olc, ong = v("tork_olculen"), v("tork_ongorulen")
    tork_orani = math.dist(olc, (0, 0, 0)) / math.dist(ong, (0, 0, 0))
    # 1 N, 30 s ve 300 s yanmalar 1A kütlesinde; kayıp oranı uygulanınca
    dv_30 = 1.0 * v("yanma_sonra_s") / v("m0_1a") * (1 - v("kayip_1a_hi"))
    dv_300 = 1.0 * v("yanma_once_s") / v("m0_1a")
    a2 = R_E + v("h_s2")
    dv_2c_min = hiz(v("h_s2")) * abs(v("da_2c")) / (2 * a2)
    return {"dv_sma_min": dv_sma_min, "dv_inc_min": dv_inc_min, "kayma_km_gun": kayma,
            "tork_orani": tork_orani, "dv_30s": dv_30, "dv_300s": dv_300, "dv_2c_min": dv_2c_min}


def itici_omru() -> dict:
    """Bertaraf planının ve tam yakıtın itici yeterlilik ömrüne oranları (dört itici eşit paylaşırsa) [H]."""
    yanma_h = v("eol_yanma_sayi") * v("eol_yanma_s") / 3600
    impuls_2c = v("yakit_2c") * v("ag_isp_nom") * G0
    impuls_1c = v("yakit_1c") * 205.0 * G0
    return {"eol_yanma_h": yanma_h, "eol_yanma_payi": yanma_h / v("ag_toplam_yanma_h"),
            "tek_yanma_payi": v("eol_yanma_s") / 3600 / v("ag_tek_yanma_h"),
            "eol_yakit_itici_basina": v("yakit_eol") / 4, "yakit_gecis_payi": v("yakit_eol") / 4 / v("ag_yakit_gecis"),
            "impuls_2c_itici_basina": impuls_2c / 4, "impuls_2c_payi": impuls_2c / 4 / v("ag_toplam_impuls"),
            "impuls_1c_toplam": impuls_1c, "impuls_1c_tip2_kati": impuls_1c / v("rafael_toplam_impuls")[1]}


def sunum_degerleri() -> dict:
    """Sunumda görünen bütün değerler (yuvarlanmamış)."""
    return {"s1a": s1a_edinim(), "butce": s2c_butce(), "yanma": s2c_yanma(), "tank": tank_ve_kapasite(),
            "bertaraf": bertaraf(), "kaynak": kaynak_tutarliligi(), "omur": itici_omru()}


if __name__ == "__main__":
    import json
    print(json.dumps(sunum_degerleri(), indent=2, ensure_ascii=False))
