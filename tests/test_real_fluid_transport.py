"""Reference transport properties against NIST WebBook tables and the source papers.

Reference values: ``tests/reference_data/nist_real_fluid_transport.json``
(NIST Chemistry WebBook SRD 69, fetched 2026-10-25, every row with its URL,
and check values printed in the source papers with their locators).

Tolerances, all relative, with the largest deviation measured on
2026-10-25 in brackets. Rows are evaluated at the WebBook ``(T, rho)``
unless stated.

* Nitrogen and oxygen (Lemmon and Jacobsen 2004, the WebBook model):
  viscosity 1e-7 [8.8e-10]; conductivity 1e-5 [2.3e-6, near-critical rows].
* Hydrogen viscosity (Muzny et al. 2013 as transcribed from CoolProp; the
  paper was not opened, so these rows are the evidence for the
  transcription): 1e-9 [4.6e-12]; through ``transport_tp`` 1e-5 [2.5e-6],
  the WebBook gas-constant density difference described in
  test_real_fluids. Hydrogen conductivity (Assael et al. 2011): 5e-5
  [2.2e-5].
* Helium viscosity (TN 1334, Eqs. 1-9): 1e-9 [5.5e-10]; conductivity
  (Hands and Arp via CoolProp, paper not opened): 1e-5 [9.0e-7]; through
  ``transport_tp`` 1e-5 [1.5e-6].
* Dense states (rho > rho_c) of the Olchowy-Sengers fluids whose
  enhancement is below 1 % of ``k`` carry 2e-3 [1.1e-3, hydrogen at 60 K,
  30 MPa]: there the WebBook differs by up to about 1.7e-4 W/(m K), in
  either direction, including rows where the correlation-length bracket of
  Eq. (10) of Lemmon and Jacobsen is negative and the paper prescribes
  ``lambda_c = 0`` (p. 27). How the WebBook treats very small correlation
  lengths is not documented in the papers and was not identified.
* Methane (Friend, Ely and Ingham 1989): the WebBook viscosity is an
  unpublished correlation, so the viscosity carries 5e-2 [4.4e-2]. The
  conductivity carries 3e-2 [1.8e-2, near the critical point] with the
  Friend-Ely-Ingham viscosity and 1e-2 [7.1e-3] with the WebBook viscosity
  supplied to Eq. (18a); the remaining offset tracks the dilute viscosity
  offset (300 K, 0.1 MPa: +0.52 % in k, +0.55 % in eta), consistent with
  the WebBook using its own viscosity in Eqs. (13a) and (18a). That is an
  inference from these numbers, not documented by NIST.
* Paper check values: Lemmon and Jacobsen Table V, 6 digits, 5e-6
  [2.4e-6] (5e-5 [1.8e-5] at the near-critical rows, where the equation of
  state enters); Assael Table 7, 5 digits, 1e-4 [4.6e-5]; Friend-Ely-Ingham
  Table A1, 3-4 digits, 1.5e-3 [1.8e-4]; Table A3 at its tabulated 4-digit
  densities, 3e-3 [1.7e-3]; Muzny erratum check values, 4-5 digits, 2e-5
  [6.1e-6].
* Dilute-gas limit against the NASA CEA fits of
  ``thermochemistry.transport`` (200-1000 K overlap, helium to 1500 K):
  independent fits, expected to agree within their combined uncertainty;
  viscosity 1.5e-2 [1.0e-2], conductivity 5e-2 [3.7e-2].
"""

from __future__ import annotations

import doctest
import json
import unittest
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.thermochemistry.transport import get_species_transport
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid, reference_fluid_model
from rocket_propulsion.thermodynamics.real_fluids import transport as transport_module
from rocket_propulsion.thermodynamics.real_fluids.transport import (
    TRANSPORT_FLUID_KEYS,
    dilute_gas_transport,
    thermal_conductivity_rho_t,
    transport_from_state,
    transport_rho_t,
    transport_tp,
    transport_validity,
    viscosity_rho_t,
)

DATA = json.loads((Path(__file__).parent / "reference_data" / "nist_real_fluid_transport.json")
                  .read_text(encoding="utf-8"))
CHECKS = DATA["primary_check_values"]
HYDROGEN = ("parahydrogen", "normal_hydrogen")
MOLAR_MASS = {"nitrogen": 28.01348e-3, "oxygen": 31.9988e-3, "methane": 16.0428e-3}
# (viscosity, conductivity, conductivity through transport_tp) tolerances by fluid.
TOL = {"nitrogen": (1e-7, 1e-5, 1e-5), "oxygen": (1e-7, 1e-5, 1e-5),
       "parahydrogen": (1e-9, 5e-5, 5e-5), "normal_hydrogen": (1e-9, 5e-5, 5e-5),
       "helium": (1e-9, 1e-5, 1e-5)}
DENSE_SMALL_ENHANCEMENT_TOL = 2e-3


def rel(a: float, b: float) -> float:
    return abs(a / b - 1.0)


def rows(fluid: str) -> list[dict]:
    return DATA["fluids"][fluid]["rows"]


def dense_small_enhancement(fluid: str, result) -> bool:
    """Olchowy-Sengers fluids, rho > rho_c, enhancement below 1 % of k."""

    if fluid in ("helium", "methane"):
        return False
    model = reference_fluid_model(fluid)
    rho_c = model.rho_reduce_mol_m3 * model.molar_mass_kg_mol
    return (result.density_kg_m3 > rho_c
            and result.conductivity_critical_w_m_k < 0.01 * result.thermal_conductivity_w_m_k)


class ReferenceDataTests(unittest.TestCase):
    def test_every_fluid_has_rows_in_every_region(self) -> None:
        self.assertEqual(set(DATA["fluids"]), set(TRANSPORT_FLUID_KEYS))
        for name in TRANSPORT_FLUID_KEYS:
            data = rows(name)
            self.assertGreaterEqual(len(data), 8, name)
            self.assertEqual({r["phase"] for r in data}, {"liquid", "vapor", "supercritical"})
            self.assertGreaterEqual(sum(r["near_critical"] for r in data), 2, name)
            self.assertTrue(all(r["url"].startswith("https://webbook.nist.gov/") for r in data))
        for name in HYDROGEN:
            self.assertGreaterEqual(len(rows(name)), 10)


class WebBookTests(unittest.TestCase):
    """Nitrogen, oxygen, both hydrogens and helium: the WebBook models."""

    def test_viscosity_and_conductivity(self) -> None:
        for name in ("nitrogen", "oxygen", "parahydrogen", "normal_hydrogen", "helium"):
            tol_eta, tol_k, tol_tp = TOL[name]
            for row in rows(name):
                eta_wb = row["eta_uPa_s"] * 1e-6
                with self.subTest(fluid=name, T=row["T_K"], p=row["p_MPa"]):
                    self.assertLess(rel(viscosity_rho_t(name, row["rho_kg_m3"], row["T_K"]),
                                        eta_wb), tol_eta)
                    if not row.get("conductivity_in_range", True):
                        with self.assertRaises(DomainError):
                            transport_rho_t(name, row["rho_kg_m3"], row["T_K"])
                        continue
                    r = transport_rho_t(name, row["rho_kg_m3"], row["T_K"])
                    loose = dense_small_enhancement(name, r)
                    tol = DENSE_SMALL_ENHANCEMENT_TOL if loose else tol_k
                    self.assertLess(rel(r.thermal_conductivity_w_m_k, row["k_W_mK"]), tol)
                    tp = transport_tp(name, row["T_K"], row["p_MPa"] * 1e6)
                    self.assertEqual(tp.phase, row["phase"])
                    self.assertLess(rel(tp.viscosity_pa_s, eta_wb), max(tol_eta, 1e-5))
                    tol = DENSE_SMALL_ENHANCEMENT_TOL if loose else tol_tp
                    self.assertLess(rel(tp.thermal_conductivity_w_m_k, row["k_W_mK"]), tol)

    def test_near_critical_rows_have_large_enhancement(self) -> None:
        for name in ("nitrogen", "oxygen", "parahydrogen", "normal_hydrogen", "helium"):
            for row in (r for r in rows(name) if r["near_critical"] and r["phase"] != "liquid"):
                r = transport_rho_t(name, row["rho_kg_m3"], row["T_K"])
                with self.subTest(fluid=name, T=row["T_K"]):
                    share = r.conductivity_critical_w_m_k / r.thermal_conductivity_w_m_k
                    self.assertGreater(share, 0.05)
                    codes = {w.code for w in r.warnings}
                    self.assertEqual("transport_near_critical" in codes, share > 0.1)

    def test_viscosity_override_is_used(self) -> None:
        row = next(r for r in rows("parahydrogen") if r["near_critical"])
        own = transport_rho_t("parahydrogen", row["rho_kg_m3"], row["T_K"])
        doubled = transport_rho_t("parahydrogen", row["rho_kg_m3"], row["T_K"],
                                  viscosity_pa_s=2.0 * own.viscosity_pa_s)
        self.assertEqual(doubled.viscosity_model, "caller-supplied")
        self.assertIsNone(doubled.viscosity_dilute_pa_s)
        self.assertAlmostEqual(doubled.conductivity_critical_w_m_k,
                               0.5 * own.conductivity_critical_w_m_k, delta=1e-12)
        self.assertIn("COOLPROP", own.viscosity_model)


class MethaneWebBookTests(unittest.TestCase):
    def test_against_webbook(self) -> None:
        for row in rows("methane"):
            eta_wb = row["eta_uPa_s"] * 1e-6
            with self.subTest(T=row["T_K"], p=row["p_MPa"]):
                supplied = thermal_conductivity_rho_t(
                    "methane", row["rho_kg_m3"], row["T_K"], viscosity_pa_s=eta_wb)
                own = thermal_conductivity_rho_t("methane", row["rho_kg_m3"], row["T_K"])
                self.assertLess(rel(supplied, row["k_W_mK"]), 1e-2)
                self.assertLess(rel(own, row["k_W_mK"]), 3e-2)
                v = transport_validity("methane")
                if row["T_K"] <= v.viscosity_t_max_k and row["p_MPa"] * 1e6 <= v.viscosity_p_max_pa:
                    r = transport_tp("methane", row["T_K"], row["p_MPa"] * 1e6)
                    self.assertLess(rel(r.viscosity_pa_s, eta_wb), 5e-2)
                    self.assertAlmostEqual(r.thermal_conductivity_w_m_k, own, delta=1e-6 * own)


class PrimaryCheckValueTests(unittest.TestCase):
    def test_lemmon_jacobsen_table_v(self) -> None:
        table = CHECKS["lemmon_jacobsen_2004_table_v"]
        for name in ("nitrogen", "oxygen"):
            for i, (t, rho_mol_dm3, eta, k) in enumerate(table[name]):
                with self.subTest(fluid=name, T=t, rho=rho_mol_dm3):
                    if rho_mol_dm3 == 0.0:
                        eta_c, k_c = dilute_gas_transport(name, t)
                    else:
                        r = transport_rho_t(name, rho_mol_dm3 * 1e3 * MOLAR_MASS[name], t)
                        eta_c, k_c = r.viscosity_pa_s, r.thermal_conductivity_w_m_k
                    tol = 5e-5 if i == len(table[name]) - 1 else 5e-6
                    self.assertLess(rel(eta_c, eta * 1e-6), tol)
                    self.assertLess(rel(k_c, k * 1e-3), tol)

    def test_assael_table_7(self) -> None:
        for t, rho, k_normal, k_para in CHECKS["assael_2011_table_7"]["rows"]:
            for name, k_ref in (("normal_hydrogen", k_normal), ("parahydrogen", k_para)):
                with self.subTest(fluid=name, T=t, rho=rho):
                    if rho == 0.0:
                        k = dilute_gas_transport(name, t)[1]
                    else:
                        r = transport_rho_t(name, rho, t)
                        self.assertEqual(r.conductivity_critical_w_m_k, 0.0)
                        k = r.thermal_conductivity_w_m_k
                    self.assertLess(rel(k, k_ref * 1e-3), 1e-4)

    def test_muzny_erratum_check_values(self) -> None:
        for t, rho, eta in CHECKS["muzny_2022_erratum_check_values"]["rows"]:
            with self.subTest(T=t, rho=rho):
                for name in HYDROGEN:
                    if rho == 0.0:
                        value = dilute_gas_transport(name, t)[0]
                    else:
                        value = viscosity_rho_t(name, rho, t)
                    self.assertLess(rel(value, eta * 1e-6), 2e-5)

    def test_friend_ely_ingham_tables(self) -> None:
        for t, eta0, k0 in CHECKS["friend_ely_ingham_1989_table_a1"]["rows"]:
            with self.subTest(table="A1", T=t):
                eta_c, k_c = dilute_gas_transport("methane", t)
                self.assertLess(rel(eta_c, eta0 * 1e-6), 1.5e-3)
                self.assertLess(rel(k_c, k0 * 1e-3), 1.5e-3)
        for t, _p, rho_mol_dm3, eta, k in CHECKS["friend_ely_ingham_1989_table_a3"]["rows"]:
            with self.subTest(table="A3", T=t, rho=rho_mol_dm3):
                r = transport_rho_t("methane", rho_mol_dm3 * 1e3 * MOLAR_MASS["methane"], t)
                self.assertLess(rel(r.viscosity_pa_s, eta * 1e-6), 3e-3)
                self.assertLess(rel(r.thermal_conductivity_w_m_k, k * 1e-3), 3e-3)


class DiluteGasConsistencyTests(unittest.TestCase):
    def test_against_cea_fits(self) -> None:
        cases = (("nitrogen", "N2", (200.0, 300.0, 500.0, 700.0, 1000.0)),
                 ("oxygen", "O2", (200.0, 250.0, 300.0)),
                 ("methane", "CH4", (200.0, 300.0, 400.0, 500.0, 625.0)),
                 ("normal_hydrogen", "H2", (200.0, 300.0, 500.0, 700.0, 1000.0)),
                 ("helium", "He", (200.0, 300.0, 500.0, 1000.0, 1500.0)))
        for name, species, temperatures in cases:
            cea = get_species_transport(species)
            for t in temperatures:
                with self.subTest(fluid=name, T=t):
                    eta, k = dilute_gas_transport(name, t)
                    self.assertLess(rel(eta, cea.viscosity_pa_s(t)), 1.5e-2)
                    self.assertLess(rel(k, cea.thermal_conductivity_w_m_k(t)), 5e-2)

    def test_low_pressure_gas_approaches_dilute_limit(self) -> None:
        for name, t in (("nitrogen", 300.0), ("normal_hydrogen", 300.0), ("helium", 300.0)):
            r = transport_tp(name, t, 10.0)
            eta0, k0 = dilute_gas_transport(name, t)
            self.assertLess(rel(r.viscosity_pa_s, eta0), 1e-5)
            self.assertLess(rel(r.thermal_conductivity_w_m_k, k0), 1e-5)


class ApiAndFailClosedTests(unittest.TestCase):
    def test_parts_and_prandtl(self) -> None:
        for name, t, p in (("oxygen", 160.0, 6.0e6), ("parahydrogen", 35.0, 1.5e6),
                           ("helium", 6.0, 0.3e6)):
            r = transport_tp(name, t, p)
            self.assertAlmostEqual(
                r.thermal_conductivity_w_m_k,
                r.conductivity_dilute_w_m_k + r.conductivity_residual_w_m_k
                + r.conductivity_critical_w_m_k, delta=1e-15)
            self.assertAlmostEqual(r.viscosity_pa_s,
                                   r.viscosity_dilute_pa_s + r.viscosity_residual_pa_s,
                                   delta=1e-18)
            self.assertAlmostEqual(r.prandtl_number,
                                   r.cp_j_kg_k * r.viscosity_pa_s / r.thermal_conductivity_w_m_k)
        far = transport_tp("nitrogen", 300.0, 1.0e5)
        self.assertEqual(far.conductivity_critical_w_m_k, 0.0)

    def test_from_state_matches_rho_t(self) -> None:
        state = reference_fluid("methane").state_tp(150.0, 5.0e6)
        a = transport_from_state("methane", state)
        b = transport_rho_t("methane", state.density_kg_m3, 150.0)
        self.assertAlmostEqual(a.thermal_conductivity_w_m_k, b.thermal_conductivity_w_m_k,
                               delta=1e-12)
        with self.assertRaises(InputError):
            transport_from_state("nitrogen", state)

    def test_fail_closed(self) -> None:
        with self.assertRaises(InputError):
            transport_tp("argon", 100.0, 1e6)
        with self.assertRaises(DomainError):
            transport_tp("oxygen", 310.0, 1e6)  # above the 300 K equation-of-state limit
        with self.assertRaises(DomainError):
            transport_tp("parahydrogen", 300.0, 150e6)  # conductivity range ends at 100 MPa
        self.assertGreater(viscosity_rho_t("parahydrogen", 50.0, 300.0), 0.0)
        with self.assertRaises(DomainError):
            viscosity_rho_t("parahydrogen", 120.0, 300.0)  # about 380 MPa > 200 MPa
        with self.assertRaises(DomainError):
            transport_tp("helium", 3.0, 3.0e4)  # conductivity range starts at 3.5 K
        with self.assertRaises(DomainError):
            transport_tp("methane", 450.0, 1e6)  # viscosity range ends at 400 K
        self.assertGreater(thermal_conductivity_rho_t("methane", 5.0, 450.0), 0.0)
        with self.assertRaises(DomainError):
            transport_rho_t("nitrogen", 300.0, 100.0)  # two-phase density at 100 K
        with self.assertRaises(InputError):
            transport_tp("nitrogen", 300.0, 1e5, viscosity_pa_s=-1.0)
        with self.assertRaises(DomainError):
            dilute_gas_transport("oxygen", 400.0)

    def test_module_doctest(self) -> None:
        failures, _ = doctest.testmod(transport_module)
        self.assertEqual(failures, 0)


if __name__ == "__main__":
    unittest.main()
