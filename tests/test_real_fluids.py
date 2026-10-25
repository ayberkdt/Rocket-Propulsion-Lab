"""Real-fluid equations of state against NIST WebBook tables and identities.

Reference values: ``tests/reference_data/nist_real_fluids.json`` (NIST
Chemistry WebBook SRD 69, fetched 2026-10-25, every row with its URL) and the
check values printed in NIST IR 8474 (helium), Tables 3-4, p. 16.

Tolerances. The WebBook evaluates the same reference equations, so agreement
is limited by constants, not by the model:

* 1e-5 relative for density, cp, cv, speed of sound, vapor pressure and
  saturated densities. The largest measured deviations (about 2.6e-6) occur
  for hydrogen and helium, whose WebBook values use R = 8.3144598 J/(mol K)
  instead of the 8.314472 J/(mol K) printed in the source papers; with the
  WebBook constant the agreement is better than 1e-10 (tested below). Oxygen
  cp and cv differ by up to 4.5e-6 because only the printed ideal-gas heat
  capacity Eq. (11) of Stewart et al. (1991) is implemented.
* Joule-Thomson coefficient: 1e-4 K/MPa absolute plus 1e-5 relative.
* Enthalpy and entropy (reference states differ): differences between
  states, 1e-5 cp T and 2e-5 cp.
"""

from __future__ import annotations

import dataclasses
import json
import math
import unittest
from pathlib import Path

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.thermodynamics.real_fluids import (
    REFERENCE_FLUID_KEYS,
    CubicResidual,
    RealFluid,
    constant_cp_ideal,
    cubic_counterpart,
    cubic_model,
    reference_fluid,
    reference_fluid_model,
)

DATA = json.loads((Path(__file__).parent / "reference_data" / "nist_real_fluids.json")
                  .read_text(encoding="utf-8"))
REL_TOL = 1e-5
WEBBOOK_R = 8.3144598  # value that reproduces the WebBook hydrogen/helium tables (see docstring)

_FLUIDS: dict[str, RealFluid] = {}


def fluid(name: str) -> RealFluid:
    if name not in _FLUIDS:
        _FLUIDS[name] = reference_fluid(name)
    return _FLUIDS[name]


def rel(a: float, b: float) -> float:
    return abs(a / b - 1.0)


class NistSinglePhaseTests(unittest.TestCase):
    def test_every_fluid_has_reference_points_in_all_regions(self) -> None:
        self.assertEqual(set(DATA["fluids"]), set(REFERENCE_FLUID_KEYS))
        for name, rows in DATA["fluids"].items():
            phases = {row["phase"] for row in rows["single_phase"]}
            self.assertEqual(phases, {"liquid", "vapor", "supercritical"}, name)
            self.assertGreaterEqual(len(rows["single_phase"]), 15, name)
            self.assertGreaterEqual(len(rows["saturation"]), 6, name)
            self.assertTrue(all(r["url"].startswith("https://webbook.nist.gov/")
                                for r in rows["single_phase"]))

    def test_density_heat_capacities_sound_speed_and_phase(self) -> None:
        for name, rows in DATA["fluids"].items():
            f = fluid(name)
            for row in rows["single_phase"]:
                with self.subTest(fluid=name, T=row["T_K"], p=row["p_MPa"]):
                    s = f.state_tp(row["T_K"], row["p_MPa"] * 1e6)
                    self.assertEqual(s.phase, row["phase"])
                    self.assertLess(rel(s.density_kg_m3, row["rho_kg_m3"]), REL_TOL)
                    self.assertLess(rel(s.cp_j_kg_k, row["cp_J_gK"] * 1e3), REL_TOL)
                    self.assertLess(rel(s.cv_j_kg_k, row["cv_J_gK"] * 1e3), REL_TOL)
                    self.assertLess(rel(s.speed_of_sound_m_s, row["w_m_s"]), REL_TOL)
                    jt = s.joule_thomson_k_pa * 1e6
                    self.assertLess(abs(jt - row["jt_K_MPa"]), 1e-4 + 1e-5 * abs(row["jt_K_MPa"]))

    def test_enthalpy_and_entropy_differences(self) -> None:
        for name, rows in DATA["fluids"].items():
            f = fluid(name)
            states = [(f.state_tp(r["T_K"], r["p_MPa"] * 1e6), r) for r in rows["single_phase"]]
            s0, r0 = states[0]
            for s, r in states[1:]:
                with self.subTest(fluid=name, T=r["T_K"], p=r["p_MPa"]):
                    dh = (s.enthalpy_j_kg - s0.enthalpy_j_kg) - (r["h_kJ_kg"] - r0["h_kJ_kg"]) * 1e3
                    ds = (s.entropy_j_kg_k - s0.entropy_j_kg_k) - (r["s_J_gK"] - r0["s_J_gK"]) * 1e3
                    self.assertLess(abs(dh), 1e-5 * s.cp_j_kg_k * r["T_K"])
                    self.assertLess(abs(ds), 2e-5 * s.cp_j_kg_k)

    def test_webbook_gas_constant_explains_hydrogen_and_helium_residual(self) -> None:
        for name in ("parahydrogen", "normal_hydrogen", "helium"):
            model = dataclasses.replace(reference_fluid_model(name),
                                        gas_constant_j_mol_k=WEBBOOK_R)
            f = RealFluid(model)
            for row in DATA["fluids"][name]["single_phase"]:
                with self.subTest(fluid=name, T=row["T_K"], p=row["p_MPa"]):
                    s = f.state_tp(row["T_K"], row["p_MPa"] * 1e6)
                    self.assertLess(rel(s.density_kg_m3, row["rho_kg_m3"]), 1e-10)
                    self.assertLess(rel(s.speed_of_sound_m_s, row["w_m_s"]), 1e-10)


class NistSaturationTests(unittest.TestCase):
    def test_saturation_pressure_densities_and_latent_heat(self) -> None:
        for name, rows in DATA["fluids"].items():
            f = fluid(name)
            for row in rows["saturation"]:
                liq, vap = row["liquid"], row["vapor"]
                with self.subTest(fluid=name, T=row["T_K"]):
                    sat = f.saturation_t(row["T_K"])
                    self.assertLess(rel(sat.pressure_pa, vap["p_MPa"] * 1e6), REL_TOL)
                    self.assertLess(rel(sat.liquid.density_kg_m3, liq["rho_kg_m3"]), REL_TOL)
                    self.assertLess(rel(sat.vapor.density_kg_m3, vap["rho_kg_m3"]), REL_TOL)
                    hvap = (vap["h_kJ_kg"] - liq["h_kJ_kg"]) * 1e3
                    self.assertLess(rel(sat.enthalpy_of_vaporization_j_kg, hvap), REL_TOL)
                    self.assertLess(rel(sat.liquid.cp_j_kg_k, liq["cp_J_gK"] * 1e3), REL_TOL)
                    self.assertLess(rel(sat.vapor.speed_of_sound_m_s, vap["w_m_s"]), REL_TOL)

    def test_ancillary_seeds_only(self) -> None:
        sat = fluid("nitrogen").saturation_t(80.0)
        self.assertIsNotNone(sat.ancillary_pressure_pa)
        self.assertLess(rel(sat.ancillary_pressure_pa, sat.pressure_pa), 1e-3)
        self.assertNotEqual(sat.ancillary_pressure_pa, sat.pressure_pa)
        self.assertIsNone(fluid("helium").saturation_t(4.0).ancillary_pressure_pa)


class HeliumReportCheckValueTests(unittest.TestCase):
    """NIST IR 8474 Tables 3-4 (p. 16) to half a unit of the last printed digit."""

    CHECK = DATA["ortiz_vega_2023_check_values"]

    @staticmethod
    def half_unit(text: str) -> float:
        decimals = len(text.split(".")[1]) if "." in text else 0
        return 0.5 * 10.0 ** -decimals

    def test_single_phase_check_values_with_si_gas_constant(self) -> None:
        # Table 1 prints R = 8.314472, but Table 3 is reproduced to every printed digit
        # with the SI 2019 value; with 8.314472 the deviation is 1.1e-6 relative.
        model = reference_fluid_model("helium")
        si = RealFluid(dataclasses.replace(model, gas_constant_j_mol_k=MOLAR_GAS_CONSTANT_J_MOL_K))
        default = RealFluid(model)
        for t, rho, p, cv, w in self.CHECK["single_phase"]:
            mass = float(rho) * 1e3 * model.molar_mass_kg_mol
            with self.subTest(T=t, rho=rho):
                s = si.state_rho_t(mass, t)
                self.assertLessEqual(abs(s.pressure_pa / 1e6 - float(p)), self.half_unit(p))
                self.assertLessEqual(abs(s.cv_j_kg_k * model.molar_mass_kg_mol - float(cv)),
                                     self.half_unit(cv))
                self.assertLessEqual(abs(s.speed_of_sound_m_s - float(w)), self.half_unit(w))
                d = default.state_rho_t(mass, t)
                self.assertLess(rel(d.pressure_pa / 1e6, float(p)), 2e-6)

    def test_saturation_check_values(self) -> None:
        f = fluid("helium")
        mm = f.model.molar_mass_kg_mol
        for t, ps, rl, rv, hl, hv in self.CHECK["saturation"]:
            with self.subTest(T=t):
                sat = f.saturation_t(float(t))
                slack = 2e-6
                self.assertLessEqual(abs(sat.pressure_pa / 1e3 - float(ps)),
                                     self.half_unit(ps) + slack * float(ps))
                self.assertLessEqual(abs(sat.liquid.molar_density_mol_m3 / 1e3 - float(rl)),
                                     self.half_unit(rl) + slack * float(rl))
                self.assertLessEqual(abs(sat.vapor.molar_density_mol_m3 / 1e3 - float(rv)),
                                     self.half_unit(rv) + slack * float(rv))
                dh = sat.enthalpy_of_vaporization_j_kg * mm
                self.assertLessEqual(abs(dh - (float(hv) - float(hl))),
                                     2.0 * self.half_unit(hv) + slack * abs(dh))


class ConsistencyTests(unittest.TestCase):
    STATES = (("nitrogen", 300.0, 5.0e6), ("nitrogen", 100.0, 5.0e6), ("methane", 150.0, 3.0e6),
              ("parahydrogen", 25.0, 2.0e6), ("oxygen", 200.0, 10.0e6), ("helium", 10.0, 1.0e6))

    def test_cp_minus_cv_and_maxwell_relation(self) -> None:
        for name, t, p in self.STATES:
            f = fluid(name)
            s = f.state_tp(t, p)
            rho = s.density_kg_m3
            h_t, h_r = 1e-4 * t, 1e-6 * rho
            dp_dt = (f.pressure_unchecked(rho, t + h_t) - f.pressure_unchecked(rho, t - h_t)) / (2 * h_t)
            dp_dr = (f.pressure_unchecked(rho + h_r, t) - f.pressure_unchecked(rho - h_r, t)) / (2 * h_r)
            with self.subTest(fluid=name, T=t, p=p):
                self.assertLess(rel(dp_dt, s.dp_dt_rho_pa_k), 1e-6)
                self.assertLess(rel(dp_dr, s.dp_drho_t_pa_m3_kg), 1e-6)
                cp_cv = t * dp_dt ** 2 / (rho ** 2 * dp_dr)
                self.assertLess(rel(s.cp_j_kg_k - s.cv_j_kg_k, cp_cv), 1e-6)
                # Maxwell relation (ds/drho)_T = -(dp/dT)_rho / rho^2
                ds_dr = (f.state_rho_t(rho + h_r, t).entropy_j_kg_k
                         - f.state_rho_t(rho - h_r, t).entropy_j_kg_k) / (2 * h_r)
                self.assertLess(rel(ds_dr, -dp_dt / rho ** 2), 1e-5)

    def test_caloric_derivatives(self) -> None:
        for name, t, p in self.STATES:
            f = fluid(name)
            s = f.state_tp(t, p)
            h_t, h_p = 1e-4 * t, 1e-4 * p
            up, dn = f.state_tp(t + h_t, p), f.state_tp(t - h_t, p)
            cp = (up.enthalpy_j_kg - dn.enthalpy_j_kg) / (2 * h_t)
            rho = s.density_kg_m3
            cv = (f.state_rho_t(rho, t + h_t).internal_energy_j_kg
                  - f.state_rho_t(rho, t - h_t).internal_energy_j_kg) / (2 * h_t)
            dh_dp = (f.state_tp(t, p + h_p).enthalpy_j_kg
                     - f.state_tp(t, p - h_p).enthalpy_j_kg) / (2 * h_p)
            with self.subTest(fluid=name, T=t, p=p):
                self.assertLess(rel(cp, s.cp_j_kg_k), 1e-6)
                self.assertLess(rel(cv, s.cv_j_kg_k), 1e-6)
                self.assertLess(abs(-dh_dp / s.cp_j_kg_k - s.joule_thomson_k_pa),
                                1e-6 * abs(s.joule_thomson_k_pa) + 1e-12)
                self.assertLess(rel(s.gibbs_energy_j_kg,
                                    s.enthalpy_j_kg - t * s.entropy_j_kg_k), 1e-12)

    def test_saturation_equal_gibbs_and_clausius_clapeyron(self) -> None:
        for name in REFERENCE_FLUID_KEYS:
            f = fluid(name)
            t = 0.5 * (f.model.t_min_k + f.model.t_crit_k)
            sat = f.saturation_t(t)
            rt_m = f.model.gas_constant_j_mol_k * t / f.model.molar_mass_kg_mol
            with self.subTest(fluid=name):
                self.assertLess(abs(sat.liquid.gibbs_energy_j_kg - sat.vapor.gibbs_energy_j_kg),
                                1e-10 * rt_m)
                self.assertLess(rel(sat.liquid.pressure_pa, sat.vapor.pressure_pa), 1e-8)
                h = 1e-4 * t
                dpdt = (f.saturation_t(t + h).pressure_pa - f.saturation_t(t - h).pressure_pa) / (2 * h)
                dv = 1.0 / sat.vapor.density_kg_m3 - 1.0 / sat.liquid.density_kg_m3
                self.assertLess(rel(dpdt, sat.enthalpy_of_vaporization_j_kg / (t * dv)), 1e-6)

    def test_saturation_by_pressure_round_trip(self) -> None:
        f = fluid("methane")
        sat = f.saturation_t(150.0)
        back = f.saturation_p(sat.pressure_pa)
        self.assertAlmostEqual(back.temperature_k, 150.0, delta=1e-8)

    def test_departure_functions(self) -> None:
        f = fluid("nitrogen")
        t, p = 150.0, 8.0e6
        dep = f.departures(t, p)
        mm = f.model.molar_mass_kg_mol
        low = f.state_tp(t, 1.0)
        here = f.state_tp(t, p)
        self.assertLess(abs(dep.enthalpy_j_mol - (here.enthalpy_j_kg - low.enthalpy_j_kg) * mm),
                        1e-4 * abs(dep.enthalpy_j_mol))
        self.assertAlmostEqual(dep.gibbs_energy_j_mol,
                               dep.enthalpy_j_mol - t * dep.entropy_j_mol_k, delta=1e-9)
        # ln(phi) = integral of (Z - 1)/p dp along the isotherm (trapezoid, 4000 steps)
        n, total, prev = 4000, 0.0, 0.0
        for i in range(1, n + 1):
            pi = p * i / n
            z = f.state_tp(t, pi).compressibility_factor
            cur = (z - 1.0) / pi
            total += 0.5 * (cur + prev) * (p / n) if i > 1 else cur * (p / n)
            prev = cur
        self.assertAlmostEqual(dep.ln_fugacity_coefficient, total, delta=2e-4)


class ValidityTests(unittest.TestCase):
    def test_out_of_range_states_are_refused(self) -> None:
        f = fluid("oxygen")
        with self.assertRaises(DomainError):
            f.state_tp(50.0, 1e6)  # below triple point 54.361 K
        with self.assertRaises(DomainError):
            f.state_tp(301.0, 1e6)  # above 300 K (Stewart et al. 1991 range)
        with self.assertRaises(DomainError):
            f.state_tp(200.0, 81e6)  # above 80 MPa
        with self.assertRaises(DomainError):
            f.state_tp(200.0, 0.0)
        with self.assertRaises(DomainError):
            f.state_tp(math.nan, 1e6)
        with self.assertRaises(DomainError):
            f.saturation_t(160.0)  # above T_c
        with self.assertRaises(DomainError):
            f.state_rho_t(500.0, 120.0)  # inside the two-phase dome
        with self.assertRaises(DomainError):
            f.saturation_p(6e6)  # above p_c
        with self.assertRaises(InputError):
            reference_fluid("water")

    def test_saturated_pressure_is_refused_for_single_phase_call(self) -> None:
        f = fluid("nitrogen")
        sat = f.saturation_t(90.0)
        with self.assertRaises(DomainError):
            f.state_tp(90.0, sat.pressure_pa)

    def test_near_critical_saturation_fails_closed(self) -> None:
        f = fluid("nitrogen")
        try:
            sat = f.saturation_t(f.model.t_crit_k * (1.0 - 1e-9))
        except ConvergenceError:
            return
        self.assertGreater(sat.liquid.density_kg_m3, sat.vapor.density_kg_m3)


class CubicTests(unittest.TestCase):
    def test_residual_derivatives_match_finite_differences(self) -> None:
        for kind in ("peng_robinson", "soave_redlich_kwong"):
            res = cubic_counterpart("methane", kind).model.residual
            tau, delta, h = 1.3, 0.8, 1e-6
            d = res.derivatives(tau, delta)
            dd = [(a - b) / (2 * h) for a, b in zip(res.derivatives(tau, delta + h),
                                                     res.derivatives(tau, delta - h), strict=True)]
            dt = [(a - b) / (2 * h) for a, b in zip(res.derivatives(tau + h, delta),
                                                     res.derivatives(tau - h, delta), strict=True)]
            with self.subTest(kind=kind):
                self.assertAlmostEqual(dd[0], d[1], delta=1e-7 * abs(d[1]) + 1e-10)
                self.assertAlmostEqual(dd[1], d[2], delta=1e-7 * abs(d[2]) + 1e-10)
                self.assertAlmostEqual(dt[0], d[3], delta=1e-7 * abs(d[3]) + 1e-10)
                self.assertAlmostEqual(dt[3], d[4], delta=1e-7 * abs(d[4]) + 1e-10)
                self.assertAlmostEqual(dt[1], d[5], delta=1e-7 * abs(d[5]) + 1e-10)

    def test_pressure_matches_pressure_explicit_form(self) -> None:
        for kind, (d1, d2) in (("peng_robinson", (1 + 2 ** 0.5, 1 - 2 ** 0.5)),
                               ("soave_redlich_kwong", (1.0, 0.0))):
            f = cubic_counterpart("nitrogen", kind)
            res: CubicResidual = f.model.residual
            t, rho_molar = 150.0, 8000.0
            v = 1.0 / rho_molar
            a, _, _ = res.attraction(t)
            b, r = res.b_m3_mol, res.gas_constant_j_mol_k
            p_direct = r * t / (v - b) - a / ((v + d1 * b) * (v + d2 * b))
            p_model = f.pressure_unchecked(rho_molar * f.model.molar_mass_kg_mol, t)
            with self.subTest(kind=kind):
                self.assertLess(rel(p_model, p_direct), 1e-12)

    def test_critical_point_is_reproduced(self) -> None:
        for kind in ("peng_robinson", "soave_redlich_kwong"):
            f = cubic_counterpart("nitrogen", kind)
            m = f.model
            tc, mm = m.t_crit_k, m.molar_mass_kg_mol
            # find the inflection of the critical isotherm (dp/drho minimum)
            best = min(((f.pressure_unchecked(r * (1 + 1e-6), tc)
                         - f.pressure_unchecked(r * (1 - 1e-6), tc)) / (2e-6 * r), r)
                       for r in (x * 0.5 for x in range(200, 1400)))
            with self.subTest(kind=kind):
                self.assertLess(abs(best[0]), 1e-3 * m.p_crit_pa / best[1])
                self.assertLess(rel(f.pressure_unchecked(best[1], tc), m.p_crit_pa), 2e-3)
                self.assertGreater(best[1] / mm, 0.0)

    def test_cubic_accuracy_against_nist_and_warning(self) -> None:
        rows = DATA["fluids"]["nitrogen"]
        for kind in ("peng_robinson", "soave_redlich_kwong"):
            f = cubic_counterpart("nitrogen", kind)
            for row in rows["single_phase"]:
                if row["phase"] == "vapor" and row["p_MPa"] <= 1.0:
                    s = f.state_tp(row["T_K"], row["p_MPa"] * 1e6)
                    self.assertLess(rel(s.density_kg_m3, row["rho_kg_m3"]), 0.02)
                    self.assertIn("cubic_eos_lower_accuracy", {w.code for w in s.warnings})
            for row in rows["saturation"]:
                sat = f.saturation_t(row["T_K"])
                self.assertLess(rel(sat.pressure_pa, row["vapor"]["p_MPa"] * 1e6), 0.05)

    def test_arbitrary_fluid_and_alpha_function_limit(self) -> None:
        ref = reference_fluid_model("nitrogen")
        model = cubic_model("soave_redlich_kwong", name="test-n2", t_crit_k=ref.t_crit_k,
                            p_crit_pa=ref.p_crit_pa, acentric_factor=0.0372,
                            molar_mass_kg_mol=ref.molar_mass_kg_mol,
                            ideal=constant_cp_ideal(3.5), t_min_k=70.0, t_max_k=3000.0,
                            p_max_pa=1e8)
        f = RealFluid(model)
        s = f.state_tp(300.0, 100.0)
        r_m = model.gas_constant_j_mol_k / model.molar_mass_kg_mol
        self.assertLess(rel(s.cp_j_kg_k, 3.5 * r_m), 1e-4)
        with self.assertRaises(DomainError):
            f.state_tp(2500.0, 1e5)  # 1 + m (1 - sqrt(T/Tc)) < 0 for this fluid
        with self.assertRaises(InputError):
            constant_cp_ideal(1.0)


if __name__ == "__main__":
    unittest.main()
