"""Reproduce the Sentinel series assignment's *ideal* calculations using SI units.

This is a report companion, not a flight propagator or a hardware qualification model.
No network, project installation or third-party packages are required.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

G0 = 9.80665
MU = 3.986004418e14
EARTH_RADIUS = 6378137.0
EARTH_RATE = 7.2921159e-5


def finite(value: float, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a finite number in the documented SI unit")
    result = float(value)
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError(f"{name} is outside its physical domain")
    return result


def propellant(mass_kg: float, dv_m_s: float, isp_s: float, axial_factor: float = 1) -> float:
    mass = finite(mass_kg, "mass_kg", positive=True)
    dv = finite(dv_m_s, "dv_m_s")
    isp = finite(isp_s, "isp_s", positive=True)
    axial = finite(axial_factor, "axial_factor", positive=True)
    if axial > 1:
        raise ValueError("axial_factor must be in (0, 1]")
    return mass * -math.expm1(-dv / (G0 * isp * axial))


def capacity(mass_kg: float, fuel_kg: float, isp_s: float) -> float:
    mass = finite(mass_kg, "mass_kg", positive=True)
    fuel = finite(fuel_kg, "fuel_kg")
    isp = finite(isp_s, "isp_s", positive=True)
    if fuel >= mass:
        raise ValueError("fuel_kg must be less than mass_kg")
    return -G0 * isp * math.log1p(-fuel / mass)


def orbit(altitude_m: float) -> dict[str, float]:
    radius = EARTH_RADIUS + finite(altitude_m, "altitude_m")
    period = 2 * math.pi * math.sqrt(radius**3 / MU)
    return {"speed_m_s": math.sqrt(MU / radius), "period_s": period,
            "beta_zero_eclipse_s": period * math.asin(EARTH_RADIUS / radius) / math.pi}


def lowering_dv(initial_m: float, perigee_m: float, apogee_m: float) -> float:
    """Two coplanar impulses to a lower ellipse, including circular final orbits."""
    r0, rp, ra = (EARTH_RADIUS + finite(v, "altitude_m")
                  for v in (initial_m, perigee_m, apogee_m))
    if not rp <= ra <= r0:
        raise ValueError("require final perigee <= final apogee <= initial altitude")
    if r0 == rp == ra:
        return 0.0
    transfer_a, final_a = (r0 + rp) / 2, (rp + ra) / 2
    first = math.sqrt(MU / r0) - math.sqrt(MU * (2 / r0 - 1 / transfer_a))
    second = math.sqrt(MU * (2 / rp - 1 / transfer_a)) - math.sqrt(
        MU * (2 / rp - 1 / final_a))
    return first + second


def ledger(mass_kg: float, fuel_kg: float, phases: list[dict],
           isp_s: float, axial_factor: float = 1) -> list[dict]:
    mass = finite(mass_kg, "mass_kg", positive=True)
    fuel = finite(fuel_kg, "fuel_kg")
    if fuel >= mass:
        raise ValueError("fuel must be less than total mass")
    rows, seen = [], set()
    for phase in phases:
        key = phase["id"]
        if key in seen:
            raise ValueError(f"double-counted phase: {key}")
        seen.add(key)
        # Signed burns consume fuel even if a later burn reverses the orbit change.
        vector = phase["dv_components_m_s"]
        if not vector or any(isinstance(v, bool) or not isinstance(v, (int, float))
                             or not math.isfinite(v) for v in vector):
            raise ValueError("dv_components_m_s must be a nonempty finite vector")
        dv = math.hypot(*vector)
        used = propellant(mass, dv, isp_s, axial_factor)
        if used > fuel:
            raise ValueError(f"insufficient propellant at phase {key}")
        row = {"id": key, "dv_m_s": dv, "mass_before_kg": mass,
               "propellant_ideal_kg": used, "mass_after_kg": mass - used,
               "fuel_remaining_ideal_kg": fuel - used}
        duration = phase.get("duration_s")
        if duration is not None:
            duration = finite(duration, "duration_s", positive=True)
            row["duration_s"] = duration
            # Constant effective exhaust speed and constant thrust approximation.
            row["inferred_average_axial_thrust_N"] = used * G0 * isp_s * axial_factor / duration
        rows.append(row)
        mass, fuel = mass - used, fuel - used
    return rows


def drift(dv_m_s: float, elapsed_s: float, altitude_m: float) -> dict[str, float]:
    dv = finite(dv_m_s, "dv_m_s")
    elapsed = finite(elapsed_s, "elapsed_s")
    speed = orbit(altitude_m)["speed_m_s"]
    phase_s = 3 * dv * elapsed / speed
    return {"along_track_time_s": phase_s,
            "ground_track_cross_track_m": EARTH_RATE * EARTH_RADIUS * phase_s}


def thermal_radiator(power_W: float, emissivity: float, temperature_K: float) -> float:
    power = finite(power_W, "power_W")
    eps = finite(emissivity, "emissivity", positive=True)
    temperature = finite(temperature_K, "temperature_K", positive=True)
    if eps > 1:
        raise ValueError("emissivity must be <= 1")
    return power / (eps * 5.670374419e-8 * temperature**4)


def calculate(data: dict) -> dict:
    if data["schema"] != "sentinel-series-v4-SI":
        raise ValueError("unrecognised data schema or unit convention")
    cases = {}
    for key, record in data["capacity_cases"].items():
        if record["mass_vehicle"] != record["fuel_vehicle"]:
            raise ValueError("mass and fuel belong to different vehicles")
        cases[key] = {str(isp): capacity(record["launch_mass_kg"], record["loaded_fuel_kg"], isp)
                      for isp in record["assumed_isp_s"]}
    c = data["capacity_cases"]["S2C"]
    rows = ledger(c["launch_mass_kg"], c["loaded_fuel_kg"], data["s2c_events"], 220)
    acquisition = sum(r["dv_m_s"] for r in rows[:3])
    tandem = sum(r["dv_m_s"] for r in rows[3:])
    s2_eol = lowering_dv(786000, 653000, 653000)
    s2_orbit = orbit(786000)
    scenarios = []
    for isp in (200, 220):
        for axial in (0.8, 1.0):
            scenario = ledger(c["launch_mass_kg"], c["loaded_fuel_kg"],
                              data["s2c_events"], isp, axial)
            scenarios.append({"isp_s_assumed": isp, "axial_factor_assumed": axial,
                              "fuel_used_ideal_kg": c["loaded_fuel_kg"]
                              - scenario[-1]["fuel_remaining_ideal_kg"]})
    return {
        "scope": "Ideal analytical estimates; no measured fuel consumption or flight certification",
        "units": "SI; kg, m, s, N; decimal GB = 1e9 bytes",
        "ideal_capacity_m_s": cases,
        "s2c_ledger_220s_axial1": rows,
        "s2c_acquisition_dv_m_s": acquisition,
        "s2c_table3_dv_m_s": tandem,
        "s2c_total_dv_m_s": acquisition + tandem,
        "s2c_fuel_sensitivity": scenarios,
        "s1_693_to_530km_circular_dv_m_s": lowering_dv(693000, 530000, 530000),
        "s1_693_to_580x610km_dv_m_s": lowering_dv(693000, 580000, 610000),
        "s2_786_to_653km_dv_m_s": s2_eol,
        "s2_73m_s_plan_excess_percent": 100 * (73 / s2_eol - 1),
        "orbit_s1_693km": orbit(693000), "orbit_s2_786km": s2_orbit,
        "illustrative_s2_constant_1400W_load_Wh": 1400 * s2_orbit["period_s"] / 3600,
        "illustrative_s2_1700W_sunlight_generation_Wh": 1700 * (
            s2_orbit["period_s"] - s2_orbit["beta_zero_eclipse_s"]) / 3600,
        "tandem_drift_1cm_s_one_day": drift(0.01, 86400, 786000),
        "tandem_drift_3_85cm_s_one_day": drift(0.0385, 86400, 786000),
        "one_N_220s_mass_flow_kg_s": 1 / (G0 * 220),
        "catalogue_0_44g_s_220s_implied_N": 0.00044 * G0 * 220,
        "rafael_heater_28V_257ohm_W_per_element": 28**2 / 257,
        "downlink_600s_520Mbps_GB": 520e6 * 600 / 8e9,
        "downlink_600s_560Mbps_GB": 560e6 * 600 / 8e9,
        "s2_300GB_emptying_at_560Mbps_s": 300e9 * 8 / 560e6,
        "illustrative_radiator_300W_eps0_8_300K_m2": thermal_radiator(300, 0.8, 300),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1]
                        / "docs/missions/sentinel/SERI_VERILERI.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate(json.loads(args.data.read_text(encoding="utf-8")))
    rendered = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
