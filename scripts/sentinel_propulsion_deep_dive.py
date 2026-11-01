"""Reproducible academic Sentinel calculations. SI; no flight performance map.

All non-published scenario parameters are explicitly recorded in the output.
The finite burn is a planar two-body educational model, not flight software.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sentinel_series_budget import EARTH_RADIUS, G0, MU, finite, ledger, orbit, propellant


def cutoff_residual(tank_m3: float, load_kg: float, density_kg_m3: float,
                    initial_bar_abs: float, inlet_min_bar_abs: float,
                    loss_bar: float, exponent: float) -> dict[str, float]:
    """PV**n = const. Loss is a specified cutoff allowance, not a pipe model."""
    values = [finite(v, name, positive=True) for v, name in (
        (tank_m3, "tank_m3"), (load_kg, "load_kg"),
        (density_kg_m3, "density_kg_m3"), (initial_bar_abs, "initial_bar_abs"),
        (inlet_min_bar_abs, "inlet_min_bar_abs"), (exponent, "exponent"))]
    volume, load, density, pressure, minimum, n = values
    loss = finite(loss_bar, "loss_bar")
    gas0 = volume - load / density
    if gas0 <= 0 or pressure <= minimum + loss or not 1 <= n <= 5 / 3:
        raise ValueError("require ullage, initial feed margin and n in [1,5/3]")
    gas_end = min(volume, gas0 * (pressure / (minimum + loss)) ** (1 / n))
    residual = max(0.0, density * (volume - gas_end))
    return {"initial_gas_L": gas0 * 1000, "cutoff_gas_L": gas_end * 1000,
            "residual_kg": residual, "accessible_kg": load - residual,
            "empty_pressure_bar_abs": pressure * (gas0 / volume) ** n}


def burn(mass_kg: float, dv_m_s: float, isp_s: float, engines: int,
         force_each_N: float, cant_deg: float, duty: float = 1) -> dict[str, float]:
    """Equal symmetric thrusters, common duty and Isp, zero transverse resultant."""
    if isinstance(engines, bool) or not isinstance(engines, int) or engines < 1:
        raise ValueError("engines must be a positive integer")
    angle = finite(cant_deg, "cant_deg")
    duty = finite(duty, "duty", positive=True)
    if angle >= 90 or duty > 1:
        raise ValueError("require cant < 90 deg and duty <= 1")
    force = finite(force_each_N, "force_each_N", positive=True)
    isp = finite(isp_s, "isp_s", positive=True)
    eta = math.cos(math.radians(angle))
    used = propellant(mass_kg, dv_m_s, isp, eta)
    scalar = engines * force * duty
    rate = scalar / (G0 * isp)
    return {"axial_factor": eta, "effective_isp_s": isp * eta,
            "scalar_thrust_N": scalar, "axial_thrust_N": scalar * eta,
            "mass_flow_kg_s": rate, "propellant_kg": used,
            "duration_s": used / rate, "total_scalar_impulse_Ns": used * G0 * isp}


def finite_burn(altitude_m: float, mass_kg: float, target_scalar_dv_m_s: float,
                force_N: float, isp_s: float, step_s: float = 1) -> dict[str, float]:
    """Velocity-antiparallel constant-thrust burn with mass loss and RK4.

    Integrate specific work and angular momentum from thrust independently of
    deriving those invariants from the final state. Includes central gravity only.
    """
    h = finite(altitude_m, "altitude_m")
    mass = finite(mass_kg, "mass_kg", positive=True)
    dv = finite(target_scalar_dv_m_s, "target_scalar_dv_m_s", positive=True)
    force = finite(force_N, "force_N", positive=True)
    isp = finite(isp_s, "isp_s", positive=True)
    step = finite(step_s, "step_s", positive=True)
    if dv >= orbit(h)["speed_m_s"] / 10 or step > 5:
        raise ValueError("screening model requires small burns and step <= 5 s")
    used = propellant(mass, dv, isp)
    rate = force / (G0 * isp)
    duration = used / rate
    count = math.ceil(duration / step)
    if count > 100000:
        raise ValueError("integration budget exceeded")
    dt = duration / count
    radius = EARTH_RADIUS + h
    speed = math.sqrt(MU / radius)
    # x, y, vx, vy, m, integral(a_thrust dot v), integral(r cross a_thrust)
    state = [radius, 0., 0., speed, mass, 0., 0.]

    def rhs(s: list[float]) -> list[float]:
        x, y, vx, vy, m = s[:5]
        r, v = math.hypot(x, y), math.hypot(vx, vy)
        ax, ay = -force / m * vx / v, -force / m * vy / v
        return [vx, vy, -MU * x / r**3 + ax, -MU * y / r**3 + ay,
                -rate, ax * vx + ay * vy, x * ay - y * ax]

    for _ in range(count):
        k1 = rhs(state)
        k2 = rhs([a + dt / 2 * b for a, b in zip(state, k1)])
        k3 = rhs([a + dt / 2 * b for a, b in zip(state, k2)])
        k4 = rhs([a + dt * b for a, b in zip(state, k3)])
        state = [a + dt / 6 * (b + 2 * c + 2 * d + e)
                 for a, b, c, d, e in zip(state, k1, k2, k3, k4)]
    x, y, vx, vy, mf, work, dh = state
    r = math.hypot(x, y)
    energy = (vx**2 + vy**2) / 2 - MU / r
    momentum = x * vy - y * vx
    a = -MU / (2 * energy)
    # Eccentricity vector avoids cancellation in 1+2*E*h*h/mu**2.
    ex, ey = vy * momentum / MU - x / r, -vx * momentum / MU - y / r
    eccentricity = math.hypot(ex, ey)
    vi = speed - dv
    impulse_a = 1 / (2 / radius - vi**2 / MU)
    impulse_e = 1 - radius * vi**2 / MU
    return {"duration_s": duration, "steps": count, "step_s": dt,
            "propellant_kg": used, "final_mass_kg": mf,
            "perigee_km": (a * (1 - eccentricity) - EARTH_RADIUS) / 1000,
            "apogee_km": (a * (1 + eccentricity) - EARTH_RADIUS) / 1000,
            "eccentricity": eccentricity,
            "impulsive_perigee_km": (impulse_a * (1 - impulse_e) - EARTH_RADIUS) / 1000,
            "impulsive_apogee_km": h / 1000,
            "energy_work_residual_J_kg": energy + MU / (2 * radius) - work,
            "angular_momentum_residual_m2_s": momentum - radius * speed - dh,
            "mass_residual_kg": mf - (mass - used),
            "burn_orbit_fraction": duration / orbit(h)["period_s"]}


def calculate() -> dict:
    data = json.loads((Path(__file__).resolve().parents[1]
                       / "docs/missions/sentinel/SERI_VERILERI.json").read_text(encoding="utf-8"))
    rows = ledger(1143, 133, data["s2c_events"], 220)
    dv_s2 = math.fsum(r["dv_m_s"] for r in rows)
    feed = [{"P0_bar_abs_ASSUMED": p, "n_ASSUMED": n,
             **cutoff_residual(.204, 154, 1004, p, 5.5, .5, n)}
            for p in (22., 24.5) for n in (1., 1.2, 5 / 3)]
    thrust_cases = [{"vehicle": name, "cant_deg_ASSUMED": cant,
                     "force_each_N_ASSUMED": force,
                     **burn(mass, 1, 220 if name == "S2C" else 205,
                            4 if name == "S2C" else 1, force, cant)}
                    for name, mass in (("S1C", 2186), ("S2C", 1143))
                    for cant in (0., 20., 30.) for force in (.5, 1.)]
    uncertainty = [{"mass_kg_ASSUMED": m, "isp_s_ASSUMED": isp,
                    "cant_deg_ASSUMED": cant,
                    "propellant_kg": propellant(m, dv_s2, isp, math.cos(math.radians(cant)))}
                   for m in (1140, 1143, 1146) for isp in (200, 220, 223)
                   for cant in (0., 20., 30.)]
    # Published forces are simulator inputs reflecting a specific S1B state,
    # not a complete mission-long map or guaranteed C/D performance.
    s1b_forces = [.57, .56, .47, .58]
    allocations = []
    for name, load, isp, cap in (("S1C", 154, 205, 60000), ("S2C", 133, 220, 135000)):
        for count in (1, 4):
            impulse_each = load * G0 * isp / count
            allocations.append({"vehicle": name, "engines_equal_share_ASSUMED": count,
                                "load_kg": load, "Isp_s_ASSUMED": isp,
                                "impulse_each_Ns": impulse_each,
                                "catalogue_limit_Ns": cap,
                                "impulse_fraction": impulse_each / cap,
                                "mass_each_kg": load / count})
    drag = [{"vehicle": name, "rho_kg_m3_ASSUMED": rho,
             "Cd_ASSUMED": 2.2, "area_m2_ASSUMED": area,
             "dv_per_year_m_s": .5 * rho * 2.2 * area / mass
             * orbit(alt)["speed_m_s"]**2 * (365.25 * 86400)}
            for name, mass, alt in (("S1C", 2186, 693000), ("S2C", 1143, 786000))
            for rho in (1e-14, 1e-13, 1e-12) for area in (10., 20.)]
    return {
        "scope": "V5 academic derivations. ASSUMED labels are scenarios, not flight data.",
        "density_basis": "NASA/TM-2013-216560 Table 1: 1004 kg/m3 at 300 K",
        "feed_scope": "S1 family 204 L and C load combined for sensitivity. Initial pressures, n and 0.5 bar loss assumed. No S2 tank volume inferred.",
        "feed_cutoff_sensitivity": feed,
        "one_m_s_burn_sensitivity": thrust_cases,
        "s2c_delta_v_sum_m_s": dv_s2,
        "s2c_bounded_scenarios": uncertainty,
        "s2c_bounded_scenarios_min_max_kg": [min(r["propellant_kg"] for r in uncertainty),
                                             max(r["propellant_kg"] for r in uncertainty)],
        "s2c_ledger_220s_eta1": rows,
        "s1b_published_simulator_forces_N": s1b_forces,
        "s1b_scalar_force_sum_N": math.fsum(s1b_forces),
        "s1b_force_spread_percent_of_mean": 100 * (max(s1b_forces)-min(s1b_forces))
        / (math.fsum(s1b_forces)/4),
        "life_allocation_screen": allocations,
        "MIB_velocity_quanta_m_s": {"S1C_0.008Ns": .008 / 2186,
                                     "S2C_0.01Ns": .01 / 1143,
                                     "S2C_0.043Ns": .043 / 1143},
        "plane_change_0.01deg_m_s": {name: 2 * orbit(alt)["speed_m_s"] * math.sin(math.radians(.01)/2)
                                    for name, alt in (("S1", 693000), ("S2", 786000))},
        "drag_scenarios": drag,
        "finite_burn_assumptions": {"altitude_m": 693000, "mass_kg": 2100,
                                     "scalar_dv_m_s": .7, "isp_s": 205,
                                     "guidance": "instantaneously anti-velocity; no J2/drag/plume/attitude"},
        "finite_burn_cases": [{"force_N_ASSUMED": f,
                               **finite_burn(693000, 2100, .7, f, 205)} for f in (1., 2.18, 4.)],
        "finite_burn_half_step": finite_burn(693000, 2100, .7, 2.18, 205, .5),
        "reserve_example": {"current_mass_kg_ASSUMED": 2100, "remaining_fuel_kg_ASSUMED": 70,
                            "withheld_fuel_kg_ASSUMED": 60, "annual_dv_m_s_ASSUMED": 1,
                            "Isp_s_ASSUMED": 205,
                            "ideal_years": -205 * G0 * math.log1p(-10 / 2100)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(calculate(), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
