"""Run a reference-labelled preliminary pressure-fed blowdown case.

References
----------
NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
"""

from __future__ import annotations

import json

from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    simulate_pressure_fed_blowdown,
)


def main() -> None:
    """Print one converged tabulated blowdown artifact.

    The numerical values are a software example, not a qualification dataset
    or a recommendation for a particular propulsion system.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    rated_point = ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=100.0,
        system_specific_impulse_s=220.0,
        chamber_pressure_pa=1.5e6,
        source="illustrative rated point",
    ).operating_point()
    tank = BlowdownTankModel(
        initial_pressure_pa=2.0e6,
        initial_ullage_volume_m3=0.01,
        initial_propellant_mass_kg=10.0,
        propellant_density_kg_m3=1000.0,
        polytropic_exponent=1.2,
    )
    performance = PressurePerformanceLaw(
        reference_supply_pressure_pa=2.0e6,
        minimum_supply_pressure_pa=1.0e6,
        thrust_pressure_exponent=1.0,
        specific_impulse_pressure_exponent=0.02,
    )
    result = simulate_pressure_fed_blowdown(tank, rated_point, performance)
    print(json.dumps(result.artifact.to_dict(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
