"""Run a sourced pressure-fed propulsion uncertainty ensemble.

References
----------
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
"""

from __future__ import annotations

from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    DistributionKind,
    PressurePerformanceLaw,
    UncertainParameter,
    UncertaintyClass,
    simulate_blowdown_uncertainty,
)


def main() -> None:
    """Print compact statistics for an illustrative seeded ensemble.

    Numerical bounds are software-demonstration values and are not a
    qualification dataset or a recommendation for flight design.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    """

    rated = ConstantPerformanceProvider.monopropellant(
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
    law = PressurePerformanceLaw(
        reference_supply_pressure_pa=2.0e6,
        minimum_supply_pressure_pa=1.0e6,
        maximum_supply_pressure_pa=2.2e6,
        thrust_pressure_exponent=1.0,
        specific_impulse_pressure_exponent=0.02,
    )
    parameters = (
        UncertainParameter(
            name="initial_pressure_pa",
            nominal=2.0e6,
            lower=1.9e6,
            upper=2.1e6,
            distribution=DistributionKind.TRUNCATED_NORMAL,
            uncertainty_class=UncertaintyClass.ALEATORY,
            standard_deviation=40_000.0,
            source="illustrative tank-pressure acceptance band",
            rationale="demonstrate bounded realization-to-realization variation",
        ),
        UncertainParameter(
            name="rated_delivered_thrust_n",
            nominal=100.0,
            lower=97.0,
            upper=103.0,
            distribution=DistributionKind.TRIANGULAR,
            uncertainty_class=UncertaintyClass.EPISTEMIC,
            source="illustrative hot-fire calibration band",
            rationale="demonstrate reducible rated-performance uncertainty",
        ),
    )
    result = simulate_blowdown_uncertainty(
        tank,
        rated,
        law,
        parameters,
        sample_count=32,
        seed=42,
        maximum_mass_step_kg=0.5,
        integration_relative_tolerance=1e-6,
    )
    for statistic in result.statistics:
        print(
            statistic.metric,
            f"mean={statistic.mean:.9g}",
            f"p05={statistic.p05:.9g}",
            f"p95={statistic.p95:.9g}",
        )
    print("sampling-stable", result.convergence.converged)


if __name__ == "__main__":
    main()
