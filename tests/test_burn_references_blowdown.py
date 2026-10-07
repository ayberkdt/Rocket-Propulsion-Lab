"""Reference, documentation, and physics checks for pressure-fed blowdown."""

from __future__ import annotations

import inspect
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    BLOWDOWN_REFERENCE_IDS,
    BURN_REFERENCES,
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    PressurePerformanceSample,
    PressurePerformanceTable,
    blowdown,
    burn_reference,
    calibration,
    cluster,
    discrepancy,
    ensemble,
    ensemble_serialization,
    feed_line,
    feed_system,
    multivariate_calibration,
    performance_surface,
    propagator,
    references,
    response,
    scenarios,
    sequence,
    simulate_pressure_fed_blowdown,
    tabulated,
    uncertainty,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2


class BurnReferenceTests(unittest.TestCase):
    """Keep the reviewed source catalog and public transient API auditable."""

    def test_reference_registry_is_unique_and_uses_primary_source_urls(self) -> None:
        identifiers = tuple(reference.reference_id for reference in BURN_REFERENCES)
        self.assertEqual(len(identifiers), len(set(identifiers)))
        self.assertTrue(set(BLOWDOWN_REFERENCE_IDS).issubset(identifiers))
        for reference in BURN_REFERENCES:
            self.assertIn(
                reference.organization,
                {
                    "NASA",
                    "NASA Glenn Research Center",
                    "CCSDS",
                    "NIST",
                    "JCGM",
                    "Orekit Project",
                },
            )
            self.assertTrue(
                reference.url.startswith("https://ntrs.nasa.gov/")
                or reference.url.startswith("https://www1.grc.nasa.gov/")
                or reference.url.startswith("https://www.grc.nasa.gov/")
                or reference.url.startswith("https://standards.nasa.gov/")
                or reference.url.startswith("https://ccsds.org/")
                or reference.url.startswith("https://www.nist.gov/")
                or reference.url.startswith("https://www.itl.nist.gov/")
                or reference.url.startswith("https://www.bipm.org/")
                or reference.url.startswith("https://www.orekit.org/")
            )
            self.assertEqual(burn_reference(reference.reference_id), reference)
        with self.assertRaises(InputError):
            burn_reference("unreviewed-source")

    def test_public_transient_types_and_functions_carry_references(self) -> None:
        modules = (
            references,
            response,
            sequence,
            cluster,
            tabulated,
            blowdown,
            calibration,
            discrepancy,
            multivariate_calibration,
            performance_surface,
            uncertainty,
            scenarios,
            ensemble,
            ensemble_serialization,
            feed_line,
            feed_system,
            propagator,
        )
        missing: list[str] = []
        for module in modules:
            for name, value in inspect.getmembers(module):
                if name.startswith("_") or getattr(value, "__module__", None) != module.__name__:
                    continue
                if (inspect.isclass(value) or inspect.isfunction(value)) and "References" not in (
                    inspect.getdoc(value) or ""
                ):
                    missing.append(f"{module.__name__}.{name}")
        self.assertEqual(missing, [], f"Public transient API lacks references: {missing}")


class PressureFedBlowdownTests(unittest.TestCase):
    """Check pressure decay, calibrated scaling, cutoff, and integration closure."""

    def setUp(self) -> None:
        self.tank = BlowdownTankModel(
            initial_pressure_pa=2_000_000.0,
            initial_ullage_volume_m3=0.01,
            initial_propellant_mass_kg=10.0,
            propellant_density_kg_m3=1_000.0,
            polytropic_exponent=1.2,
        )
        self.point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=220.0,
            chamber_pressure_pa=1_500_000.0,
            source="unit-test reference point",
        ).operating_point()
        self.law = PressurePerformanceLaw(
            reference_supply_pressure_pa=2_000_000.0,
            minimum_supply_pressure_pa=1_000_000.0,
            maximum_supply_pressure_pa=2_000_000.0,
            thrust_pressure_exponent=1.0,
            specific_impulse_pressure_exponent=0.02,
            chamber_pressure_exponent=1.0,
        )

    def test_polytropic_pressure_law_and_inverse(self) -> None:
        self.assertAlmostEqual(self.tank.pressure_pa(10.0), 2_000_000.0)
        empty_pressure = 2_000_000.0 * (0.01 / 0.02) ** 1.2
        self.assertAlmostEqual(self.tank.pressure_pa(0.0), empty_pressure)
        remaining = self.tank.remaining_mass_at_pressure(1_000_000.0)
        self.assertAlmostEqual(self.tank.pressure_pa(remaining), 1_000_000.0)

    def test_pressure_scaling_preserves_system_isp_mass_identity(self) -> None:
        reference = self.law.scale_operating_point(self.point, 2_000_000.0)
        lower = self.law.scale_operating_point(self.point, 1_000_000.0)
        self.assertAlmostEqual(reference.delivered_thrust_n, 100.0)
        self.assertAlmostEqual(lower.delivered_thrust_n, 50.0)
        self.assertLess(lower.system_specific_impulse_s, reference.system_specific_impulse_s)
        self.assertAlmostEqual(
            lower.delivered_thrust_n,
            STANDARD_GRAVITY_M_S2
            * lower.system_specific_impulse_s
            * lower.total_tank_flow_kg_s,
        )
        self.assertAlmostEqual(
            sum(stream.mass_flow_kg_s for stream in lower.streams),
            lower.total_tank_flow_kg_s,
        )
        with self.assertRaises(DomainError):
            self.law.scale_operating_point(self.point, 999_999.0)
        with self.assertRaisesRegex(DomainError, "Reference supply pressure"):
            PressurePerformanceLaw(
                reference_supply_pressure_pa=2_000_000.0,
                minimum_supply_pressure_pa=1_000_000.0,
                maximum_supply_pressure_pa=1_500_000.0,
            )

    def test_measured_pressure_map_interpolates_without_extrapolation(self) -> None:
        table = PressurePerformanceTable(
            samples=(
                PressurePerformanceSample(1_000_000.0, 48.0, 215.0, 700_000.0),
                PressurePerformanceSample(1_500_000.0, 76.0, 218.0, 1_100_000.0),
                PressurePerformanceSample(2_000_000.0, 100.0, 220.0, 1_500_000.0),
            ),
            source_id="hot-fire-series-A",
        )
        point = table.scale_operating_point(self.point, 1_250_000.0)
        self.assertAlmostEqual(point.delivered_thrust_n, 62.0)
        self.assertAlmostEqual(point.system_specific_impulse_s, 216.5)
        self.assertAlmostEqual(point.chamber_pressure_pa, 900_000.0)
        self.assertAlmostEqual(
            point.delivered_thrust_n,
            STANDARD_GRAVITY_M_S2
            * point.system_specific_impulse_s
            * point.total_tank_flow_kg_s,
        )
        with self.assertRaises(DomainError):
            table.scale_operating_point(self.point, 2_000_001.0)

    def test_blowdown_builds_monotonic_reference_labelled_history(self) -> None:
        result = simulate_pressure_fed_blowdown(
            self.tank,
            self.point,
            self.law,
            maximum_mass_step_kg=0.5,
            relative_tolerance=1e-6,
        )
        self.assertEqual(result.cutoff_reason, "minimum_supply_pressure")
        self.assertEqual(result.artifact.reference_ids, BLOWDOWN_REFERENCE_IDS)
        self.assertEqual(result.reference_ids[-1], "NASA-20040000363")
        self.assertGreater(len(result.samples), 2)
        for previous, current in zip(result.samples, result.samples[1:], strict=False):
            self.assertLess(previous.time_s, current.time_s)
            self.assertGreater(previous.remaining_propellant_mass_kg, current.remaining_propellant_mass_kg)
            self.assertGreater(previous.supply_pressure_pa, current.supply_pressure_pa)
            self.assertGreater(previous.delivered_thrust_n, current.delivered_thrust_n)
        self.assertAlmostEqual(result.samples[-1].supply_pressure_pa, 1_000_000.0)
        self.assertLessEqual(
            result.evidence.mass_closure_error_kg,
            max(1e-11, result.evidence.consumed_propellant_kg * 1e-6),
        )
        self.assertLessEqual(
            result.evidence.impulse_closure_error_n_s,
            max(1e-11, result.evidence.delivered_impulse_n_s * 1e-6),
        )
        self.assertTrue(result.artifact.artifact_hash)


if __name__ == "__main__":
    unittest.main()
