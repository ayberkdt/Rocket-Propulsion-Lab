"""Hot-fire calibration fields over operating-condition grids.

Oracles: exact reproduction of the node calibrations (``calibrate_hot_fire_scale``
is tested separately in ``test_burn_calibration``); multilinear weights from
the coordinates; the linear-envelope and independent-node uncertainty forms of
JCGM 100:2008 5.2.2 NOTE 1 and equation 10, written out here by hand.
"""

from __future__ import annotations

import json
import math
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns.calibration import (
    CalibrationScope,
    HotFireScaleMeasurement,
    calibrate_hot_fire_scale,
)
from rocket_propulsion.propulsion.burns.calibration_field import (
    CalibrationFieldNode,
    HotFireCalibrationField,
    apply_calibration_fields,
    calibrate_hot_fire_field,
)
from rocket_propulsion.propulsion.burns.ensemble import ScenarioConditionalUncertainty
from rocket_propulsion.propulsion.burns.performance_surface import (
    PerformanceAxis,
    PerformanceScalePoint,
    RectilinearPerformanceSurface,
)
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256
from rocket_propulsion.propulsion.burns.uncertainty import UncertaintyClass

PRESSURES = (1_000_000.0, 2_000_000.0, 3_000_000.0)
RATIOS = (5.5, 6.0)
AXES = (
    PerformanceAxis("chamber_pressure_pa", "Pa", PRESSURES),
    PerformanceAxis("mixture_ratio", "1", RATIOS),
)


def condition_id(i: int, j: int) -> str:
    return f"pc{i}-mr{j}"


def stratum(
    parameter: str, i: int, j: int, *, condition: str | None = None
) -> tuple[HotFireScaleMeasurement, ...]:
    """Three engines, unbalanced, with offsets that grow with pressure index."""
    base = (0.96 if parameter == "thrust_scale" else 0.98) + 0.01 * i + 0.004 * j
    spread = 0.004 + 0.002 * i + 0.001 * j
    runs = {"E1": (-0.001, 0.001), "E2": (-0.0015, 0.0, 0.0015), "E3": (0.0005, -0.0005)}
    offsets = {"E1": -spread, "E2": 0.3 * spread, "E3": spread}
    return tuple(
        HotFireScaleMeasurement(
            engine_id=engine,
            run_id=f"r{k}",
            condition_id=condition or condition_id(i, j),
            value=base + offsets[engine] + deviation,
            standard_uncertainty=0.0004 + 0.0001 * k,
            source="calibration-field fixture",
        )
        for engine, deviations in runs.items()
        for k, deviation in enumerate(deviations)
    )


def all_measurements(parameter: str = "thrust_scale") -> tuple[HotFireScaleMeasurement, ...]:
    return tuple(
        item
        for i in range(len(PRESSURES))
        for j in range(len(RATIOS))
        for item in stratum(parameter, i, j)
    )


def mapping() -> dict[str, tuple[float, float]]:
    return {
        condition_id(i, j): (PRESSURES[i], RATIOS[j])
        for i in range(len(PRESSURES))
        for j in range(len(RATIOS))
    }


def make_field(parameter: str = "thrust_scale") -> HotFireCalibrationField:
    return calibrate_hot_fire_field(
        parameter,
        AXES,
        mapping(),
        all_measurements(parameter),
        source="qualification campaign field fixture",
        rationale="per-condition random-effects calibration",
    )


def node(parameter: str, i: int, j: int, coordinates=None) -> CalibrationFieldNode:
    calibration = calibrate_hot_fire_scale(
        parameter, stratum(parameter, i, j), source="fixture", rationale="fixture"
    )
    return CalibrationFieldNode(coordinates or (PRESSURES[i], RATIOS[j]), calibration)


def point(pc: float, mr: float) -> dict[str, float]:
    return {"chamber_pressure_pa": pc, "mixture_ratio": mr}


class CalibrationFieldInterpolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.field = make_field()

    def test_nodes_are_reproduced_exactly(self) -> None:
        for pc in PRESSURES:
            for mr in RATIOS:
                calibration = self.field.node_calibration((pc, mr))
                result = self.field.evaluate(point(pc, mr))
                self.assertEqual(result.consensus_scale, calibration.consensus_scale)
                self.assertEqual(
                    result.consensus_standard_uncertainty,
                    calibration.consensus_standard_uncertainty,
                )
                self.assertEqual(
                    result.independent_nodes_standard_uncertainty,
                    calibration.consensus_standard_uncertainty,
                )
                self.assertEqual(
                    result.between_engine_process_variance,
                    calibration.between_engine_process_variance,
                )
                self.assertEqual(
                    result.within_engine_process_variance,
                    calibration.within_engine_process_variance,
                )
                self.assertEqual(
                    result.population_process_variance, calibration.population_process_variance
                )
                self.assertEqual(result.node_weights, (((pc, mr), 1.0),))
                self.assertEqual(set(result.boundary_axes), {"chamber_pressure_pa", "mixture_ratio"})
                self.assertGreater(calibration.between_engine_process_variance, 0.0)

    def test_edge_midpoint_and_cell_centre_follow_the_stated_rule(self) -> None:
        a = self.field.node_calibration((1_000_000.0, 5.5))
        b = self.field.node_calibration((2_000_000.0, 5.5))
        mid = self.field.evaluate(point(1_500_000.0, 5.5))
        self.assertAlmostEqual(mid.consensus_scale, (a.consensus_scale + b.consensus_scale) / 2, 15)
        # JCGM 100 5.2.2 NOTE 1 (r = +1): linear sum with sensitivities 1/2.
        self.assertAlmostEqual(
            mid.consensus_standard_uncertainty,
            (a.consensus_standard_uncertainty + b.consensus_standard_uncertainty) / 2,
            15,
        )
        # JCGM 100 equation 10 (independent nodes).
        self.assertAlmostEqual(
            mid.independent_nodes_standard_uncertainty,
            math.hypot(a.consensus_standard_uncertainty, b.consensus_standard_uncertainty) / 2,
            15,
        )
        self.assertAlmostEqual(
            mid.between_engine_process_variance,
            (a.between_engine_process_variance + b.between_engine_process_variance) / 2,
            15,
        )
        corners = [
            self.field.node_calibration((pc, mr))
            for pc in (2_000_000.0, 3_000_000.0)
            for mr in RATIOS
        ]
        centre = self.field.evaluate(point(2_500_000.0, 5.75))
        self.assertAlmostEqual(
            centre.consensus_scale, sum(item.consensus_scale for item in corners) / 4, 15
        )
        self.assertAlmostEqual(
            centre.population_process_variance,
            sum(item.population_process_variance for item in corners) / 4,
            15,
        )
        self.assertEqual(len(centre.node_weights), 4)
        self.assertAlmostEqual(sum(weight for _, weight in centre.node_weights), 1.0, 15)
        self.assertEqual(centre.cell_lower_coordinates, (2_000_000.0, 5.5))
        self.assertEqual(centre.cell_upper_coordinates, (3_000_000.0, 6.0))
        # Off-centre weights: t = 0.25 along pressure, t = 0.2 along mixture ratio.
        weights = dict(self.field.evaluate(point(1_250_000.0, 5.6)).node_weights)
        expected = {
            (1_000_000.0, 5.5): 0.75 * 0.8,
            (1_000_000.0, 6.0): 0.75 * 0.2,
            (2_000_000.0, 5.5): 0.25 * 0.8,
            (2_000_000.0, 6.0): 0.25 * 0.2,
        }
        self.assertEqual(set(weights), set(expected))
        for key, value in expected.items():
            # 5.6 - 5.5 is not exactly 0.1 in binary64.
            self.assertAlmostEqual(weights[key], value, delta=1e-14)
        off = self.field.evaluate(point(1_250_000.0, 5.6))
        self.assertAlmostEqual(
            off.consensus_standard_uncertainty,
            sum(
                value * self.field.node_calibration(key).consensus_standard_uncertainty
                for key, value in expected.items()
            ),
            15,
        )

    def test_bounds_ordering_identity_and_monotonicity(self) -> None:
        previous = None
        for step in range(41):
            pc = PRESSURES[0] + step * (PRESSURES[-1] - PRESSURES[0]) / 40
            for mr in (5.5, 5.6, 5.83, 6.0):
                result = self.field.evaluate(point(pc, mr))
                corners = [
                    self.field.node_calibration((x, y))
                    for x in {result.cell_lower_coordinates[0], result.cell_upper_coordinates[0]}
                    for y in {result.cell_lower_coordinates[1], result.cell_upper_coordinates[1]}
                ]
                mus = [item.consensus_scale for item in corners]
                us = [item.consensus_standard_uncertainty for item in corners]
                tolerance = 1e-15
                self.assertLessEqual(min(mus) - tolerance, result.consensus_scale)
                self.assertLessEqual(result.consensus_scale, max(mus) + tolerance)
                self.assertLessEqual(min(us) - tolerance, result.consensus_standard_uncertainty)
                self.assertLessEqual(result.consensus_standard_uncertainty, max(us) + tolerance)
                self.assertLessEqual(
                    result.independent_nodes_standard_uncertainty,
                    result.consensus_standard_uncertainty + tolerance,
                )
                self.assertAlmostEqual(
                    result.population_process_variance,
                    result.within_engine_process_variance + result.between_engine_process_variance,
                    17,
                )
                self.assertGreaterEqual(result.between_engine_process_variance, 0.0)
            # Node consensus increases with pressure, so the field is monotone in it.
            value = self.field.evaluate(point(pc, 5.83)).consensus_scale
            if previous is not None:
                self.assertGreaterEqual(value, previous - 1e-15)
            previous = value

    def test_gradient_matches_central_difference(self) -> None:
        conditions = point(1_700_000.0, 5.7)
        result = self.field.evaluate(conditions)
        gradient = dict(result.consensus_scale_gradient)
        for name, step in (("chamber_pressure_pa", 100.0), ("mixture_ratio", 1e-4)):
            up = dict(conditions)
            down = dict(conditions)
            up[name] += step
            down[name] -= step
            estimate = (
                self.field.evaluate(up).consensus_scale - self.field.evaluate(down).consensus_scale
            ) / (2 * step)
            self.assertAlmostEqual(gradient[name], estimate, delta=abs(estimate) * 1e-6 + 1e-15)

    def test_extrapolation_is_refused_and_domain_events_are_exposed(self) -> None:
        for bad in (point(3_000_100.0, 5.8), point(999_000.0, 5.8), point(2e6, 6.01)):
            with self.assertRaises(DomainError):
                self.field.evaluate(bad)
        with self.assertRaises(DomainError):
            self.field.evaluate({"chamber_pressure_pa": 2e6})
        with self.assertRaises(DomainError):
            self.field.evaluate({**point(2e6, 5.8), "throttle": 1.0})
        # Within the snapping tolerance of an endpoint is accepted and snapped.
        snapped = self.field.evaluate(point(3_000_000.0 + 1e-7, 6.0))
        self.assertEqual(dict(snapped.conditions)["chamber_pressure_pa"], 3_000_000.0)
        events = self.field.validity_surfaces()
        self.assertEqual(len(events), 4)
        for event in events:
            self.assertTrue(event.stops_burn)
            self.assertTrue(event.event_id.startswith("calibration-field:thrust_scale:"))
            self.assertGreater(event.value(point(2e6, 5.8)), 0.0)
        self.assertLess(min(event.value(point(3.2e6, 5.8)) for event in events), 0.0)
        self.assertEqual(
            min(event.value(point(3e6, 5.8)) for event in events), 0.0
        )
        self.assertEqual(
            [item.coordinate for item in self.field.grid_line_surfaces()], [2_000_000.0]
        )


class CalibrationFieldContractTests(unittest.TestCase):
    def nodes(self) -> list[CalibrationFieldNode]:
        return [node("thrust_scale", i, j) for i in range(3) for j in range(2)]

    def build(self, nodes) -> HotFireCalibrationField:
        return HotFireCalibrationField("thrust_scale", AXES, tuple(nodes), "fixture", "fixture")

    def test_complete_field_is_accepted(self) -> None:
        self.assertEqual(len(self.build(self.nodes()).nodes), 6)

    def test_missing_node_is_refused(self) -> None:
        with self.assertRaisesRegex(DomainError, "missing=1"):
            self.build(self.nodes()[:-1])

    def test_duplicate_node_is_refused(self) -> None:
        nodes = self.nodes()
        duplicate = CalibrationFieldNode(
            nodes[0].coordinates,
            calibrate_hot_fire_scale(
                "thrust_scale",
                stratum("thrust_scale", 0, 0, condition="other"),
                source="fixture",
                rationale="fixture",
            ),
        )
        with self.assertRaisesRegex(DomainError, "Duplicate"):
            self.build([*nodes, duplicate])

    def test_off_grid_and_outside_grid_nodes_are_refused(self) -> None:
        nodes = self.nodes()
        with self.assertRaisesRegex(DomainError, "not on a grid coordinate"):
            self.build([*nodes[:-1], node("thrust_scale", 2, 1, (2_500_000.0, 6.0))])
        with self.assertRaisesRegex(DomainError, "outside the grid"):
            self.build([*nodes[:-1], node("thrust_scale", 2, 1, (4_000_000.0, 6.0))])
        with self.assertRaises(DomainError):
            self.build([*nodes[:-1], node("thrust_scale", 2, 1, (3_000_000.0,))])

    def test_condition_reuse_and_parameter_mismatch_are_refused(self) -> None:
        nodes = self.nodes()
        reused = CalibrationFieldNode(nodes[-1].coordinates, nodes[0].calibration)
        with self.assertRaisesRegex(DomainError, "two nodes"):
            self.build([*nodes[:-1], reused])
        with self.assertRaises(DomainError):
            self.build([*nodes[:-1], node("specific_impulse_scale", 2, 1)])

    def test_measurement_grouping_fails_closed(self) -> None:
        measurements = all_measurements()
        coordinates = mapping()
        kwargs = {"source": "fixture", "rationale": "fixture"}
        unmapped = dict(coordinates)
        del unmapped[condition_id(1, 1)]
        with self.assertRaisesRegex(DomainError, "unmapped"):
            calibrate_hot_fire_field("thrust_scale", AXES, unmapped, measurements, **kwargs)
        shared = dict(coordinates)
        shared[condition_id(1, 1)] = shared[condition_id(0, 0)]
        with self.assertRaisesRegex(DomainError, "Duplicate"):
            calibrate_hot_fire_field("thrust_scale", AXES, shared, measurements, **kwargs)
        empty = {key: value for key, value in coordinates.items() if key != condition_id(2, 1)}
        empty["unused"] = (PRESSURES[2], RATIOS[1])
        without = tuple(item for item in measurements if item.condition_id != condition_id(2, 1))
        with self.assertRaisesRegex(DomainError, "no measurements"):
            calibrate_hot_fire_field("thrust_scale", AXES, empty, without, **kwargs)
        outside = dict(coordinates)
        outside[condition_id(2, 1)] = (3_500_000.0, 6.0)
        with self.assertRaisesRegex(DomainError, "outside the grid"):
            calibrate_hot_fire_field("thrust_scale", AXES, outside, measurements, **kwargs)
        missing = tuple(item for item in measurements if item.condition_id != condition_id(2, 1))
        reduced = {key: value for key, value in coordinates.items() if key != condition_id(2, 1)}
        with self.assertRaisesRegex(DomainError, "missing=1"):
            calibrate_hot_fire_field("thrust_scale", AXES, reduced, missing, **kwargs)


class CalibrationFieldPersistenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.field = make_field()

    def test_round_trip_reproduces_hash_and_evaluations(self) -> None:
        text = self.field.to_json()
        reopened = HotFireCalibrationField.from_json(text)
        self.assertEqual(reopened.field_hash, self.field.field_hash)
        self.assertEqual(reopened.to_json(), text)
        conditions = point(1_234_567.0, 5.66)
        self.assertEqual(reopened.evaluate(conditions), self.field.evaluate(conditions))

    def test_tampering_is_refused(self) -> None:
        payload = json.loads(self.field.to_json())
        payload["nodes"][0]["calibration"]["measurements"][0]["value"] += 0.01
        with self.assertRaises(InputError):
            HotFireCalibrationField.from_dict(payload)

        forged = json.loads(self.field.to_json())
        forged["nodes"][0]["calibration"]["consensus_scale"] += 0.01
        del forged["field_hash"]
        forged["field_hash"] = canonical_sha256(forged)
        with self.assertRaises(InputError):
            HotFireCalibrationField.from_dict(forged)

        policy = json.loads(self.field.to_json())
        policy["uncertainty_numerics"] = "variance_average"
        del policy["field_hash"]
        policy["field_hash"] = canonical_sha256(policy)
        with self.assertRaises(InputError):
            HotFireCalibrationField.from_dict(policy)

        schema = json.loads(self.field.to_json())
        schema["schema"] = "other"
        with self.assertRaises(InputError):
            HotFireCalibrationField.from_dict(schema)
        with self.assertRaises(InputError):
            HotFireCalibrationField.from_json("[]")


class CalibrationFieldCompositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.thrust = make_field("thrust_scale")
        cls.isp = make_field("specific_impulse_scale")
        cls.surface = RectilinearPerformanceSurface(
            axes=AXES,
            points=tuple(
                PerformanceScalePoint((pc, mr), 0.8 + 0.1 * i + 0.01 * j, 0.9 + 0.05 * i - 0.02 * j)
                for i, pc in enumerate(PRESSURES)
                for j, mr in enumerate(RATIOS)
            ),
            source_id="map-A",
            source_sha256="e" * 64,
        )

    def test_corrected_nodes_and_uncertainties(self) -> None:
        calibrated = apply_calibration_fields(
            self.surface, thrust_field=self.thrust, specific_impulse_field=self.isp
        )
        base = {item.coordinates: item for item in self.surface.points}
        for item in calibrated.surface.points:
            mu_f = self.thrust.node_calibration(item.coordinates).consensus_scale
            mu_i = self.isp.node_calibration(item.coordinates).consensus_scale
            self.assertAlmostEqual(item.force_scale, base[item.coordinates].force_scale * mu_f, 15)
            self.assertAlmostEqual(
                item.mass_flow_scale, base[item.coordinates].mass_flow_scale * mu_f / mu_i, 15
            )
        self.assertEqual(calibrated.base_surface_hash, self.surface.surface_hash)
        conditions = point(2_250_000.0, 5.9)
        result = calibrated.evaluate(conditions)
        thrust = self.thrust.evaluate(conditions)
        isp = self.isp.evaluate(conditions)
        a = thrust.consensus_standard_uncertainty / thrust.consensus_scale
        b = isp.consensus_standard_uncertainty / isp.consensus_scale
        self.assertEqual(result.force_scale_relative_standard_uncertainty, a)
        self.assertEqual(result.specific_impulse_relative_standard_uncertainty, b)
        self.assertEqual(result.mass_flow_scale_relative_standard_uncertainty, a + b)
        reopened = RectilinearPerformanceSurface.from_json(calibrated.surface.to_json())
        self.assertEqual(reopened.surface_hash, calibrated.surface.surface_hash)
        self.assertEqual(calibrated.validity_surfaces(), calibrated.surface.validity_surfaces())
        with self.assertRaises(DomainError):
            calibrated.evaluate(point(3.5e6, 5.9))

    def test_thrust_only_correction_keeps_isp(self) -> None:
        calibrated = apply_calibration_fields(self.surface, thrust_field=self.thrust)
        for item, base in zip(
            sorted(calibrated.surface.points, key=lambda p: p.coordinates),
            sorted(self.surface.points, key=lambda p: p.coordinates),
            strict=True,
        ):
            self.assertAlmostEqual(
                item.force_scale / item.mass_flow_scale, base.force_scale / base.mass_flow_scale, 14
            )

    def test_composition_refuses_mismatched_inputs(self) -> None:
        with self.assertRaises(DomainError):
            apply_calibration_fields(self.surface)
        with self.assertRaises(DomainError):
            apply_calibration_fields(self.surface, thrust_field=self.isp)
        other = RectilinearPerformanceSurface(
            axes=(
                PerformanceAxis("chamber_pressure_pa", "Pa", PRESSURES),
                PerformanceAxis("mixture_ratio", "1", (5.5, 6.5)),
            ),
            points=tuple(
                PerformanceScalePoint((pc, mr), 1.0, 1.0) for pc in PRESSURES for mr in (5.5, 6.5)
            ),
            source_id="map-B",
            source_sha256="f" * 64,
        )
        with self.assertRaises(DomainError):
            apply_calibration_fields(other, thrust_field=self.thrust)

    def test_ensemble_parameters_are_scenario_ready(self) -> None:
        calibrated = apply_calibration_fields(
            self.surface, thrust_field=self.thrust, specific_impulse_field=self.isp
        )
        conditions = point(1_600_000.0, 5.5)
        parameters = calibrated.ensemble_parameters(conditions)
        self.assertEqual(
            [item.name for item in parameters],
            [
                "thrust_scale::aleatory:hot-fire-field",
                "thrust_scale::epistemic:hot-fire-field-consensus",
                "specific_impulse_scale::aleatory:hot-fire-field",
                "specific_impulse_scale::epistemic:hot-fire-field-consensus",
            ],
        )
        self.assertEqual(
            [item.uncertainty_class for item in parameters],
            [UncertaintyClass.ALEATORY, UncertaintyClass.EPISTEMIC] * 2,
        )
        definition = ScenarioConditionalUncertainty(
            scenario_id="field-point",
            parameters=parameters,
            sample_count=16,
            source="calibration field",
            rationale="interpolated per-condition calibration",
        )
        self.assertEqual(len(definition.parameters), 4)
        evaluation = self.thrust.evaluate(conditions)
        residual = parameters[0]
        self.assertAlmostEqual(
            residual.standard_deviation,
            math.sqrt(evaluation.population_process_variance) / evaluation.consensus_scale,
            15,
        )

    def test_residual_parameter_at_node_matches_node_calibration(self) -> None:
        coordinates = (2_000_000.0, 6.0)
        calibration = self.thrust.node_calibration(coordinates)
        for scope in CalibrationScope:
            expected = calibration.residual_uncertain_parameter(scope=scope)
            actual = self.thrust.residual_uncertain_parameter(point(*coordinates), scope=scope)
            self.assertEqual(actual.name, "thrust_scale")
            self.assertEqual(actual.standard_deviation, expected.standard_deviation)
            self.assertEqual((actual.lower, actual.upper), (expected.lower, expected.upper))
        with self.assertRaises(DomainError):
            self.thrust.residual_uncertain_parameter(point(*coordinates), coverage_factor=0.0)


if __name__ == "__main__":
    unittest.main()
