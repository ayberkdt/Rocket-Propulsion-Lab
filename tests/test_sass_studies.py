"""Independent checks for the sourced UZM 508E preparation package.

Oracles: NASA Glenn ideal rocket equation / Isp; MIT 16.522 L5-L6 circular
orbital energy; NASA ideal-gas law; NASA SE Handbook 6.8 and Appendix D.
The mission residual oracle uses bisection on independently propagated mass,
not the production closed-form wet-mass expression.
"""

from __future__ import annotations

import ast
import copy
import json
import math
import re
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2 as G0
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256
from rocket_propulsion.studies import sass
from rocket_propulsion.studies.sass import (
    Citation,
    CitedScalar,
    Criterion,
    DeltaVItem,
    DesignCase,
    Maneuver,
    Requirement,
    assumption,
    bipropellant_split,
    build_project_html,
    burn_duration_screen,
    check_requirements,
    circular_spiral_delta_v,
    delta_v_budget,
    design_sensitivity,
    hohmann_transfer,
    ideal_cold_gas_storage,
    ideal_regulated_pressurant,
    liquid_tank_volume,
    nominal_operating_point,
    perigee_lowering_delta_v,
    power_resources,
    project_from_dict,
    project_from_json,
    size_mission,
    weighted_trade,
)
from rocket_propulsion.studies.sass import catalog as catalog_module
from rocket_propulsion.studies.sass.__main__ import main
from rocket_propulsion.studies.sass.evidence import derived

ROOT = Path(__file__).resolve().parents[1]


def q(value: float, unit: str = "1") -> CitedScalar:
    """Explicit synthetic test quantity; NASA-SE-6 6.8 assumption discipline."""
    return assumption(value, unit, "Synthetic oracle input; not catalogue or course data.")


def maneuver(name: str, dv: float, isp: float = 200, thrust: float = 10,
             axial: float = 1, duty: float = 1) -> Maneuver:
    """Construct an ideal test burn; NASA-ISP and NASA-ROCKET identities."""
    return Maneuver(name, q(dv, "m/s"), q(thrust, "N"), q(isp, "s"), q(axial), q(duty))


def example_request() -> dict:
    """Read the explicit synthetic scenario; NASA-SE-APP reproducible inputs."""
    return json.loads((ROOT / "examples" / "sass_mission.json").read_text(encoding="utf-8"))


def catalogue_study():
    """Run the normal sizing entry point with a fixed in-memory test case."""
    return sass.analyze_project(
        "Catalogue review fixture", q(50, "kg"), (DeltaVItem("burn", q(5, "m/s"), q(0), 1),),
        (sass.Candidate("nitrogen", "moog-058e151", q(4, "kg"), 4, 2, q(1), q(1)),),
        residual_fraction=q(0), fixed_reserve=q(0, "kg"), requirements=(),
        criteria=(Criterion("wet_mass", q(55, "kg"), q(100, "kg"), q(1)),),
    )


def coherent_life_catalogue():
    """Isolate life screening with a synthetic envelope around the Ariane 1 N study point."""
    record = sass.thruster("ariane-1n")
    flow = record.field("reported_mass_flow_range")
    upper = replace(flow.upper, value=0.0005,
                    citations=(Citation("ARIANE-1N", "Synthetic upper flow endpoint for life-test isolation"),))
    fields = tuple(replace(field, upper=upper) if field is flow else field for field in record.fields)
    catalogue = tuple(replace(item, fields=fields) if item is record else item
                      for item in catalog_module.THRUSTERS)
    return patch.object(catalog_module, "THRUSTERS", catalogue)


def life_study(on_time, *, installed=1, firing=1, count=1, duty=1.0, axial=1.0, reserve=0.0):
    """Choose delta-v from the independently specified firing time and retained mass."""
    exhaust_speed = 220.0 * 9.80665
    dry_mass = 500.0 + installed * 0.29
    propellant = firing * on_time / exhaust_speed
    delta_v = axial * exhaust_speed * math.log1p(propellant / (dry_mass + reserve))
    with coherent_life_catalogue():
        return sass.analyze_project(
            "Life screening fixture", q(500, "kg"), (DeltaVItem("burn", q(delta_v / count, "m/s"), q(0), count),),
            (sass.Candidate("hydrazine", "ariane-1n", q(0, "kg"), installed, firing, q(axial), q(duty)),),
            residual_fraction=q(0), fixed_reserve=q(reserve, "kg"), requirements=(),
            criteria=(Criterion("wet_mass", q(500, "kg"), q(3000, "kg"), q(1)),),
        )


class EvidenceTests(unittest.TestCase):
    def test_no_uncited_value_unknown_source_or_implicit_assumption(self) -> None:
        with self.assertRaises(DomainError):
            CitedScalar(1, "N", ())
        with self.assertRaises(DomainError):
            Citation("invented", "page 1")
        with self.assertRaises(DomainError):
            Citation("NASA-ISP", "")
        with self.assertRaises(DomainError):
            CitedScalar(1, "N", (Citation("NASA-ISP", "Isp equation"),), "assumption")
        for value in (True, None, "1", float("nan"), float("inf"), 10**400):
            with self.subTest(value=str(value)[:25]), self.assertRaises(DomainError):
                q(value)

    def test_derived_keeps_original_and_method_sources(self) -> None:
        published = CitedScalar(65, "s", (Citation("MOOG-COLD", "058E151 Isp"),))
        result = derived(1, "N", (published,), (Citation("NASA-ISP", "F=mdot*g0*Isp"),))
        self.assertEqual({x.source_id for x in result.citations}, {"MOOG-COLD", "NASA-ISP"})
        self.assertEqual(result.status, "derived")
        with self.assertRaises(DomainError):
            derived(1, "N", (published,), ())

    def test_every_shipped_function_and_class_has_registered_references(self) -> None:
        missing = []
        keys = {ref.key for ref in sass.REFERENCES}
        for path in Path(sass.__file__).parent.glob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    doc = ast.get_docstring(node) or ""
                    if "References\n" not in doc or not any(key in doc for key in keys):
                        missing.append(f"{path.name}:{node.lineno}:{node.name}")
        self.assertEqual(missing, [])

    def test_catalog_citations_and_spec_semantics_are_complete(self) -> None:
        self.assertEqual(len({x.key for x in sass.REFERENCES}), len(sass.REFERENCES))
        for item in sass.list_thrusters():
            self.assertTrue(item.citations)
            self.assertTrue(item.mass_scope)
            for field in item.fields:
                self.assertTrue(field.upper.citations)
                if field.lower is not None:
                    self.assertTrue(field.lower.citations)
        moog = sass.thruster("moog-058e151")
        self.assertEqual(moog.field("thrust").upper.value, 0.12)
        self.assertEqual(moog.field("isp").upper.value, 65)
        self.assertEqual(moog.field("dry_mass").upper.value, 0.07)
        self.assertEqual(moog.field("opening_power").kind, "upper_bound")

    def test_source_registry_is_present_in_the_reference_index(self) -> None:
        index = (ROOT / "docs" / "REFERENCES.md").read_text(encoding="utf-8")
        for source in sass.REFERENCES:
            self.assertIn(source.url, index, source.key)
            self.assertIn(f"**{source.key}**", index)
            self.assertTrue(source.locator.strip())
            self.assertTrue(source.url.startswith("https://"))

    def test_bounds_and_source_disagreement_do_not_create_operating_points(self) -> None:
        for key in ("ecaps-1n-hpgp", "ecaps-1n-hpgp-bradford", "busek-bit3-page"):
            with self.assertRaisesRegex(DomainError, "nominal"):
                nominal_operating_point(key)
        bit = sass.thruster("busek-bit3-page")
        self.assertEqual({x.source_id for x in bit.citations}, {"BUSEK-BIT3", "BUSEK-RF"})
        self.assertEqual(bit.field("dry_mass").upper.citations[0].source_id, "BUSEK-BIT3-PDF")
        with self.assertRaises(DomainError):
            sass.thruster("moog-058e151").field("tank_mass")

    def test_manufacturer_flow_disagreement_is_exposed_not_reconciled_silently(self) -> None:
        point = nominal_operating_point("ariane-1n")
        self.assertEqual(point.get("reported_mass_flow").value, 0.00044)
        self.assertAlmostEqual(point.get("closed_mass_flow").value, 1 / (220 * G0))
        self.assertGreater(abs(point.get("reported_flow_relative_difference").value), 0.01)

    def test_ecaps_sources_preserve_distinct_isp_and_preheat_power_scope(self) -> None:
        survey = sass.thruster("ecaps-1n-hpgp")
        supplier = sass.thruster("ecaps-1n-hpgp-bradford")
        self.assertEqual(survey.field("isp").upper.value, 235)
        self.assertEqual(supplier.field("isp").upper.value, 231)
        self.assertEqual(survey.field("survey_power").upper.value, 10)
        self.assertIn("reactor_preheat_power", {field.name for field in supplier.fields})
        preheat = supplier.field("reactor_preheat_power")
        self.assertEqual((preheat.lower.value, preheat.upper.value), (8, 10))
        self.assertEqual(preheat.upper.citations[0].source_id, "ECAPS-1N-BRADFORD")
        for record in (survey, supplier):
            with self.assertRaises(DomainError):
                record.field("power")

    def test_nominal_labels_do_not_pair_different_operating_points_or_sources(self) -> None:
        moog = sass.thruster("moog-058e151")
        isp = moog.field("isp")
        foreign = replace(isp.upper, citations=(Citation("ARIANE-1N", "different thruster"),))
        for changed in (replace(isp, nominal_point_id=None),
                        replace(isp, nominal_point_id="different-inlet-condition"),
                        replace(isp, lower=foreign, upper=foreign)):
            altered = replace(moog, fields=tuple(changed if field.name == "isp" else field
                                                for field in moog.fields))
            catalogue = (altered, *catalog_module.THRUSTERS[1:])
            with self.subTest(isp=changed), patch.object(catalog_module, "THRUSTERS", catalogue):
                with self.assertRaisesRegex(DomainError, "paired nominal"):
                    nominal_operating_point(moog.key)
                study = catalogue_study()
                self.assertIsNone(study.candidates[0].sizing)
                self.assertIsNone(study.scores[0].score)
        with self.assertRaisesRegex(DomainError, "Only nominal fields"):
            replace(moog.field("opening_power"), nominal_point_id=isp.nominal_point_id)

    def test_study_point_cannot_publish_a_flow_range_endpoint_as_nominal(self) -> None:
        original = sass.thruster("ariane-1n")
        flow = original.field("reported_mass_flow")
        interval = replace(flow, kind="range", lower=replace(flow.lower, value=0.00042))
        altered = replace(original, fields=tuple(interval if field.name == flow.name else field
                                                 for field in original.fields))
        with patch.object(catalog_module, "THRUSTERS", (altered,)):
            result = nominal_operating_point(original.key)
            self.assertAlmostEqual(result.get("closed_mass_flow").value, 1 / (220 * G0))
            with self.assertRaises(DomainError):
                result.get("reported_mass_flow")

    def test_nominal_flow_preserves_finite_result_across_intermediate_overflow(self) -> None:
        original = sass.thruster("moog-058e151")
        fields = []
        for field in original.fields:
            if field.name in ("thrust", "isp"):
                value = replace(field.upper, value=5e307 if field.name == "thrust" else 1e308)
                field = replace(field, lower=value, upper=value)
            fields.append(field)
        altered = replace(original, fields=tuple(fields))
        with patch.object(catalog_module, "THRUSTERS", (altered,)):
            result = nominal_operating_point(original.key)
        self.assertAlmostEqual(result.get("closed_mass_flow").value, 0.5 / 9.80665, places=15)


class ManeuverTests(unittest.TestCase):
    def test_hohmann_matches_independent_dimensioned_circular_example(self) -> None:
        result = hohmann_transfer(q(1, "m"), q(4, "m"), q(1, "m^3/s^2"))
        self.assertAlmostEqual(result.get("departure_delta_v").value, math.sqrt(1.6) - 1)
        self.assertAlmostEqual(result.get("arrival_delta_v").value, 0.5 - math.sqrt(0.1))
        self.assertAlmostEqual(result.get("transfer_coast_time").value, math.pi * 2.5**1.5)
        reverse = hohmann_transfer(q(4, "m"), q(1, "m"), q(1, "m^3/s^2"))
        self.assertAlmostEqual(result.get("total_delta_v").value, reverse.get("total_delta_v").value)
        self.assertLess(reverse.get("departure_delta_v").value, 0)
        self.assertLess(reverse.get("arrival_delta_v").value, 0)

    def test_equal_radius_zero_and_low_thrust_distinction(self) -> None:
        result = hohmann_transfer(q(1, "m"), q(1, "m"), q(1, "m^3/s^2"))
        self.assertEqual(result.get("total_delta_v").value, 0)
        self.assertEqual(result.get("transfer_coast_time").value, 0)
        spiral = circular_spiral_delta_v(q(1, "m"), q(4, "m"), q(1, "m^3/s^2"))
        self.assertEqual(spiral.get("total_delta_v").value, 0.5)
        impulse = hohmann_transfer(q(1, "m"), q(4, "m"), q(1, "m^3/s^2"))
        self.assertGreater(spiral.get("total_delta_v").value, impulse.get("total_delta_v").value)

    def test_perigee_lowering_uses_only_first_impulse(self) -> None:
        lowered = perigee_lowering_delta_v(q(4, "m"), q(1, "m"), q(1, "m^3/s^2"))
        transfer = hohmann_transfer(q(4, "m"), q(1, "m"), q(1, "m^3/s^2"))
        self.assertEqual(lowered.get("total_delta_v").value,
                         abs(transfer.get("departure_delta_v").value))
        self.assertLess(lowered.get("total_delta_v").value, transfer.get("total_delta_v").value)
        with self.assertRaises(DomainError):
            perigee_lowering_delta_v(q(1, "m"), q(4, "m"), q(1, "m^3/s^2"))

    def test_units_domains_and_analyst_threshold(self) -> None:
        with self.assertRaises(DomainError):
            hohmann_transfer(q(1, "km"), q(4, "m"), q(1, "m^3/s^2"))
        with self.assertRaises(DomainError):
            hohmann_transfer(q(0, "m"), q(4, "m"), q(1, "m^3/s^2"))
        result = burn_duration_screen(q(3, "s"), q(100, "s"), q(0.01))
        self.assertLess(result.get("fraction_margin").value, 0)


class SizingTests(unittest.TestCase):
    def test_delta_v_budget_applies_counts_and_margins_exactly_once(self) -> None:
        result = delta_v_budget((DeltaVItem("a", q(10, "m/s"), q(0.2), 3),
                                 DeltaVItem("b", q(4, "m/s"), q(0), 1)))
        self.assertEqual(result.get("nominal_delta_v").value, 34)
        self.assertEqual(result.get("allocated_delta_v").value, 40)
        for count in (0, -1, True, 1.5):
            with self.assertRaises(DomainError):
                DeltaVItem("a", q(10, "m/s"), q(0), count)

    def test_variable_isp_residual_closure_against_independent_bisection(self) -> None:
        burns = (maneuver("first", 300, 200, 10, 0.8, 0.4),
                 maneuver("second", 90, 310, 2, 0.95, 0.5))
        result = size_mission(q(80, "kg"), burns, residual_fraction=q(0.03),
                              fixed_reserve=q(0.5, "kg"))
        low, high = 0.0, 100.0
        for _ in range(100):
            loaded = 0.5 * (low + high)
            mass = 80 + loaded
            for burn in burns:
                mass *= math.exp(-burn.delta_v.value / (G0 * burn.system_isp.value
                                                        * burn.axial_efficiency.value))
            if mass > 80 + 0.03 * loaded + 0.5:
                high = loaded
            else:
                low = loaded
        self.assertAlmostEqual(result.summary.get("loaded_propellant").value, loaded, places=10)
        self.assertAlmostEqual(result.summary.get("mass_closure_residual").value, 0, places=11)
        for burn, (_, outputs) in zip(burns, result.burns, strict=True):
            m0, mf = outputs.get("initial_mass").value, outputs.get("final_mass").value
            self.assertAlmostEqual(G0 * burn.system_isp.value * burn.axial_efficiency.value
                                   * math.log(m0 / mf), burn.delta_v.value, places=10)
            used = outputs.get("consumed_propellant").value
            self.assertAlmostEqual(outputs.get("tank_mass_flow").value
                                   * outputs.get("on_time").value, used)
            self.assertAlmostEqual(outputs.get("scalar_impulse").value,
                                   used * G0 * burn.system_isp.value)

    def test_small_delta_v_and_zero_burn_do_not_lose_mass_precision(self) -> None:
        result = size_mission(q(100, "kg"), (maneuver("micro", 1e-9),),
                              residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        expected = 100 * math.expm1(1e-9 / (G0 * 200))
        self.assertAlmostEqual(result.summary.get("loaded_propellant").value / expected, 1)
        zero = size_mission(q(100, "kg"), (maneuver("zero", 0),),
                            residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        self.assertEqual(zero.summary.get("on_time").value, 0)

    def test_cant_changes_useful_performance_without_hiding_drain(self) -> None:
        aligned = size_mission(q(100, "kg"), (maneuver("a", 100),),
                               residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        canted = size_mission(q(100, "kg"), (maneuver("a", 100, axial=0.5),),
                              residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        self.assertGreater(canted.summary.get("loaded_propellant").value,
                           aligned.summary.get("loaded_propellant").value)
        self.assertEqual(canted.burns[0][1].get("tank_mass_flow").value,
                         aligned.burns[0][1].get("tank_mass_flow").value)

    def test_infeasible_residual_units_and_operating_domains_fail(self) -> None:
        with self.assertRaises(DomainError):
            size_mission(q(100, "kg"), (maneuver("a", G0 * 200),),
                         residual_fraction=q(0.5), fixed_reserve=q(0, "kg"))
        with self.assertRaises(DomainError):
            maneuver("a", 10, duty=0)
        with self.assertRaises(DomainError):
            maneuver("a", 10, axial=1.1)

    def test_extreme_exhaust_velocity_and_flow_fail_with_domain_errors(self) -> None:
        for isp, thrust in ((1e308, 1), (5e-324, 1), (200, 5e-324)):
            with self.subTest(isp=isp, thrust=thrust), self.assertRaises(DomainError):
                maneuver("numeric-domain", 1, isp=isp, thrust=thrust)

    def test_split_and_resource_ledger(self) -> None:
        split = bipropellant_split(q(10, "kg"), q(3))
        self.assertEqual(split.get("fuel_mass").value, 2.5)
        self.assertEqual(split.get("oxidizer_mass").value, 7.5)
        resources = power_resources(q(20, "W"), q(10, "s"), q(0.5), q(2, "W"),
                                     q(50, "W"), q(4, "s"))
        self.assertEqual(resources.get("energy").value, 420)
        self.assertEqual(resources.get("elapsed_time").value, 24)
        self.assertEqual(resources.get("peak_bus_power").value, 50)
        self.assertEqual(resources.get("average_bus_power").value, 17.5)


class StorageTests(unittest.TestCase):
    def test_ullage_is_fraction_of_tank_not_multiplier_on_liquid(self) -> None:
        result = liquid_tank_volume(q(8, "kg"), q(1000, "kg/m^3"), q(0.2))
        self.assertAlmostEqual(result.get("tank_volume").value, 0.01)
        self.assertAlmostEqual(result.get("ullage_volume").value, 0.002)
        with self.assertRaises(DomainError):
            liquid_tank_volume(q(8, "kg"), q(1000, "kg/m^3"), q(1))

    def test_cold_gas_matches_endpoint_inventories(self) -> None:
        result = ideal_cold_gas_storage(q(2, "kg"), q(3000000, "Pa"), q(1000000, "Pa"),
                                        q(300, "K"), q(300, "J/(kg*K)"))
        volume = result.get("tank_volume").value
        self.assertAlmostEqual(volume, 0.09)
        self.assertAlmostEqual(3000000 * volume / (300 * 300), 3)
        self.assertAlmostEqual(result.get("remaining_gas").value, 1)
        with self.assertRaises(DomainError):
            ideal_cold_gas_storage(q(2, "kg"), q(1, "Pa"), q(2, "Pa"),
                                   q(300, "K"), q(300, "J/(kg*K)"))

    def test_regulated_pressurant_has_independent_mass_closure(self) -> None:
        result = ideal_regulated_pressurant(q(0.01, "m^3"), q(2000000, "Pa"),
                                            q(20000000, "Pa"), q(3000000, "Pa"),
                                            q(300, "K"), q(2000, "J/(kg*K)"))
        self.assertAlmostEqual(result.get("loaded_pressurant").value
                               - result.get("remaining_pressurant").value,
                               result.get("delivered_pressurant").value)
        self.assertAlmostEqual(result.get("delivered_pressurant").value, 20000 / 600000)


class TradeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.criteria = (Criterion("mass", q(1, "kg"), q(3, "kg"), q(1)),)
        self.requirements = (Requirement("M", "Mass ceiling", "mass", "at_most", q(2, "kg")),)

    def test_missing_metrics_are_unknown_and_failed_candidates_are_not_ranked(self) -> None:
        rows = weighted_trade({"good": {"mass": q(1, "kg")}, "bad": {"mass": q(4, "kg")},
                               "unknown": {}}, self.criteria, self.requirements)
        self.assertEqual([r.eligibility for r in rows], ["eligible", "infeasible", "incomplete"])
        self.assertEqual(rows[0].score.value, 1)
        self.assertIsNone(rows[1].score)
        self.assertIsNone(rows[2].score)
        self.assertEqual(rows[2].checks[0].status, "unknown")

    def test_fixed_anchors_prevent_candidate_dependent_rescaling(self) -> None:
        before = weighted_trade({"a": {"mass": q(1.5, "kg")}}, self.criteria, ())
        after = weighted_trade({"a": {"mass": q(1.5, "kg")}, "dominated": {"mass": q(100, "kg")}},
                                self.criteria, ())
        self.assertEqual(before[0].score.value, after[0].score.value)
        self.assertEqual(after[1].score.value, 0)
        with self.assertRaises(DomainError):
            weighted_trade({"a": {}}, (replace(self.criteria[0], weight=q(0.5)),), ())

    def test_constraint_units_and_duplicate_ids_are_rejected(self) -> None:
        with self.assertRaises(DomainError):
            check_requirements({"mass": q(1, "m")}, self.requirements)
        with self.assertRaises(DomainError):
            check_requirements({}, self.requirements * 2)


class ProjectTests(unittest.TestCase):
    def test_example_reproduces_and_preserves_unknown_architectures(self) -> None:
        first = project_from_dict(example_request())
        second = project_from_dict(example_request())
        self.assertEqual(first.to_dict(), second.to_dict())
        self.assertEqual(project_from_dict(first.to_dict()["inputs"]).to_dict(), first.to_dict())
        self.assertEqual(len(first.candidates), 5)
        self.assertEqual(first.scores[3].eligibility, "incomplete")
        self.assertEqual(first.scores[4].eligibility, "incomplete")
        self.assertIsNone(first.candidates[4].sizing)
        payload = first.to_dict()
        expected_hash = payload.pop("study_hash")
        self.assertEqual(expected_hash, canonical_sha256(payload))

    def test_dry_mass_counts_installed_heads_once_not_only_firing_heads(self) -> None:
        study = project_from_dict(example_request())
        first = study.candidates[0]
        self.assertAlmostEqual(first.sizing.summary.get("dry_mass").value, 50 + 4 + 4 * 0.07)
        self.assertAlmostEqual(first.sizing.burns[0][1].get("tank_mass_flow").value,
                               2 * 0.12 / (65 * G0))

    def test_input_copy_cannot_rewrite_the_calculation_provenance(self) -> None:
        study = project_from_dict(example_request())
        before = study.to_dict()
        inputs = study.input_data
        inputs["items"][0]["nominal_delta_v"]["value"] = 999
        self.assertEqual(study.to_dict(), before)

    def test_old_study_cannot_export_a_newer_catalogue_record(self) -> None:
        study = catalogue_study()
        original_document = study.to_dict()
        moog = sass.thruster("moog-058e151")
        mass = moog.field("dry_mass")
        updated_mass = replace(mass.upper, value=0.17)
        newer_mass = replace(mass, lower=updated_mass, upper=updated_mass)
        newer_record = replace(moog, fields=tuple(newer_mass if field.name == "dry_mass" else field
                                                   for field in moog.fields))
        with patch.object(catalog_module, "THRUSTERS", (newer_record, *catalog_module.THRUSTERS[1:])):
            self.assertEqual(study.to_dict(), original_document)
            rerun = catalogue_study()
            self.assertAlmostEqual(rerun.candidates[0].sizing.summary.get("dry_mass").value,
                                   50 + 4 + 4 * 0.17)
            self.assertNotEqual(rerun.to_dict()["study_hash"], original_document["study_hash"])

    def test_catalogue_mass_with_wrong_units_or_sign_cannot_enter_bus_budget(self) -> None:
        moog = sass.thruster("moog-058e151")
        mass = moog.field("dry_mass")
        for invalid in (replace(mass.upper, unit="g"), replace(mass.upper, value=-0.07)):
            bad_mass = replace(mass, lower=invalid, upper=invalid)
            bad_record = replace(moog, fields=tuple(bad_mass if field.name == "dry_mass" else field
                                                     for field in moog.fields))
            with self.subTest(mass=invalid), patch.object(
                catalog_module, "THRUSTERS", (bad_record, *catalog_module.THRUSTERS[1:]),
            ), self.assertRaises(DomainError):
                catalogue_study()

    def test_continuous_burn_limit_is_not_applied_to_aggregated_or_rotating_heads(self) -> None:
        continuous = life_study(45_000, installed=2, firing=2)
        check = {row.metric: row for row in continuous.candidates[0].catalogue_life_checks}
        self.assertAlmostEqual(check["single_burn_life"].demand.value, 45_000)
        self.assertEqual(check["single_burn_life"].published_limit.value, 43_200)
        self.assertEqual(check["single_burn_life"].status, "exceeds_reference")
        self.assertEqual(continuous.scores[0].eligibility, "incomplete")
        self.assertIsNone(continuous.scores[0].score)
        self.assertIn("exceeds reference", build_project_html(continuous))
        for options in (
            {"installed": 2, "firing": 2, "count": 2},
            {"installed": 2, "firing": 2, "duty": 0.4},
            {"installed": 3, "firing": 2},
        ):
            with self.subTest(options=options):
                grouped = life_study(45_000, **options)
                rows = {row.metric: row for row in grouped.candidates[0].catalogue_life_checks}
                self.assertIsNone(rows["single_burn_life"].demand)
                self.assertEqual(rows["single_burn_life"].status, "unassessed")
                self.assertEqual(rows["cycle_life"].published_limit.value, 59_000)
                self.assertIsNone(rows["cycle_life"].demand)
                self.assertEqual(rows["cycle_life"].status, "unassessed")
                self.assertIsNotNone(grouped.scores[0].score)

    def test_accumulated_burn_reference_uses_simultaneous_and_installed_counts(self) -> None:
        for firing, expected_demand, expected_status in (
            (1, 95_000, "no_exceedance_identified"), (2, 190_000, "exceeds_reference"),
        ):
            with self.subTest(firing=firing):
                study = life_study(190_000, installed=2, firing=firing, duty=0.25)
                checks = {row.metric: row for row in study.candidates[0].catalogue_life_checks}
                accumulated = checks["accumulated_burn_life"]
                self.assertAlmostEqual(accumulated.demand.value, expected_demand)
                self.assertEqual(accumulated.published_limit.value, 180_000)
                self.assertEqual(accumulated.status, expected_status)
                self.assertEqual(study.scores[0].score is None, firing == 2)
                self.assertAlmostEqual(checks["propellant_throughput"].demand.value,
                                       expected_demand / (9.80665 * 220.0))

    def test_throughput_uses_consumed_propellant_and_impulse_uses_scalar_load(self) -> None:
        # The reserve dwarfs the 67 kg reference but does not pass through the head.
        # Axial efficiency makes delivered axial impulse 14,000 N s, while the
        # head still delivers 140,000 N s of scalar impulse.
        study = life_study(140_000, duty=0.4, axial=0.1, reserve=1000.0)
        rows = {row.metric: row for row in study.candidates[0].catalogue_life_checks}
        self.assertAlmostEqual(rows["propellant_throughput"].demand.value, 140_000 / (9.80665 * 220))
        self.assertEqual(rows["propellant_throughput"].published_limit.value, 67)
        self.assertEqual(rows["propellant_throughput"].status, "no_exceedance_identified")
        self.assertAlmostEqual(rows["total_impulse_life"].demand.value, 140_000)
        self.assertEqual(rows["total_impulse_life"].published_limit.value, 135_000)
        self.assertEqual(rows["total_impulse_life"].status, "exceeds_reference")
        self.assertIsNone(study.scores[0].score)
        s10 = sass.thruster("ariane-s10-13")
        self.assertEqual(s10.field("single_burn_life").upper.value, 28_800)
        self.assertEqual(s10.field("accumulated_burn_life").upper.value, 248_400)
        with self.assertRaises(DomainError):
            s10.field("propellant_throughput")

    def test_report_keeps_citations_and_escapes_hostile_titles_and_notes(self) -> None:
        request = example_request()
        request["title"] = "<script>alert(1)</script>"
        html = build_project_html(project_from_dict(request))
        self.assertNotIn("<script>alert(1)</script>", html)
        self.assertIn("&lt;script&gt;", html)
        self.assertIn("incomplete", html)
        self.assertIn("BUSEK-RF", html)
        embedded = re.search(r'<script type="application/json" id="study-document">(.*?)</script>', html)
        self.assertIsNotNone(embedded)
        self.assertEqual(json.loads(embedded.group(1))["title"], request["title"])

    def test_unknown_fields_bare_numbers_and_extrapolation_are_refused(self) -> None:
        for transform in ("unknown", "bare", "extrapolated", "uncited"):
            raw = example_request()
            if transform == "unknown":
                raw["extrapolate"] = True
            elif transform == "bare":
                raw["residual_fraction"] = 0.1
            elif transform == "extrapolated":
                raw["residual_fraction"]["status"] = "extrapolated"
            else:
                raw["residual_fraction"]["citations"] = []
            with self.subTest(transform=transform), self.assertRaises(DomainError):
                project_from_dict(raw)

    def test_hostile_json_fails_with_typed_errors(self) -> None:
        for raw in ("null", "[]", "{", "NaN", "9" * 5000, "[" * 10000,
                    None, b"{}", 10):
            with self.subTest(value=str(raw)[:20]), self.assertRaises((InputError, DomainError)):
                project_from_json(raw)
        raw = example_request()
        raw["residual_fraction"]["value"] = float("nan")
        with self.assertRaises(DomainError):
            project_from_json(json.dumps(raw))

    def test_sensitivity_reruns_from_baseline_and_never_changes_catalogue(self) -> None:
        raw = example_request()
        original = copy.deepcopy(raw)
        cases = (DesignCase("baseline", q(1), q(1)), DesignCase("more-dv", q(2), q(1)))
        first, second = design_sensitivity(raw, cases)
        self.assertEqual(raw, original)
        self.assertGreater(second.candidates[0].sizing.summary.get("loaded_propellant").value,
                           first.candidates[0].sizing.summary.get("loaded_propellant").value)
        self.assertEqual(second.candidates[0].operating_point, first.candidates[0].operating_point)

    def test_cli_matches_python_and_refuses_input_overwrite(self) -> None:
        source = ROOT / "examples" / "sass_mission.json"
        parent = ROOT / 'output' / 'assurance' / 'test_artifacts' / 'sass_cli'
        parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as temp:
            out, report = Path(temp) / "study.json", Path(temp) / "study.html"
            self.assertEqual(main([str(source), "--output", str(out), "--report", str(report)]), 0)
            self.assertEqual(json.loads(out.read_text(encoding="utf-8")),
                             project_from_dict(example_request()).to_dict())
            self.assertIn("MOOG-COLD", report.read_text(encoding="utf-8"))
            self.assertEqual(main([str(source), "--output", str(source), "--report", str(report)]), 2)


if __name__ == "__main__":
    unittest.main()
