"""Create a source-bound inventory and an explicitly unachieved hardware evidence ledger.

References: NASA-7009B data pedigree, verification/validation domains and reporting;
NASA ESTO hardware TRL5. This is a gap ledger, not an automated readiness assessor.
"""

import ast
import hashlib
import json
import platform
import sys
from dataclasses import asdict
from pathlib import Path

from rocket_propulsion.studies.sass import REFERENCES, list_thrusters, nominal_operating_point

ROOT = Path(__file__).resolve().parents[1]
FEED_REVIEW = {
    "feed_system.py": {"IdealPressurantGas.__post_init__", "IdealPressurantGas.specific_heat_cv_j_kg_k",
        "IdealPressurantGas.specific_heat_cp_j_kg_k", "RegulatedFeedSystem._regulator_flow",
        "RegulatedFeedSystem._ullage_volume", "RegulatedFeedSystem._pressures",
        "RegulatedFeedSystem.initial_state", "RegulatedFeedSystem.evaluate",
        "RegulatedFeedDerivative.__post_init__", "FeedPressurePartials.__post_init__",
        "RegulatedFeedEvaluation.__post_init__"},
    "feed_line.py": {"DynamicFeedLine.evaluate", "DynamicFeedLinePartials.__post_init__",
                     "DynamicFeedLineDerivative.__post_init__",
                     "DynamicFeedLineEvaluation.__post_init__",
                     "evaluate_dynamic_feed_propulsion",
                     "DynamicFeedPropulsionEvaluation.reverse_flow_event_extension"},
    "bipropellant.py": {"BipropellantDynamicFeedCoupling.evaluate",
                       "BipropellantDynamicFeedEvaluation.reverse_flow_event_extension",
                       "BipropellantDynamicFeedConfiguration.__post_init__",
                       "BipropellantDynamicFeedState.__post_init__",
                       "BipropellantMixtureRatioEventSurface.__post_init__",
                       "BipropellantMixtureRatioEventSurface.value",
                       "BipropellantDynamicFeedCoupling.__post_init__",
                       "_finite_quotient", "_mixture_ratio_partial",
                       "_finite_columns", "_finite_pairs",
                       "BipropellantCoupledPartials.__post_init__",
                       "BipropellantDynamicFeedEvaluation.__post_init__"},
    "numeric_outputs.py": {"_fraction", "_finite_rational"},
}
BRIDGE_FLOW_REVIEW = {"_mass_profile", "_mass_profile.channel",
                     "PropagatorPropulsionBridge._evaluate",
                     "PropagatorPropulsionBridge._evaluate.state_parameter_partials",
                     "PropagatorPropulsionBridge._evaluate.state_condition_partials"}
PULSE_REVIEW = {"PulseModeThruster.plan_delta_v" + suffix for suffix in (
    "", ".delta_v", ".widest", ".low_end", ".high_end")}
PULSE_REVIEW.add("PulseModeThruster._solve_budget_width")
API_REVIEW = {
    "api/capabilities.py": {"_failed_probe", "SideraCapabilityProbe", "SideraCapabilityProbe.__init__",
                            "SideraCapabilityProbe._run", "SideraCapabilityProbe.snapshot"},
    "api/routes.py": {"sidera_capabilities_route"},
}
DIRECTION_REVIEW = {
    "propagator.py": {"PropagatorPropulsionBridge._evaluate"},
    "direction.py": {"_DirectionGeometry", "_DirectionGeometry.project", "_DirectionGeometry.jacobian",
                     "_direction_geometry"},
    "throttleable.py": {"ThrottleableEngineModel.evaluate"},
    "electric.py": {"PowerLimitedPropulsionModel.evaluate"},
}
THROTTLE_FLOW_REVIEW = {"ThrottleableEngineModel.__post_init__", "ThrottleableEngineModel.evaluate",
                        "ThrottleableEngineModel.mass_closure"}
NUMERICAL_POLICY_REVIEW = {
    "serialization.py": {"_require_numerical_policies"},
    "propagator.py": {"PropagatorPropulsionBridge.from_dict"},
    "electric.py": {"PowerLimitedPropulsionModel.from_dict"},
    "throttleable.py": {"ThrottleableEngineModel.from_dict"},
}
CUBIC_REVIEW = {
    "numeric_outputs.py": {"_finite_ratio"},
    "interpolation.py": {"pchip_slopes", "_pchip_rational_slopes",
                         "_pchip_rational_slopes.end_slope", "_hermite_coefficients", "_hermite_segment"},
    "throttleable.py": {"MonotoneCubicCurve.__post_init__", "MonotoneCubicCurve.evaluate",
                        "MonotoneCubicCurve._evaluate_checked", "MonotoneCubicCurve._value_checked",
                        "ThrottleableEngineModel.manifest"},
    "pulse_mode.py": {"PulseModeThruster._raw", "PulseModeThruster._raw_values", "PulseModeThruster.manifest",
                      "PulseModeThruster.from_dict"},
    "performance_surface.py": {"RectilinearPerformanceSurface._hermite_derivatives",
                               "RectilinearPerformanceSurface._evaluate_hermite",
                               "RectilinearPerformanceSurface.evaluate",
                               "RectilinearPerformanceSurface.to_dict",
                               "RectilinearPerformanceSurface.from_dict"},
}
MULTILINEAR_REVIEW = {
    "interpolation.py": {"_multilinear_terms"},
    "performance_surface.py": {"RectilinearPerformanceSurface.evaluate",
                               "RectilinearPerformanceSurface.to_dict",
                               "RectilinearPerformanceSurface.from_dict"},
    "stream_performance.py": {"RectilinearStreamPerformanceSurface.evaluate",
                              "RectilinearStreamPerformanceSurface.to_dict",
                              "RectilinearStreamPerformanceSurface.from_dict"},
}
INVENTORY_CLOSURE_REVIEW = {
    "mass_closure.py": {"_inventory_ledger", "_closure_ratio"},
    "propagator.py": {"PropellantMassClosure.__post_init__", "PropellantMassClosure.is_closed",
                      "PropagatorPropulsionBridge.mass_closure"},
    "electric.py": {"PowerLimitedPropulsionModel.mass_closure"},
    "throttleable.py": {"ThrottleableEngineModel.mass_closure"},
}
CONSTANT_BURN_REVIEW = {"_negative_exponential", "_mass_log", "_target_duration", "simulate_constant_burn"}
PROFILE_BURN_REVIEW = {
    "models.py": {"ThrottleSegment.__post_init__", "ThrottleSegment.exposure_s",
                  "ThrottleSchedule.__post_init__", "ThrottleSchedule.total_duration_s", "ThrottleSchedule.total_exposure_s"},
    "profiles.py": {"ScheduleInterval.__post_init__", "ScheduleInterval.duration_s", "ScheduleInterval.exposure_s",
                    "ScheduleInterval.throttle_first_moment_s2", "_ScheduleSlice", "_ScheduleSlice.exposure",
                    "_ScheduleSlice.moment", "_schedule_totals", "_schedule_slices", "_schedule_integrals",
                    "_sqrt_fraction", "_inverse_exposure", "resolved_schedule_intervals", "integrated_throttle",
                    "throttle_first_moment", "time_for_integrated_throttle", "ProfiledPerformanceProvider.discontinuities_s"},
    "profile_integration.py": {"_required_exposure", "_tank_rates", "_tank_masses_at_exposure", "_residual_upper",
                               "_point", "_build_profile", "_profile_cutoff", "simulate_profiled_burn"},
    "integration.py": {"_negative_exponential", "_mass_log", "_target_duration"},
}
ADAPTER_DIRECTION_REVIEW = {
    "integrations/direction.py": {"_directed_vector", "_stored_unit_direction"},
    "integrations/sidera.py": {"_normalized_direction", "_validated_export_inputs", "export_sidera_artifact",
                               "to_native_finite_burn", "to_native_maneuver_plan"},
    "integrations/sidera_ensemble.py": {"_unit_direction", "SideraEnsembleBurnJob.__post_init__",
                                       "SideraEnsembleHandoff.to_dict", "SideraEnsembleHandoff.from_dict",
                                       "export_sidera_ensemble_handoff"},
    "integrations/lunaris.py": {"lunaris_impulsive_equivalents", "LunarisImpulsiveEquivalent.to_lunaris"},
    "propulsion/burns/direction.py": {"_direction_geometry", "_DirectionGeometry.project"},
}


RECORD_CONTRACT_CLASSES = {
    'PropellantFlow', 'EngineOperatingPoint', 'TankDefinition', 'BurnTarget', 'BurnDefinition',
    'PropulsionState', 'BurnProfilePoint', 'ThrottleSegment', 'ThrottleSchedule',
    'ProfiledBurnPoint', 'BurnSummary', 'IntegrationEvidence', 'BurnResult', 'ProfiledBurnResult',
}
RECORD_CONTRACT_REVIEW = {
    'models.py': RECORD_CONTRACT_CLASSES | {name+'.__post_init__' for name in RECORD_CONTRACT_CLASSES}
        | {'PropellantRole', 'FlowDestination', 'BurnTargetKind', 'CutoffReason',
           '_finite_nonnegative', '_finite_positive', '_boolean', '_enum', '_utf8',
           '_model_fields', '_model_strings', '_inventory_pairs', '_point_contract', '_result_contract'},
    'serialization.py': {'burn_result_to_dict', 'burn_result_to_json', '_record_data', '_array',
                         'burn_result_from_dict', 'burn_result_from_dict.record',
                         '_unique_object', '_reject_constant', 'burn_result_from_json'},
}

FREE_SPACE_REVIEW = {
    'gauss_legendre', '_gauss_legendre_rule', '_gauss_legendre_rule.legendre', '_decimal',
    '_mass_intervals', '_mass_intervals.mass_at', '_quadrature', '_bernstein_enclosure',
    '_mass_log_enclosure', '_upper_output', 'FreeSpaceBurnMetrics',
    '_reference_integrals', '_exact_integrals', '_constant_exhaust',
    '_range_output', '_range_output.rounded',
    'FreeSpaceBurnMetrics.__post_init__', 'free_space_burn_metrics',
}
TIME_COORDINATE_REVIEW = {
    'propulsion/burns/time_coordinates.py': {'_epoch_table','_epoch_profile_time'},
    'propulsion/burns/tabulated.py': {'TabulatedBurnArtifact.__post_init__',
        'TabulatedBurnArtifact.boundary_tolerance_s','TabulatedBurnArtifact.segment_index_at',
        'TabulatedBurnArtifact._segment_index_at','TabulatedBurnArtifact._snap_to_boundary',
        'TabulatedBurnArtifact.snap_to_boundary','TabulatedBurnArtifact.is_boundary_time',
        'TabulatedBurnArtifact._is_boundary_time','TabulatedBurnSegment._fraction_at',
        'TabulatedBurnArtifact._history_dict','TabulatedBurnArtifact.from_dict'},
    'propulsion/burns/propagator.py': {'PropagatorPropulsionBridge.evaluate',
        'PropagatorPropulsionBridge.evaluate_at_epoch','PropagatorPropulsionBridge._evaluate',
        'PropagatorPropulsionBridge.time_event_jump_at_epoch','PropagatorPropulsionBridge._jump_from_evaluations',
        'PropagatorPropulsionBridge.time_event_jump',
        'PropagatorPropulsionBridge._segment_at','PropagatorPropulsionBridge.manifest'},
    'propulsion/burns/timeline.py': {'ScheduledBurn.cutoff_time_s', 'PropulsionTimeline.__post_init__',
        'PropulsionTimeline.events','PropulsionTimeline.evaluate','PropulsionTimeline.manifest',
        'PropulsionTimeline.from_dict'},
    'integrations/lunaris.py': {'lunaris_burn_windows','lunaris_segmentation_times','_merged_times',
        'LunarisFiniteBurnDynamics.__init__','LunarisFiniteBurnDynamics.predicted_mass_kg',
        'LunarisFiniteBurnDynamics._stopped','LunarisFiniteBurnDynamics._wrap',
        'LunarisFiniteBurnDynamics._wrap.rhs'},
}
TABULATED_LEDGER_REVIEW = {
    '_segment_ledger', '_artifact_ledger', 'TabulatedBurnSegment.__post_init__',
    'TabulatedBurnSegment._validate_stream_closure', 'TabulatedBurnSegment._integral',
    'TabulatedBurnSegment._local_first_moment', 'TabulatedBurnSegment.delivered_impulse_n_s',
    'TabulatedBurnSegment.axial_impulse_n_s', 'TabulatedBurnSegment.propellant_mass_kg',
    'TabulatedBurnSegment.delivered_thrust_first_moment_n_s2', 'TabulatedBurnSegment._fraction_at',
    'TabulatedBurnSegment._axial_at', 'TabulatedBurnSegment.state_at',
    'TabulatedBurnSegment.state_at.interpolate', 'TabulatedBurnArtifact.__post_init__',
    'TabulatedBurnArtifact.delivered_total_impulse_n_s', 'TabulatedBurnArtifact.axial_total_impulse_n_s',
    'TabulatedBurnArtifact.consumed_propellant_kg', 'TabulatedBurnArtifact.thrust_centroid_time_s',
    'TabulatedBurnArtifact.system_equivalent_specific_impulse_s', 'TabulatedBurnArtifact.artifact_hash',
    'TabulatedBurnArtifact._history_dict', 'TabulatedBurnArtifact.to_dict', 'TabulatedBurnArtifact.to_json',
    'TabulatedBurnArtifact.from_dict', 'scale_tabulated_artifact', 'scale_tabulated_artifact.scaled',
    'scale_tabulated_artifact.scaled_streams', 'tabulated_artifact_from_profiled_burn',
    'tabulated_artifact_from_profiled_burn.flows', 'tabulated_artifact_from_csv',
}


def symbols(nodes, prefix=""):
    """Visit explicit declarations/lambdas through all AST children and control flow."""
    for node in nodes:
        child_prefix = prefix
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = prefix + node.name
            yield name, node
            child_prefix = name + '.'
        elif isinstance(node, ast.Lambda):
            name = prefix + f'<lambda@{node.lineno}:{node.col_offset}>'
            yield name, node
            child_prefix = name + '.'
        yield from symbols(ast.iter_child_nodes(node), child_prefix)


def _review_scopes(path, name):
    """Retain overlapping documented scopes instead of hiding earlier reviews."""
    scopes = (["source_trace_and_software_verification_only"]
              if "studies/sass/" in path.as_posix() else [])
    groups = (
        (TIME_COORDINATE_REVIEW.get(path.relative_to(ROOT / 'src/rocket_propulsion').as_posix(),set()),
         'targeted_time_coordinate_ownership_and_epoch_event_mapping_only'),
        (FREE_SPACE_REVIEW if path.name == 'burn_metrics.py' else set(),
         'targeted_free_space_quadrature_and_directed_numerical_enclosure_only'),
        (TABULATED_LEDGER_REVIEW if path.name == 'tabulated.py' else
         {'PropagatorPropulsionBridge._evaluate', 'PropagatorPropulsionBridge.manifest'}
         if path.name == 'propagator.py' else set(),
         'targeted_tabulated_ledgers_producers_and_numeric_mode_boundaries_only'),
        (RECORD_CONTRACT_REVIEW.get(path.name, set()),
         'targeted_public_burn_record_and_artifact_boundary_verification_only'),
        ({"_role", "_checked_streams", "_checked_tanks", "_group_rates", "_checked_inventory",
          "flow_by_tank", "validate_tank_coverage", "tank_duration_limit", "tank_masses_at"}
         if path.name == "inventory.py" else set(),
         "targeted_direct_constant_inventory_helper_verification_only"),
        (FEED_REVIEW.get(path.name, set()), "targeted_feed_physics_and_software_verification_only"),
        (INVENTORY_CLOSURE_REVIEW.get(path.name, set()), "targeted_inventory_residual_and_decision_verification_only"),
        (CONSTANT_BURN_REVIEW if path.name == "integration.py" else set(), "targeted_constant_burn_numerical_verification_only"),
        (PROFILE_BURN_REVIEW.get(path.name, set()), "targeted_profile_burn_numerical_verification_only"),
        (ADAPTER_DIRECTION_REVIEW.get(path.relative_to(ROOT / "src/rocket_propulsion").as_posix(), set()),
         "targeted_adapter_direction_and_native_boundary_verification_only"),
        (MULTILINEAR_REVIEW.get(path.name, set()), "targeted_multilinear_map_arithmetic_verification_only"),
        (CUBIC_REVIEW.get(path.name, set()), "targeted_pchip_and_cubic_interpolation_verification_only"),
        (THROTTLE_FLOW_REVIEW if path.name == "throttleable.py" else set(),
         "targeted_throttle_flow_ledger_and_inventory_verification_only"),
        (NUMERICAL_POLICY_REVIEW.get(path.name, set()), "targeted_numerical_policy_document_verification_only"),
        (BRIDGE_FLOW_REVIEW if path.name == "propagator.py" else set(),
         "targeted_bridge_mass_flow_and_sensitivity_verification_only"),
        (DIRECTION_REVIEW.get(path.name, set()), "targeted_direction_geometry_and_projection_verification_only"),
        (PULSE_REVIEW if path.name == "pulse_mode.py" else set(), "targeted_pulse_search_and_software_verification_only"),
        (API_REVIEW.get(path.relative_to(ROOT / "src/rocket_propulsion").as_posix(), set()),
         "targeted_optional_dependency_http_verification_only"),
    )
    scopes.extend(label for names, label in groups if name in names)
    return scopes


def main():
    """Bind the audit inventory to exact source bytes and report unclosed gates."""
    inventory = []
    for path in sorted((ROOT / "src/rocket_propulsion").rglob("*.py")):
        data = path.read_bytes()
        tree = ast.parse(data.decode("utf-8"))
        for name, node in symbols(tree.body):
            doc = '' if isinstance(node, ast.Lambda) else ast.get_docstring(node) or ''
            scopes = _review_scopes(path, name)
            inventory.append({"file": path.relative_to(ROOT).as_posix(), "line": node.lineno,
                              "symbol": name, "kind": type(node).__name__, "column": node.col_offset,
                              "source_sha256": hashlib.sha256(data).hexdigest(),
                              "registered_sass_references": [r.key for r in REFERENCES
                                                            if r.key in doc],
                              "review_status": scopes[0] if scopes else "not_reaudited_in_this_change",
                              "review_scopes": scopes,
                              "empirical_validation": "not_demonstrated"})
    catalog = []
    for item in list_thrusters():
        row = {"record": asdict(item), "joint_measurement_uncertainty": "missing",
               "motor_and_feed_validation": "missing", "qualification_use": "prohibited"}
        if item.field("thrust").kind == item.field("isp").kind == "nominal":
            point = nominal_operating_point(item.key)
            values = dict(point.values)
            row["flow_relative_discrepancy"] = (
                values["reported_flow_relative_difference"].value
                if "reported_flow_relative_difference" in values else None)
        catalog.append(row)
    dossier = {
        "scope": "motor and feed system; design and mission requirements undefined (user-confirmed)",
        "hardware_readiness": "not_assessed", "trl5_achieved": False,
        "qualification_evidence": False,
        "runtime": {"python": sys.version, "platform": platform.platform()},
        "inventory_summary": {"symbols_total": len(inventory),
                              "scope": "Explicit Python AST class/function/async-function/lambda declarations under src/rocket_propulsion, including control-flow bodies and signatures. Lambda identities use static containment plus line/column; generated/runtime callables, scripts outside this source tree and other languages are not enumerated.",
                              "sass_symbols": sum("studies/sass/" in x["file"] for x in inventory)},
        "inventory": inventory, "catalog": catalog,
        "gates": [
            {"id": identifier, "status": "unmet", "needed": needed}
            for identifier, needed in (
                ("CONFIGURATION", "selected engine, propellants, feed schematic, controlled BOM and interfaces"),
                ("RELEVANT_ENVIRONMENT", "mission-derived pressures, thermal, vacuum, duty and coupled loads"),
                ("ACCEPTANCE", "approved requirements, critical performance, limits and uncertainty policy"),
                ("TEST_ARTICLE", "identified brassboard and realistic support hardware with configuration records"),
                ("MEASUREMENT", "calibrated sensors, DAQ timing and traceable measurement uncertainty budget"),
                ("DEMONSTRATION", "integrated relevant-environment operation with archived measured results"),
                ("MODEL_VALIDATION", "independent holdout comparison within an explicitly bounded validation domain"),
                ("SCALING", "documented analytical agreement and scaling requirements for subsequent phases"),
                ("REVIEW", "independent technical review, discrepancy closure and readiness assessment"),
            )],
        "references": [asdict(r) for r in REFERENCES if r.key in {"NASA-7009B", "NASA-TRL"}],
    }
    native_manifest = ROOT / "docs/assurance/lunaris_impulse_patch.json"
    if native_manifest.exists():
        native = json.loads(native_manifest.read_text())
        live_digest = hashlib.sha256(Path(native["live_path"]).read_bytes()).hexdigest()
        dossier["external_native_reviews"] = [{
            "file": native["live_path"],
            "symbols": ["_impulse_norm_and_mass", "apply_impulsive_maneuver"],
            "review_scope": "targeted_native_impulse_arithmetic_and_output_contract_only",
            "live_source_sha256": live_digest,
            "reviewed_candidate_sha256": native["candidate_source_sha256"],
            "live_source_matches_reviewed_candidate": live_digest == native["candidate_source_sha256"],
            "empirical_validation": "not_demonstrated",
            "package_wide_correctness": "not_established",
        }]
    output = ROOT / "output/assurance/audit.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dossier, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(dossier["inventory_summary"]))
    print("Hardware readiness: not assessed; TRL5 evidence gates unmet: 9")


if __name__ == "__main__":
    main()
