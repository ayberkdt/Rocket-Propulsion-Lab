"""Reintroduce specific semantic defects in memory and challenge the regression tests.

No repository source is overwritten. A killed mutant must cause an assertion
failure; an unrelated exception is an invalid experiment, not a successful kill.
These curated mutations are not an exhaustive mutation-coverage claim.
References: NASA-7009B, implementation verification and documented limitations.
"""

import hashlib
import importlib
import json
import sys
from pathlib import Path

from scripts.run_bounded_regression import ROOT, verification_inputs

MUTATIONS = (
    ("perturb_final_mass", "mass_numerics", "_at_precision", '("final_mass", mass, "kg")',
     '("final_mass", mass * Decimal("1.00000000005"), "kg")',
     "test_large_mass_ratio_keeps_positive_dry_mass", "adversarial"),
    ("erase_oxidizer", "sizing", "bipropellant_split", "total * (ratio / (1 + ratio))",
     "total - fuel", "test_small_oxidizer_inventory_does_not_disappear", "adversarial"),
    ("erase_remaining_gas", "storage", "ideal_cold_gas_storage", "remaining = m * pf / (p0 - pf)",
     "remaining = Fraction(float(loaded) - float(m))",
     "test_small_gas_residual_does_not_disappear", "adversarial"),
    ("subtract_adjacent_speeds", "maneuvers", "hohmann_transfer",
     "dv0 = v0 * x / (sqrt(r1 / a) + 1)", "dv0 = v0 * sqrt(r1 / a) - v0",
     "test_neighbouring_orbits_preserve_both_small_signed_impulses", "adversarial"),
    ("trust_catalogue_as_validation", "catalog", "nominal_operating_point",
     "if require_validated:", "if False and require_validated:",
     "test_no_catalog_entry_can_claim_validated_joint_performance", "adversarial"),
    ("discard_uncertainty", "validation", "compare_validation_sample",
     "bound = Fraction(sample.measurement_bound.value) + Fraction(sample.prediction_bound.value)",
     "bound = Fraction(0)",
     "test_uncertainty_cannot_make_a_bad_prediction_pass", "validation"),
    ("ignore_domain_review", "validation", "compare_validation_sample",
     "if not sample.in_validation_domain:", "if False and not sample.in_validation_domain:",
     "test_missing_bounds_and_outside_domain_cannot_pass", "validation"),
    ("ignore_failed_holdout", "validation", "assess_validation",
     '"fail" if "fail" in states else', '"pass" if "fail" in states else',
     "test_one_failed_sample_dominates_unknown_and_passed_samples", "validation"),
    ("subtract_rounded_phase_times", "power_numerics", "solve_power_budget",
     "standby_time = t * (1 - d) / d",
     "standby_time = Fraction(float(t / d) - float(t))",
     "test_near_continuous_duty_preserves_small_standby_energy", "power_boundaries"),
    ("round_phase_energies_before_sum", "power_numerics", "solve_power_budget",
     "energy = p * t + idle * standby_time + warm * tw",
     "energy = Fraction(float(p * t)) + Fraction(float(idle * standby_time)) + Fraction(float(warm * tw))",
     "test_phase_energy_survives_underflow_of_individual_float_products", "power_boundaries"),
    ("round_regulated_density_before_displacement", "storage", "ideal_regulated_pressurant",
     "used = pt * displaced / thermal", "used = Fraction(float(pt / thermal)) * displaced",
     "test_regulated_transfer_survives_underflow_of_density_before_displacement", "storage_boundaries"),
    ("allow_zero_for_positive_final_channel", "rational_numerics", "binary_ledger",
     "if not isfinite(binary) or (value != 0 and binary == 0):", "if not isfinite(binary):",
     "test_positive_unrepresentable_ullage_cannot_be_reported_as_empty", "storage_boundaries"),
    ("round_weighted_terms_before_score", "trade", "weighted_trade",
     "float(sum(terms, Fraction(0)))", "fsum(float(x) for x in terms)",
     "test_equal_affine_scores_cannot_gain_a_false_order_from_term_rounding", "decision_boundaries"),
    ("erase_positive_scaled_requirement", "sensitivity", "_scaled_quantity",
     'binary = binary_ledger((("scaled_quantity", scaled, value.unit),), "Design sensitivity")[0][1]',
     "binary = float(scaled)",
     "test_positive_scaled_delta_v_cannot_silently_become_zero", "decision_boundaries"),
    ("round_interval_upper_inward", "validation", "_outward_endpoint",
     "if upward and Fraction(binary) < value:", "if False and upward and Fraction(binary) < value:",
     "test_sampled_bounds_and_decisions_are_sound_against_exact_interval_corners", "validation"),
    ("discard_numerical_review_guard", "validation", "compare_validation_sample",
     "guard = _outward_endpoint(guard_value, upward=True)", "guard = 0.0",
     "test_declared_numerical_review_guard_is_not_silently_removed", "validation"),
)


def main():
    """Record assertion kills, survivors and invalid mutation experiments separately."""
    fingerprints = verification_inputs(ROOT)
    outcomes = []
    for name, module_name, symbol, before, after, test_name, suite_name in MUTATIONS:
        module = importlib.import_module("rocket_propulsion.studies.sass." + module_name)
        function = getattr(module, symbol)
        from scripts.mutation_runtime import MUTATION_RUNTIME_POLICY, compile_fault, evaluate_fault
        source,replacement = compile_fault(function,before,after,f"<mutation:{name}>")
        tests = importlib.import_module("tests.test_sass_" + suite_name)
        cls = getattr(tests, {"adversarial": "HostileRegressionTests", "validation": "ValidationAttacks",
                              "power_boundaries": "PowerBoundaryAttacks",
                              "storage_boundaries": "StorageBoundaryAttacks",
                              "decision_boundaries": "DecisionBoundaryAttacks"}[suite_name])
        # Patch every imported binding of the exact original callable. Patching
        # only its defining module misses helpers imported into a consumer's
        # globals, creating a false survivor without exercising the mutant.
        def aliases(*, function=function):return [(loaded, key) for name, loaded in list(sys.modules.items())
                   if name.startswith(("rocket_propulsion.studies.sass", "tests.test_sass_"))
                   and loaded is not None for key, value in list(vars(loaded).items())
                   if value is function]
        trial=evaluate_fault(aliases,function,replacement,cls,test_name)
        outcomes.append({"mutation": name, **trial,
                         "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                         "test": cls.__module__ + "." + cls.__name__ + "." + test_name})
    output = Path("output/assurance/mutations.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    if verification_inputs(ROOT) != fingerprints:
        raise RuntimeError("Mutation inputs changed during the run.")
    output.write_text(json.dumps({"scope": "sixteen curated semantic defects; not exhaustive",
                                  "source_sha256": fingerprints,
                                  "runtime_policy": MUTATION_RUNTIME_POLICY,
                                  "mutations": outcomes}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(outcomes, indent=2))
    return int(any(item["status"] != "killed" for item in outcomes))


if __name__ == "__main__":
    raise SystemExit(main())
