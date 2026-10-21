"""Source-bound native impulse faults; only assertion failures count as kills."""

import hashlib
import json
from pathlib import Path
import sys

from lunaris.core.propagation import plans

from scripts.run_bounded_regression import ROOT, dependency_inputs, verification_inputs
from tests.test_native_lunaris_impulse_numerics import NativeLunarisImpulseAttacks


MUTATIONS = (
    (
        "native_impulse_norm_gain_error",
        "_impulse_norm_and_mass",
        "norm_decimal = sum((value * value for value in coordinates), Decimal(0)).sqrt()",
        'norm_decimal = Decimal("1.01") * sum((value * value for value in coordinates), Decimal(0)).sqrt()',
        "test_ordinary_ric_rotation_and_mass_law_are_preserved",
    ),
    (
        "native_impulse_early_g0_isp_product",
        "_impulse_norm_and_mass",
        "Decimal.from_float(STANDARD_GRAVITY) * Decimal.from_float(isp_s)",
        "Decimal.from_float(STANDARD_GRAVITY * isp_s)",
        "test_large_isp_product_is_not_rounded_to_infinity",
    ),
    (
        "native_impulse_accept_unrepresentable_norm",
        "_impulse_norm_and_mass",
        "if not math.isfinite(dv_norm):",
        "if False and not math.isfinite(dv_norm):",
        "test_genuine_scalar_norm_overflow_has_an_explicit_error",
    ),
    (
        "native_impulse_accept_erased_final_mass",
        "_impulse_norm_and_mass",
        "if not math.isfinite(mass_after) or mass_after <= 0.0:",
        "if False and (not math.isfinite(mass_after) or mass_after <= 0.0):",
        "test_genuine_final_mass_underflow_has_an_explicit_error",
    ),
    (
        "native_impulse_accept_velocity_overflow",
        "apply_impulsive_maneuver",
        "if not np.all(np.isfinite(velocity_after)):",
        "if False and not np.all(np.isfinite(velocity_after)):",
        "test_genuine_velocity_update_overflow_has_an_explicit_error",
    ),
    (
        "native_impulse_hide_mass_change_in_record",
        "apply_impulsive_maneuver",
        '"mass_after_kg": mass_after,',
        '"mass_after_kg": mass_before,',
        "test_records_match_state_and_remain_strict_json",
    ),
)


def main():
    inputs = verification_inputs(ROOT)
    dependency_roots = [Path(plans.__file__).resolve().parents[2]]
    external = dependency_inputs(dependency_roots)
    outcomes = []
    for name, symbol, before, after, test in MUTATIONS:
        function = getattr(plans, symbol)
        from scripts.mutation_runtime import compile_fault, evaluate_fault, MUTATION_RUNTIME_POLICY
        source,replacement = compile_fault(function,before,after,f"<native-mutation:{name}>")
        def bindings():return [
            (module, key)
            for module_name, module in list(sys.modules.items())
            if module_name.startswith(("lunaris.", "rocket_propulsion.", "tests."))
            and module is not None
            for key, value in list(vars(module).items())
            if value is function
        ]
        trial=evaluate_fault(bindings,function,replacement,NativeLunarisImpulseAttacks,test)
        outcomes.append(
            {
                "mutation": name,
                **trial,
                "symbol": f"{plans.__name__}.{symbol}",
                "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                "test": f"NativeLunarisImpulseAttacks.{test}",
            }
        )
    if inputs != verification_inputs(ROOT) or external != dependency_inputs(dependency_roots):
        raise ValueError("Native mutation inputs changed during evaluation.")
    output = ROOT / "output/assurance/native_lunaris_mutations.json"
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(
        json.dumps(
            {
                "scope": "six curated native impulse numerical/trace faults; not exhaustive",
                "source_sha256": inputs,
                "runtime_policy": MUTATION_RUNTIME_POLICY,
                "external_dependency_roots": [str(path) for path in dependency_roots],
                "external_dependency_sha256": external,
                "mutations": outcomes,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(outcomes, indent=2))
    return int(any(row["status"] != "killed" for row in outcomes))


if __name__ == "__main__":
    raise SystemExit(main())
