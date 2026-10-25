"""Curated shared-root faults with assertion kills and passing before/after controls.

Uses the existing live-global mutation runtime. This campaign neither replaces
nor renews the historical feed/SASS/native campaigns on a different source tree.
"""

import hashlib
import json
import sys
from pathlib import Path

from rocket_propulsion.core import numerics as model
from scripts.mutation_runtime import MUTATION_RUNTIME_POLICY, compile_fault, evaluate_fault
from scripts.run_bounded_regression import verification_inputs
from tests.test_root_contract import RootAttacks

ROOT = Path(__file__).resolve().parents[1]
MUTATIONS = (
    ('boolean_target', model, '_controls', "real(target,'target')",
     "real(int(target) if isinstance(target,bool) else target,'target')",
     'test_invalid_controls_fail_before_any_callback'),
    ('boolean_tolerance', model, '_controls', "real(tolerance,'tolerance')",
     "real(int(tolerance) if isinstance(tolerance,bool) else tolerance,'tolerance')",
     'test_invalid_controls_fail_before_any_callback'),
    ('truncate_iteration_budget', model, '_controls', "integer(max_iterations,'max_iterations')",
     "integer(int(max_iterations) if isinstance(max_iterations,float) else max_iterations,'max_iterations')",
     'test_invalid_controls_fail_before_any_callback'),
    ('overflow_midpoint', model, '_midpoint', 'float((F(lower)+F(upper))/2)',
     '(lower+upper)/2', 'test_huge_same_sign_interval_does_not_overflow_its_midpoint'),
    ('rounded_residual_decision', model, '_accepted',
     '_residual(value,target) <= F(tolerance)*max(F(1),abs(F(target)))',
     'abs(value-target) <= tolerance*max(1,abs(target))',
     'test_exact_residual_decision_is_not_made_from_rounded_subtraction'),
    ('double_residual_allowance', model, '_accepted', 'F(tolerance)*max(F(1),abs(F(target)))',
     '2*F(tolerance)*max(F(1),abs(F(target)))',
     'test_exact_residual_decision_is_not_made_from_rounded_subtraction'),
    ('tolerance_expands_bracket', model, '_run',
     'min(low_value,high_value) <= target <= max(low_value,high_value)',
     'min(low_value,high_value)-tolerance <= target <= max(low_value,high_value)+tolerance',
     'test_target_outside_actual_bracket_is_not_admitted_by_tolerance'),
    ('always_increasing', model, '_run', 'increasing=high_value>=low_value', 'increasing=True',
     'test_increasing_and_decreasing_exact_affine_oracles'),
    ('ignore_observed_monotonicity', model, '_run',
     'if not min(low_value,high_value) <= value <= max(low_value,high_value):',
     'if False and not min(low_value,high_value) <= value <= max(low_value,high_value):',
     'test_monotonicity_violating_observation_is_refused'),
    ('halve_reported_residual', model, '_pack', 'residual=_residual(observed,target)',
     'residual=_residual(observed,target)/2',
     'test_increasing_and_decreasing_exact_affine_oracles'),
    ('overwrite_upper_observation', model, '_pack', "'upper_function_value':high_value",
     "'upper_function_value':low_value",
     'test_adjacent_float_resolution_failure_preserves_signs_without_relaxing_tolerance'),
    ('hide_float_resolution', model, '_run', "_failure('root_float_resolution',iteration-1",
     "_failure('non_convergence',iteration-1",
     'test_adjacent_float_resolution_failure_preserves_signs_without_relaxing_tolerance'),
    ('invent_failure_iteration', model, '_failure', '_pack(value,observed,iterations,',
     '_pack(value,observed,iterations+1,',
     'test_budget_is_actual_integer_and_failure_counts_real_observations'),
    ('alias_trace_records', model, '_owned_trace',
     'return tuple(RootEvaluation(row.x,row.value) for row in rows)', 'return rows',
     'test_result_replay_does_not_call_the_external_function_again'),
    ('unvalidated_evaluation_record', model.RootEvaluation, '__post_init__',
     "coerce_fields(self,'Root evaluation')", 'return None',
     'test_public_evaluation_records_refuse_nonfinite_and_boolean_fields'),
    ('ignore_derived_field_replay', model.RootResult, '__post_init__',
     "any(getattr(self,name)!=value for name,value in expected.items() if name!='residual_ratio')",
     "any(False for name,value in expected.items() if name!='residual_ratio')",
     'test_derived_fields_and_trace_coordinates_cannot_be_forged'),
    ('allow_qualification_flag', model.RootResult, '__post_init__',
     'if type(getattr(self,name)) is not bool or getattr(self,name):',
     'if False and (type(getattr(self,name)) is not bool or getattr(self,name)):',
     'test_false_evidence_flags_and_coordinate_bound_cannot_be_promoted'),
    ('invent_coordinate_bound', model.RootResult, '__post_init__',
     "self.coordinate_error_bound is not None or self.method!='bisection'",
     "self.method!='bisection'",
     'test_false_evidence_flags_and_coordinate_bound_cannot_be_promoted'),
    ('falsify_owned_observation', model, 'solve_bisect_monotonic',
     'rows.append(RootEvaluation(x,value))', 'rows.append(RootEvaluation(x,value+1))',
     'test_result_replay_does_not_call_the_external_function_again'),
    ('mislabel_failure_as_tolerance_success', model, '_failure', 'termination=code',
     "termination='residual_tolerance'",
     'test_budget_is_actual_integer_and_failure_counts_real_observations'),
    ('claim_failure_met_acceptance', model, '_failure', 'acceptance_satisfied=False',
     'acceptance_satisfied=True',
     'test_adjacent_float_resolution_failure_preserves_signs_without_relaxing_tolerance'),
)


def main():
    """Keep errors/survivors visible; never count them as assertion kills."""
    fingerprints = verification_inputs(ROOT)
    outcomes = []
    for name, owner, symbol, before, after, test_name in MUTATIONS:
        original = getattr(owner, symbol)
        source, replacement = compile_fault(original, before, after, f'<root-fault:{name}>')

        def bindings(*, owner=owner, symbol=symbol, original=original):
            targets = [(owner, symbol)]
            if owner is model:
                targets.extend((module, key)
                               for module_name, module in list(sys.modules.items())
                               if module is not None and module_name.startswith(
                                   ('rocket_propulsion.', 'tests.'))
                               for key, value in list(vars(module).items())
                               if value is original and (module is not owner or key != symbol))
            return targets

        trial = evaluate_fault(bindings, original, replacement, RootAttacks, test_name)
        outcomes.append({'mutation': name, 'symbol': owner.__name__ + '.' + symbol,
                         'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
                         'test': 'RootAttacks.' + test_name, **trial})
        print(name, trial['status'], flush=True)
    unchanged = verification_inputs(ROOT) == fingerprints
    passed = unchanged and all(row['status'] == 'killed' for row in outcomes)
    report = {'scope': f'{len(MUTATIONS)} curated shared-root faults; not exhaustive coverage',
              'runtime_policy': MUTATION_RUNTIME_POLICY, 'source_sha256': fingerprints,
              'source_unchanged_during_run': unchanged, 'campaign_passed': passed,
              'mutations': outcomes}
    output = ROOT / 'output/assurance/root_mutations.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
