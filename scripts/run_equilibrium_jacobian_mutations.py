"""Five targeted residual/Jacobian faults; no global equilibrium qualification inferred."""

import hashlib
import json
from pathlib import Path

from rocket_propulsion.thermochemistry import equilibrium as model
from scripts.mutation_runtime import MUTATION_RUNTIME_POLICY, compile_fault, evaluate_fault
from scripts.run_bounded_regression import verification_inputs
from tests.test_equilibrium_jacobian import EquilibriumJacobianAttacks

ROOT=Path(__file__).resolve().parents[1]
MUTATIONS=(
    ('repeat_first_column','shifted_columns[column][row]','shifted_columns[0][row]'),
    ('transpose_jacobian','shifted_columns[column][row]','shifted_columns[row][column]'),
    ('reverse_difference_direction','shifted[column] += step_size','shifted[column] -= step_size'),
    ('halve_jacobian','(shifted_columns[column][row] - values[row]) / step_size',
     '(shifted_columns[column][row] - values[row]) / (2*step_size)'),
    ('reverse_newton_rhs','[-value for value in values]','[value for value in values]'),
)


def main():
    fingerprints=verification_inputs(ROOT);outcomes=[]
    owner=model.LocalEquilibriumProvider;symbol='_composition_at_tp'
    original=getattr(owner,symbol)
    for name,before,after in MUTATIONS:
        source,replacement=compile_fault(original,before,after,f'<jacobian-fault:{name}>')
        trial=evaluate_fault(lambda:[(owner,symbol)],original,replacement,EquilibriumJacobianAttacks,
            'test_initial_jacobian_matches_independent_high_precision_residuals')
        outcomes.append({'mutation':name,'symbol':'LocalEquilibriumProvider._composition_at_tp',
                         'source_sha256':hashlib.sha256(source.encode()).hexdigest(),**trial})
        print(name,trial['status'],flush=True)
    unchanged=verification_inputs(ROOT)==fingerprints
    passed=unchanged and all(row['status']=='killed' for row in outcomes)
    report={'scope':'Five element-potential residual/finite-difference faults only; not exhaustive',
            'runtime_policy':MUTATION_RUNTIME_POLICY,'source_sha256':fingerprints,
            'source_unchanged_during_run':unchanged,'campaign_passed':passed,'mutations':outcomes}
    output=ROOT/'output/assurance/equilibrium_jacobian_mutations.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return int(not passed)


if __name__=='__main__':
    raise SystemExit(main())
