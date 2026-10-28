"""Independent high-precision element-potential finite-difference Jacobians."""

import math
import unittest
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.thermochemistry import equilibrium as model
from rocket_propulsion.thermochemistry.species import get_species


class EquilibriumJacobianAttacks(unittest.TestCase):
    def test_initial_jacobian_matches_independent_high_precision_residuals(self):
        candidates=tuple(get_species(key) for key in ('H2','O2','H2O','OH','H','O','CO','CO2','CH4'))
        temperature=3200.0;pressure=7e6
        gibbs=[]
        for species in candidates:
            state=species.properties(temperature)
            gibbs.append((state.enthalpy_molar_j_mol-temperature*state.entropy_molar_j_mol_k)
                         /(model.UNIVERSAL_GAS_CONSTANT_J_MOL_K*temperature))
        atoms=[dict(species.formula) for species in candidates]
        pressure_term=math.log(model.REFERENCE_PRESSURE_PA/pressure)
        inventory={'C':1.0,'H':4.0,'O':4.0}
        linear=model._solve_linear
        for elements in (('C','H','O'),('H','O','C'),('O','C','H')):
            with self.subTest(elements=elements),mp.workdps(100):
                def residual(potentials,*,elements=elements):
                    weights=[mp.exp(mp.mpf(pressure_term)-g+mp.fsum(
                        mp.mpf(record.get(e,0))*p for e,p in zip(elements,potentials,strict=True)))
                        for record,g in zip(atoms,gibbs,strict=True)]
                    total=mp.fsum(weights)
                    amounts={e:mp.fsum(w*record.get(e,0) for w,record in zip(weights,atoms,strict=True))
                             for e in elements}
                    return [mp.log(total)]+[mp.log(amounts[e]/amounts[elements[0]])
                        -mp.log(mp.mpf(inventory[e])/inventory[elements[0]]) for e in elements[1:]]
                zero=[mp.mpf(0)]*3;baseline=residual(zero);step=mp.mpf(1e-5)
                expected=[]
                for column in range(3):
                    shifted=zero.copy();shifted[column]+=step
                    expected.append([(value-base)/step for value,base in zip(residual(shifted),baseline,strict=True)])
                calls=[]

                def capture(matrix,rhs,*,expected=expected,calls=calls,baseline=baseline):
                    if not calls:
                        for row in range(3):
                            self.assertAlmostEqual(rhs[row],-float(baseline[row]),delta=1e-10)
                            for column in range(3):
                                self.assertAlmostEqual(matrix[row][column],float(expected[column][row]),delta=2e-7)
                    calls.append(1)
                    return linear(matrix,rhs)

                provider=model.LocalEquilibriumProvider()
                with patch.object(model,'_solve_linear',side_effect=capture):
                    provider._composition_at_tp(inventory={e:inventory[e] for e in elements},
                        candidates=candidates,temperature_k=temperature,pressure_pa=pressure)
                self.assertTrue(calls)
