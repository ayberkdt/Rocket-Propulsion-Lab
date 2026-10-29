"""Optional Cantera 3.2.0 cross-implementation reference; not physical data.

Run with the isolated Cantera toolchain on PYTHONPATH. The unit suite reads the
frozen JSON and does not import Cantera. The independent source is the pinned
NASA fixed-column extraction, not the production property evaluator.
"""
import hashlib
import json
import platform
from fractions import Fraction as F
from pathlib import Path

import cantera as ct

ROOT = Path(__file__).resolve().parents[1]


def main():
    if ct.__version__ != '3.2.0':
        raise RuntimeError('Reference generation requires Cantera 3.2.0.')
    path = ROOT/'docs/assurance/nasa9_catalogue_reference.json'
    data = json.loads(path.read_text(encoding='utf8'))
    species = []
    masses = {}
    for row in data['species']:
        item = ct.Species(row['key'],row['formula'])
        coefficients = [2]
        for interval in row['intervals']:
            coefficients.extend((interval['minimum_temperature_k'],interval['maximum_temperature_k'],
                                 *map(float,interval['coefficients'])))
        item.thermo = ct.Nasa9PolyMultiTempRegion(200,6000,100_000,coefficients)
        species.append(item)
        masses[row['key']] = F(float(F(float(row['molecular_mass_g_mol']))/1000))
    gas = ct.Solution(thermo='ideal-gas',species=species)
    mixtures = [(v,b) for v in ({'N2':.79,'O2':.21},{'H2O':.7,'CO2':.3},{'H2':.2,'O2':.3,'OH':.5}) for b in ('mole','mass')]
    mixtures.extend(({k:1},'mole') for k in ('H2','O2','H2O'))
    mixtures.append(({s.name:.1 for s in species},'mole'))
    rows = []
    for composition,basis in mixtures:
        amount = {k:F(v) if basis=='mole' else F(v)/masses[k] for k,v in composition.items()}
        total = sum(amount.values())
        fractions = {k:float(v/total) for k,v in amount.items()}
        for temperature in (200,298.15,1000,3500,6000):
            for pressure in (100,100_000,2e7):
                gas.TPX = temperature,pressure,fractions
                rows.append({'input_fractions':composition,'basis':basis,'temperature_k':temperature,
                             'pressure_pa':pressure,'cantera_mole_fractions':fractions,
                             'cp_over_R':gas.cp_mole/ct.gas_constant,
                             'enthalpy_over_R_k':gas.enthalpy_mole/ct.gas_constant,
                             'entropy_over_R':gas.entropy_mole/ct.gas_constant})
    native = Path(ct._cantera.__file__)
    report = ROOT/'output/assurance/cantera_install_report.json'
    installation = json.loads(report.read_text(encoding='utf8')) if report.exists() else None
    package = (next(row for row in installation['install'] if row['metadata']['name'].lower()=='cantera')
               if installation else None)
    document = {'scope':'Independent Cantera implementation using the same pinned NASA fits. Not independent physical data, gas/hardware validation, equilibrium or certification.',
                'cantera_version':ct.__version__,'python':platform.python_version(),
                'native_binary_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),
                'wheel_sha256':package['download_info']['archive_info']['hashes']['sha256'] if package else None,
                'source_fixture_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'constant_and_mass_convention':'Dimensionless molar properties remove Cantera versus NASA defining-R differences; mass-basis inputs use the pinned historical NASA molar masses. No output fitting.',
                'tolerances':{'relative':2e-12,'cp_over_R_absolute':5e-12,'enthalpy_over_R_k_absolute':3e-8,'entropy_over_R_absolute':5e-11},
                'references':['https://cantera.org/3.2/cxx/d7/dfa/classCantera_1_1IdealGasPhase.html',
                              'https://ntrs.nasa.gov/citations/20020085330'],
                'cases':rows}
    output = ROOT/'docs/assurance/mixture_cantera_reference.json'
    output.write_text(json.dumps(document,indent=2)+'\n',encoding='utf8')
    print(f'Wrote {len(rows)} independent Cantera states.')


if __name__ == '__main__':
    main()
