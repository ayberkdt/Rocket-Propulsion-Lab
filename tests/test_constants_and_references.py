"""Constants, atmosphere and reference discipline, checked like an auditor.

These tests try to catch (1) a physical constant copied as a bare literal and
drifting from its registered source, (2) a model that uses a constant other
than the one its defining document uses, and (3) documentation that states
equations without saying where they come from.

Oracles and their sources:

* SI exact constants: N_A = 6.02214076e23 /mol, k = 1.380649e-23 J/K, so
  R = N_A k; g_n = 9.80665 m/s^2; T(degC) = T(K) - 273.15.
  BIPM SI Brochure, 9th ed., 2019: https://www.bipm.org/en/publications/si-brochure
* U.S. Standard Atmosphere 1976: R* = 8.31432 J/(mol K), M0 = 0.0289644
  kg/mol, P0 = 101325 Pa, T0 = 288.15 K, layer bases H_b, T_b, L_b, and the
  hydrostatic equation dP/dH = -g0 M0 P / (R* T(H)).  Layer base pressures
  printed in the report: 22632.06, 5474.889, 868.0187, 110.9063, 66.93887,
  3.956420 Pa at 11, 20, 32, 47, 51, 71 km'.
  NOAA-S/T 76-1562 / NASA-TM-X-74335: https://ntrs.nasa.gov/citations/19770009539
* The hydrostatic oracle is integrated here with classical RK4, independent of
  the closed forms in the library (any Runge-Kutta text; e.g. NIST DLMF
  §3.7(v), https://dlmf.nist.gov/3.7#v).
"""

import ast
import inspect
import re
import unittest
from itertools import pairwise
from math import exp, isclose, log
from pathlib import Path
from typing import ClassVar

from rocket_propulsion.compressible.off_design_nozzle import standard_atmosphere_pressure_pa
from rocket_propulsion.core import constants
from rocket_propulsion.core.constants import REGISTERED_CONSTANTS
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.integrations import lunaris
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "rocket_propulsion"
DOCS = ROOT / "docs"

# USSA 1976 layers: (H_b [m'], T_b [K], L_b [K/m']) — defining values.
USSA_LAYERS = (
    (0.0, 288.15, -0.0065),
    (11_000.0, 216.65, 0.0),
    (20_000.0, 216.65, 0.0010),
    (32_000.0, 228.65, 0.0028),
    (47_000.0, 270.65, 0.0),
    (51_000.0, 270.65, -0.0028),
    (71_000.0, 214.65, -0.0020),
)
USSA_PRINTED_BASE_PRESSURES = (101_325.0, 22_632.06, 5_474.889, 868.0187, 110.9063, 66.93887,
                               3.956420)
USSA_TOP_M = 84_852.0


def ussa_temperature(height: float) -> float:
    base, temperature, lapse = USSA_LAYERS[0]
    for layer in USSA_LAYERS:
        if layer[0] <= height:
            base, temperature, lapse = layer
    return temperature + lapse * (height - base)


def hydrostatic_rk4(height: float, steps_per_km: int = 40) -> float:
    """Integrate ln P from sea level with RK4, splitting at every layer base."""

    g0 = 9.80665
    m0 = 0.0289644
    r_star = 8.31432

    def slope(h: float) -> float:
        return -g0 * m0 / (r_star * ussa_temperature(h))

    bounds = sorted({0.0, height, *[layer[0] for layer in USSA_LAYERS if layer[0] < height]})
    log_p = log(101_325.0)
    for start, end in pairwise(bounds):
        n = max(1, int((end - start) / 1000.0 * steps_per_km))
        h = (end - start) / n
        for k in range(n):
            a = start + k * h
            # d(ln P)/dh depends on h only, so RK4 reduces to Simpson's rule;
            # T(h) is continuous across layer bases.
            log_p += h / 6 * (slope(a) + 4 * slope(a + h / 2) + slope(a + h))
    return exp(log_p)


class ConstantRegistryTests(unittest.TestCase):
    def test_values_match_their_defining_sources(self) -> None:
        self.assertEqual(constants.STANDARD_GRAVITY_M_S2, 9.80665)
        self.assertEqual(STANDARD_GRAVITY_M_S2, constants.STANDARD_GRAVITY_M_S2)
        # R = N_A k exactly in the 2019 SI (to one rounding of the product).
        self.assertTrue(isclose(constants.MOLAR_GAS_CONSTANT_J_MOL_K,
                                6.02214076e23 * 1.380649e-23, rel_tol=2e-16))
        self.assertEqual(constants.USSA1976_GAS_CONSTANT_J_MOL_K, 8.31432)
        self.assertEqual(constants.USSA1976_MOLAR_MASS_KG_MOL, 0.0289644)
        # The rounded air constant must stay within its rounding of R*/M0.
        exact_air = constants.USSA1976_GAS_CONSTANT_J_MOL_K / constants.USSA1976_MOLAR_MASS_KG_MOL
        self.assertLess(abs(constants.AIR_SPECIFIC_GAS_CONSTANT_J_KG_K - exact_air), 0.005)
        # The two gas constants really differ (by 1.7e-5): mixing them is a bug.
        ratio = constants.MOLAR_GAS_CONSTANT_J_MOL_K / constants.USSA1976_GAS_CONSTANT_J_MOL_K
        self.assertAlmostEqual(ratio - 1.0, 1.7e-5, delta=0.1e-5)

    def test_every_constant_names_a_source_and_unit(self) -> None:
        names = set()
        for constant in REGISTERED_CONSTANTS:
            self.assertTrue(constant.source.strip(), constant.name)
            self.assertTrue(constant.unit.strip(), constant.name)
            self.assertTrue(constant.url.startswith("https://"), constant.name)
            self.assertNotIn(constant.name, names)
            names.add(constant.name)

    def test_no_registered_value_is_copied_as_a_bare_literal(self) -> None:
        values = {constant.value: constant.name for constant in REGISTERED_CONSTANTS}
        offenders = []
        for path in PACKAGE.rglob("*.py"):
            if path.name == "constants.py" and path.parent.name == "core":
                continue
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if (
                    isinstance(node, ast.Constant)
                    and isinstance(node.value, (int, float))
                    and not isinstance(node.value, bool)
                    and float(node.value) in values
                ):
                    offenders.append(f"{path.relative_to(ROOT)}:{node.lineno} "
                                     f"{node.value!r} ({values[float(node.value)]})")
        self.assertEqual(offenders, [], "Use rocket_propulsion.core.constants instead")


class StandardAtmosphereTests(unittest.TestCase):
    def test_printed_base_pressures_are_reproduced_from_both_sides(self) -> None:
        for (base, _, _), printed in zip(USSA_LAYERS, USSA_PRINTED_BASE_PRESSURES, strict=True):
            at_base = standard_atmosphere_pressure_pa(base)
            self.assertEqual(at_base, printed, base)
            if base > 0.0:
                below = standard_atmosphere_pressure_pa(base - 1e-7)
                # Continuity is limited only by the 7 printed digits.
                self.assertLess(abs(below - printed) / printed, 3e-7, base)

    def test_closed_forms_match_an_independent_hydrostatic_integration(self) -> None:
        for height in (0.0, 500.0, 10_999.0, 11_000.0, 15_432.1, 20_000.0, 28_000.0,
                       32_000.0, 40_001.0, 47_000.0, 49_000.0, 51_000.0, 60_000.0,
                       71_000.0, 80_000.0, USSA_TOP_M):
            oracle = hydrostatic_rk4(height)
            value = standard_atmosphere_pressure_pa(height)
            # Oracle starts from P0 and accumulates the printed-digit rounding
            # of the layer bases used by the library: 3e-7 relative.
            self.assertLess(abs(value - oracle) / oracle, 3e-7, height)

    def test_pressure_is_strictly_decreasing_and_bounded(self) -> None:
        previous = float("inf")
        for k in range(8486):
            height = min(10.0 * k, USSA_TOP_M)
            value = standard_atmosphere_pressure_pa(height)
            self.assertLess(value, previous)
            previous = value
        for bad in (-1e-9, USSA_TOP_M + 1e-6, float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(DomainError):
                standard_atmosphere_pressure_pa(bad)


class DocumentationReferenceTests(unittest.TestCase):
    # Index, plan and changelog pages carry no model equations of their own.
    EXEMPT: ClassVar[frozenset[str]] = frozenset({
        "README.md", "RELEASE_NOTES.md", "ROADMAP.md", "DEVELOPMENT_PLAN.md", "GAP_ANALYSIS.md",
        "API.md", "ARCHITECTURE.md", "USER_GUIDE.md", "REFERENCES.md",
        # A design plan; its sources are cited inline and indexed in REFERENCES.md.
        "SIDERA_BURN_INTEGRATION_PLAN.md",
    })

    def test_every_model_page_has_a_reference_section_with_sources(self) -> None:
        missing = []
        for path in sorted(DOCS.glob("*.md")):
            if path.name in self.EXEMPT:
                continue
            text = path.read_text(encoding="utf-8")
            match = re.search(r"^#{2,3} .*references?\b.*$", text,
                              flags=re.MULTILINE | re.IGNORECASE)
            if match is None:
                missing.append(f"{path.name}: no References heading")
                continue
            section = text[match.end():]
            if "https://" not in section and "(bibliographic" not in section:
                missing.append(f"{path.name}: References section cites nothing")
        self.assertEqual(missing, [])

    def test_every_reference_url_in_docs_is_listed_in_the_reference_index(self) -> None:
        index = (DOCS / "REFERENCES.md").read_text(encoding="utf-8")
        listed = set(re.findall(r"https://[^\s>)]+", index))
        unlisted = set()
        for path in DOCS.glob("*.md"):
            for url in re.findall(r"https://[^\s>)]+", path.read_text(encoding="utf-8")):
                url = url.rstrip(".,;")
                host = url.split("/")[2]
                # Code-hosting and API pages are not engineering sources.
                if host in {"github.com", "www.orekit.org", "dlmf.nist.gov"}:
                    continue
                if url not in listed:
                    unlisted.add(f"{path.name}: {url}")
        self.assertEqual(sorted(unlisted), [])

    def test_every_registered_source_is_listed_in_the_reference_index(self) -> None:
        from rocket_propulsion.propulsion.burns import BURN_REFERENCES

        index = (DOCS / "REFERENCES.md").read_text(encoding="utf-8")
        missing = [f"{item.reference_id}: {item.url}" for item in BURN_REFERENCES
                   if item.url not in index]
        missing += [f"{item.name}: {item.url}" for item in REGISTERED_CONSTANTS
                    if item.url not in index]
        self.assertEqual(missing, [])

    def test_astrodynamics_adapter_api_carries_references(self) -> None:
        missing = []
        for module in (lunaris, constants):
            for name, value in inspect.getmembers(module):
                if name.startswith("_") or getattr(value, "__module__", None) != module.__name__:
                    continue
                if (inspect.isclass(value) or inspect.isfunction(value)) and "References" not in (
                    inspect.getdoc(value) or ""
                ):
                    missing.append(f"{module.__name__}.{name}")
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
