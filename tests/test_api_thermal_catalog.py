"""Nozzle-thermal route with catalog coolants and wall conduction."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.errors import InputError

INDEX = Path(__file__).resolve().parents[1] / "src" / "rocket_propulsion" / "web" / "index.html"


def _form_defaults() -> dict[str, str]:
    """Read the browser form's default payload, as the page submits it."""
    text = INDEX.read_text(encoding="utf-8")
    start = text.index('id="nozzle-thermal-form"')
    form = text[start:text.index("</form>", start)]
    payload = {m.group(1): m.group(2) for m in re.finditer(r'<input name="([^"]+)" type="[^"]+" value="([^"]*)"', form)}
    for match in re.finditer(r'<select name="([^"]+)"[^>]*>(.*?)</select>', form, re.DOTALL):
        chosen = re.search(r'<option value="([^"]*)" selected', match.group(2))
        payload[match.group(1)] = chosen.group(1) if chosen else re.search(r'<option value="([^"]*)"', match.group(2)).group(1)
    return payload


class ThermalCatalogRouteTests(unittest.TestCase):
    def test_default_form_keeps_hand_entered_properties_and_no_wall(self) -> None:
        cooling = dispatch_calculation("/api/v1/nozzle-thermal", _form_defaults())["cooling"]
        self.assertIsNone(cooling["coolant_selection"])
        self.assertIsNone(cooling["wall_conduction"])

    def test_catalog_coolant_fills_properties_and_wall_limit_is_screened(self) -> None:
        payload = _form_defaults()
        for key in ("coolant_density_kg_m3", "coolant_dynamic_viscosity_pa_s", "coolant_specific_heat_j_kg_k"):
            payload.pop(key)
        payload.update(coolant="liquid-hydrogen", wall_material="grcop-84")
        cooling = dispatch_calculation("/api/v1/nozzle-thermal", payload)["cooling"]
        self.assertEqual(set(cooling["coolant_selection"]["filled_fields"]),
                         {"coolant_specific_heat_j_kg_k", "coolant_density_kg_m3", "coolant_dynamic_viscosity_pa_s"})
        wall = cooling["wall_conduction"]
        self.assertGreater(wall["hot_wall_temperature_k"], wall["coolant_side_wall_temperature_k"])
        self.assertEqual(wall["wall_limit_exceeded"], wall["hot_wall_temperature_k"] > wall["wall_limit_temperature_k"])

    def test_non_string_catalog_key_is_refused(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/nozzle-thermal", {**_form_defaults(), "coolant": 3})


if __name__ == "__main__":
    unittest.main()
