"""Catalog values against their printed sources, wall-conduction identities and extension contracts."""

import re
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.thermodynamics import (
    CitedProperty,
    PropertyCitation,
    calculate_regenerative_cooling,
    calculate_wall_conduction,
    coolant_catalog,
    get_coolant,
    get_wall_material,
    wall_material_catalog,
)
from rocket_propulsion.thermodynamics.materials import ACCESSED, MATERIAL_SOURCES

DOCS = Path(__file__).resolve().parents[1] / "docs"

# Printed source digits (value, unit factor) transcribed in this review; see docs/REFERENCES.md.
NIST_ROWS = {
    "liquid-hydrogen": ("72.292", "9.1868", "1.4779e-05", "0.10243", "32.938", "1.2858", 20.0),
    "liquid-methane": ("424.16", "3.4660", "0.00011994", "0.18563", "190.56", "4.5992", 111.0),
    "water": ("996.96", "4.1781", "0.00085366", "0.61000", "647.10", "22.064", 300.0),
}


def exact(text, factor="1", offset="0", divisor="1"):
    return float((Fraction(text) + Fraction(offset)) * Fraction(factor) / Fraction(divisor))


def channel(**changes):
    return dict({"absorbed_heat_w": 2.0e5, "coolant_mass_flow_kg_s": 2.0,
                 "coolant_inlet_temperature_k": 111.0, "channel_count": 40,
                 "channel_flow_area_m2": 4e-6, "hydraulic_diameter_m": 2e-3,
                 "channel_length_m": 0.3, "roughness_m": 1e-6,
                 "maximum_coolant_temperature_k": 300.0}, **changes)


class CatalogValues(unittest.TestCase):
    def test_nist_rows_match_printed_digits_and_reference_condition(self):
        for key, (rho, cp, mu, k, tc, pc, temperature) in NIST_ROWS.items():
            with self.subTest(coolant=key):
                record = get_coolant(key)
                self.assertEqual(record.density.value, exact(rho))
                self.assertEqual(record.specific_heat.value, exact(cp, "1000"))
                self.assertEqual(record.dynamic_viscosity.value, exact(mu))
                self.assertEqual(record.thermal_conductivity.value, exact(k))
                self.assertEqual(record.critical_temperature.value, exact(tc))
                self.assertEqual(record.critical_pressure.value, exact(pc, "1000000"))
                for prop in (record.density, record.specific_heat, record.dynamic_viscosity,
                             record.thermal_conductivity):
                    self.assertEqual((prop.temperature_k, prop.pressure_pa), (temperature, 1.0e6))
                    self.assertEqual(prop.citation.source_id, "NIST-WEBBOOK-FLUID")

    def test_rp1_values_and_unverified_fields_are_absent(self):
        record = get_coolant("rp-1")
        self.assertEqual(record.density.value, 799.01)
        self.assertEqual(record.density.temperature_k, exact("23.283", offset="273.15"))
        self.assertEqual(record.density.pressure_pa, 83720.0)
        self.assertEqual(record.kinematic_viscosity.value, exact("1.867", "1e-6"))
        self.assertEqual(record.thermal_conductivity.value, 0.11122)
        self.assertEqual(record.thermal_conductivity.temperature_k, 302.028)
        for name in ("specific_heat", "dynamic_viscosity", "critical_temperature"):
            self.assertIsNone(getattr(record, name), name)

    def test_wall_values_and_unit_conversions_follow_cited_factors(self):
        inconel = get_wall_material("inconel-718")
        self.assertEqual(inconel.density.value, exact("0.297", "27679.90"))
        self.assertEqual(inconel.thermal_conductivity.value, exact("79", "0.1442279"))
        self.assertEqual(inconel.thermal_conductivity.temperature_k,
                         exact("70", offset="459.67", divisor="1.8"))
        self.assertEqual(inconel.maximum_service_temperature.value,
                         exact("1300", offset="459.67", divisor="1.8"))
        self.assertEqual((inconel.solidus_temperature.value, inconel.liquidus_temperature.value),
                         (1533.15, 1609.15))
        copper = get_wall_material("ofhc-copper")
        self.assertEqual((copper.thermal_conductivity.value, copper.density.value,
                          copper.solidus_temperature.value), (391.0, 8900.0, 1356.15))
        grcop = get_wall_material("grcop-84")
        self.assertEqual((grcop.thermal_conductivity.value, grcop.thermal_conductivity.upper_value),
                         (305.0, 320.0))
        self.assertIsNone(grcop.thermal_conductivity.temperature_k)
        self.assertEqual(grcop.maximum_service_temperature.value, 973.0)
        self.assertIsNone(grcop.density)
        steel = get_wall_material("stainless-304")
        self.assertEqual((steel.thermal_conductivity.value, steel.density.value,
                          steel.maximum_service_temperature.value), (16.2, 8030.0, 1172.15))

    def test_every_citation_is_registered_dated_and_indexed(self):
        index = (DOCS / "REFERENCES.md").read_text(encoding="utf-8")
        heading = "## Thermal hardware and turbopump sources (verified 2026-10-25)"
        self.assertIn(heading, index)
        for key, (_, url) in MATERIAL_SOURCES.items():
            self.assertIn(f"**{key}**", index)
            self.assertIn(url, index)
        for record in coolant_catalog() + wall_material_catalog():
            self.assertTrue(record.limitations)
            for prop in record.properties():
                with self.subTest(record=record.key, field=prop.name):
                    self.assertIn(f"accessed {ACCESSED}", prop.citation.locator)
                    self.assertIn(prop.citation.source_id, MATERIAL_SOURCES)
                    self.assertTrue(re.search(r"\d", prop.source_value))
                    if prop.conversion is not None:
                        self.assertTrue(prop.conversion.source_id.startswith("NIST-SP811"))

    def test_temperature_limit_prefers_sourced_service_temperature(self):
        self.assertEqual(get_wall_material("ofhc-copper").temperature_limit()[1], "solidus_temperature")
        for key in ("grcop-84", "inconel-718", "stainless-304"):
            self.assertEqual(get_wall_material(key).temperature_limit()[1],
                             "maximum_service_temperature")


class CatalogErrors(unittest.TestCase):
    def test_lookup_normalizes_and_rejects_unknown_or_non_text_keys(self):
        self.assertIs(get_coolant(" Liquid_Methane "), get_coolant("liquid-methane"))
        self.assertIs(get_wall_material("GRCOP_84"), get_wall_material("grcop-84"))
        for lookup, code in ((get_coolant, "unknown_coolant"),
                             (get_wall_material, "unknown_wall_material")):
            with self.assertRaises(InputError) as caught:
                lookup("unobtainium")
            self.assertEqual(caught.exception.code, code)
            with self.assertRaises(InputError):
                lookup(42)

    def test_property_records_reject_unregistered_sources_and_bad_values(self):
        good = PropertyCitation("NIST-WEBBOOK-FLUID", "row")
        with self.assertRaises(InputError):
            PropertyCitation("MADE-UP-SOURCE", "row")
        with self.assertRaises(DomainError):
            PropertyCitation("NIST-WEBBOOK-FLUID", "")
        for changes in ({"value": -1.0}, {"value": 0.0}, {"value": float("nan")},
                        {"temperature_k": 0.0}, {"upper_value": 0.5}, {"value": True}):
            fields = dict(name="density", value=1.0, unit="kg/m^3", citation=good,
                          source_value="1 kg/m3", condition="test")
            fields.update(changes)
            with self.subTest(changes=changes), self.assertRaises(DomainError):
                CitedProperty(**fields)


class WallConduction(unittest.TestCase):
    def test_hot_wall_is_coolant_wall_plus_flux_thickness_over_conductivity(self):
        r = calculate_wall_conduction(wall_material="stainless-304", wall_thickness_m=1e-3,
                                      heat_flux_w_m2=2e6, coolant_side_wall_temperature_k=400.0)
        drop = Fraction(2e6) * Fraction(1e-3) / Fraction(16.2)
        self.assertEqual(r.wall_temperature_drop_k, float(drop))
        self.assertEqual(r.hot_wall_temperature_k, float(Fraction(400) + drop))
        self.assertEqual(r.wall_limit_temperature_k, 1172.15)
        self.assertEqual(r.wall_limit_margin_k, float(Fraction(1172.15) - Fraction(400) - drop))
        self.assertFalse(r.wall_limit_exceeded)
        self.assertFalse(r.physical_model_qualified)
        self.assertEqual(r.conductivity_source, "catalog")

    def test_film_coefficient_mode_adds_coolant_film_drop(self):
        r = calculate_wall_conduction(wall_material="ofhc-copper", wall_thickness_m=8e-4,
                                      heat_flux_w_m2=5e7,
                                      coolant_side_heat_transfer_coefficient_w_m2_k=1e5,
                                      coolant_bulk_temperature_k=150.0)
        self.assertEqual(r.coolant_side_wall_temperature_k, 650.0)
        self.assertEqual(r.hot_wall_temperature_k,
                         float(Fraction(650) + Fraction(5e7) * Fraction(8e-4) / Fraction(391)))
        self.assertTrue(any(w.code == "wall_limit_is_melting_point" for w in r.warnings))

    def test_limit_exceedance_is_flagged_critical_and_exact_equality_is_not(self):
        hot = calculate_wall_conduction(wall_material="inconel-718", wall_thickness_m=1e-3,
                                        heat_flux_w_m2=1e7, coolant_side_wall_temperature_k=300.0)
        self.assertTrue(hot.wall_limit_exceeded)
        self.assertLess(hot.wall_limit_margin_k, 0)
        self.assertTrue(any(w.code == "wall_temperature_limit_exceeded" and w.severity == "critical"
                            for w in hot.warnings))
        edge = calculate_wall_conduction(wall_material="grcop-84", wall_thickness_m=1.0,
                                         heat_flux_w_m2=305.0, coolant_side_wall_temperature_k=972.0,
                                         wall_thermal_conductivity_w_m_k=305.0)
        self.assertEqual(edge.hot_wall_temperature_k, 973.0)
        self.assertFalse(edge.wall_limit_exceeded)
        self.assertEqual(edge.wall_limit_margin_k, 0.0)

    def test_range_lower_bound_and_user_override_are_reported(self):
        catalog = calculate_wall_conduction(wall_material="grcop-84", wall_thickness_m=1e-3,
                                            heat_flux_w_m2=1e7, coolant_side_wall_temperature_k=300.0)
        self.assertEqual(catalog.wall_thermal_conductivity_w_m_k, 305.0)
        self.assertIn("wall_conductivity_range_lower_bound", [w.code for w in catalog.warnings])
        user = calculate_wall_conduction(wall_material="grcop-84", wall_thickness_m=1e-3,
                                         heat_flux_w_m2=1e7, coolant_side_wall_temperature_k=300.0,
                                         wall_thermal_conductivity_w_m_k=320.0)
        self.assertEqual(user.conductivity_source, "user")
        self.assertLess(user.hot_wall_temperature_k, catalog.hot_wall_temperature_k)

    def test_ambiguous_or_incomplete_coolant_side_is_rejected(self):
        base = dict(wall_material="grcop-84", wall_thickness_m=1e-3, heat_flux_w_m2=1e7)
        for extra in ({}, {"coolant_side_wall_temperature_k": 300.0,
                           "coolant_side_heat_transfer_coefficient_w_m2_k": 1e4},
                      {"coolant_side_heat_transfer_coefficient_w_m2_k": 1e4},
                      {"coolant_side_wall_temperature_k": 300.0, "coolant_bulk_temperature_k": 100.0}):
            with self.subTest(extra=extra), self.assertRaises(InputError):
                calculate_wall_conduction(**base, **extra)
        for name in ("wall_thickness_m", "heat_flux_w_m2"):
            with self.subTest(field=name), self.assertRaises(DomainError):
                calculate_wall_conduction(**dict(base, **{name: 0.0}),
                                          coolant_side_wall_temperature_k=300.0)

    def test_record_replays_and_rejects_forged_flags(self):
        r = calculate_wall_conduction(wall_material="inconel-718", wall_thickness_m=1e-3,
                                      heat_flux_w_m2=1e7, coolant_side_wall_temperature_k=300.0)
        self.assertEqual(type(r)(**asdict(r)), r)
        for change in ({"wall_limit_exceeded": False}, {"hot_wall_temperature_k": 1.0},
                       {"physical_model_qualified": True}, {"warnings": ()}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                replace(r, **change)


class RegenerativeExtension(unittest.TestCase):
    def test_legacy_call_is_unchanged_and_has_no_extension(self):
        r = calculate_regenerative_cooling(**channel(), coolant_specific_heat_j_kg_k=3466.0,
                                           coolant_density_kg_m3=424.16,
                                           coolant_dynamic_viscosity_pa_s=0.00011994)
        self.assertIsNone(r.coolant_selection)
        self.assertIsNone(r.wall_conduction)
        self.assertEqual(r.all_warnings, r.warnings)

    def test_catalog_fill_matches_explicit_inputs_and_warns_about_reference_condition(self):
        explicit = calculate_regenerative_cooling(**channel(), coolant_specific_heat_j_kg_k=3466.0,
                                                  coolant_density_kg_m3=424.16,
                                                  coolant_dynamic_viscosity_pa_s=0.00011994)
        filled = calculate_regenerative_cooling(**channel(), coolant="liquid-methane")
        self.assertEqual(filled.inputs, explicit.inputs)
        self.assertEqual(filled.outlet_temperature_k, explicit.outlet_temperature_k)
        self.assertEqual(filled.coolant_selection.filled_fields,
                         ("coolant_specific_heat_j_kg_k", "coolant_density_kg_m3",
                          "coolant_dynamic_viscosity_pa_s"))
        codes = [w.code for w in filled.all_warnings]
        self.assertIn("coolant_catalog_reference_condition", codes)

    def test_explicit_values_override_and_missing_catalog_property_is_an_error(self):
        r = calculate_regenerative_cooling(**channel(), coolant="liquid-methane",
                                           coolant_density_kg_m3=430.0)
        self.assertEqual(r.inputs.coolant_density_kg_m3, 430.0)
        self.assertNotIn("coolant_density_kg_m3", r.coolant_selection.filled_fields)
        with self.assertRaises(InputError) as caught:
            calculate_regenerative_cooling(**channel(), coolant="rp-1")
        self.assertEqual(caught.exception.code, "coolant_property_unavailable")
        rp1 = calculate_regenerative_cooling(**channel(coolant_inlet_temperature_k=300.0,
                                                       maximum_coolant_temperature_k=600.0),
                                             coolant="rp-1", coolant_specific_heat_j_kg_k=2000.0,
                                             coolant_dynamic_viscosity_pa_s=1.5e-3)
        self.assertEqual(rp1.coolant_selection.filled_fields, ("coolant_density_kg_m3",))
        self.assertEqual(rp1.inputs.coolant_density_kg_m3, 799.01)

    def test_supercritical_outlet_is_flagged_critical(self):
        r = calculate_regenerative_cooling(**channel(absorbed_heat_w=1.0e6,
                                                     coolant_inlet_temperature_k=20.0,
                                                     maximum_coolant_temperature_k=400.0),
                                           coolant="liquid-hydrogen")
        self.assertGreater(r.outlet_temperature_k, 32.938)
        self.assertTrue(any(w.code == "coolant_outlet_above_critical_temperature"
                            and w.severity == "critical" for w in r.all_warnings))

    def test_wall_step_uses_bulk_outlet_with_film_coefficient(self):
        r = calculate_regenerative_cooling(**channel(), coolant="liquid-methane",
                                           wall_material="grcop-84", wall_thickness_m=1e-3,
                                           wall_heat_flux_w_m2=3e7,
                                           coolant_side_heat_transfer_coefficient_w_m2_k=5e4)
        wall = r.wall_conduction
        self.assertEqual(wall.inputs.coolant_bulk_temperature_k, r.outlet_temperature_k)
        expected = (Fraction(r.outlet_temperature_k) + Fraction(3e7) / Fraction(5e4)
                    + Fraction(3e7) * Fraction(1e-3) / Fraction(305))
        self.assertEqual(wall.hot_wall_temperature_k, float(expected))
        self.assertEqual(type(r)(**asdict(r)), r)

    def test_incomplete_wall_request_and_forged_extensions_are_rejected(self):
        with self.assertRaises(InputError) as caught:
            calculate_regenerative_cooling(**channel(), coolant="water", wall_material="grcop-84")
        self.assertEqual(caught.exception.code, "incomplete_wall_conduction")
        r = calculate_regenerative_cooling(**channel(), coolant="liquid-methane",
                                           wall_material="inconel-718", wall_thickness_m=1e-3,
                                           wall_heat_flux_w_m2=1e6,
                                           coolant_side_heat_transfer_coefficient_w_m2_k=5e4)
        forged_wall = calculate_wall_conduction(
            wall_material="inconel-718", wall_thickness_m=1e-3, heat_flux_w_m2=1e6,
            coolant_side_heat_transfer_coefficient_w_m2_k=5e4, coolant_bulk_temperature_k=111.0)
        with self.assertRaises(DomainError):
            replace(r, wall_conduction=forged_wall)
        with self.assertRaises(DomainError):
            replace(r, coolant_selection=replace(r.coolant_selection, warnings=()))
        other = calculate_regenerative_cooling(**channel(), coolant_specific_heat_j_kg_k=3466.0,
                                               coolant_density_kg_m3=430.0,
                                               coolant_dynamic_viscosity_pa_s=0.00011994)
        with self.assertRaises(DomainError):
            replace(other, coolant_selection=r.coolant_selection)


if __name__ == "__main__":
    unittest.main()
