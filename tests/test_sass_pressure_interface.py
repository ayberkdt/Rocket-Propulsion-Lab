"""Pressure reference and port-substitution attacks through SASS sizing."""

from dataclasses import asdict, replace
from fractions import Fraction
from json import dumps
from math import nextafter
from unittest.mock import patch

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies import sass
from rocket_propulsion.studies.sass import catalog
from tests.test_sass_studies import q


def window(low, high=None, *, reference="absolute", port="inlet", condition="synthetic flow window",
           datum=None):
    return sass.InletPressureWindow(
        port, reference, q(low, "Pa"), q(low if high is None else high, "Pa"), condition,
        None if datum is None else q(datum[0], "Pa"),
        None if datum is None else q(datum[1], "Pa"),
        None if datum is None else "synthetic reference calibration",
    )


def study(*windows, key="moog-058e151"):
    item = sass.Candidate("screen", key, q(3, "kg"), 1, 1, q(1), q(1),
                          inlet_pressures=windows)
    return sass.analyze_project(
        "Synthetic pressure comparison", q(500, "kg"),
        (sass.DeltaVItem("short burn", q(5, "m/s"), q(0), 1),), (item,),
        residual_fraction=q(0), fixed_reserve=q(0, "kg"), requirements=(),
        criteria=(sass.Criterion("wet_mass", q(500, "kg"), q(600, "kg"), q(1)),),
    )


def test_unreferenced_gauge_and_differential_do_not_match_absolute_nominal():
    for reference in ("gauge", "differential"):
        result = study(window(690000, reference=reference))
        check = result.candidates[0].catalogue_pressure_checks[0]
        assert check.status == "unassessed"
        assert check.absolute_window_pa is None
        assert result.scores[0].score is None


def test_actual_gauge_reference_is_added_at_both_endpoints():
    same = study(window(580000, reference="gauge", datum=(110000, 110000)))
    check = same.candidates[0].catalogue_pressure_checks[0]
    assert check.absolute_window_pa == (690000, 690000)
    assert check.status == "matches_nominal"
    assert same.scores[0].eligibility == "eligible"

    uncertain = study(window(580000, reference="gauge", datum=(110000, 111000)))
    check = uncertain.candidates[0].catalogue_pressure_checks[0]
    assert check.absolute_window_pa == (690000, 691000)
    assert check.status == "overlaps_reference"
    assert uncertain.scores[0].score is None


def test_nearby_window_is_not_rounded_onto_a_nominal_point():
    # A fraction of one binary64 step disappears in nearest-rounded addition.
    tiny = (nextafter(690000, float("inf")) - 690000) / 4
    assert 690000 + tiny == 690000
    result = study(window(690000, reference="gauge", datum=(tiny, tiny)))
    check = result.candidates[0].catalogue_pressure_checks[0]
    assert check.status == "outside_reference"
    lo, hi = map(Fraction, check.absolute_window_pa)
    assert lo <= Fraction(690000) + Fraction(tiny) <= hi
    assert result.candidates[0].sizing is not None
    assert result.scores[0].score is None


def test_unstated_manufacturer_datum_and_bipropellant_ports_stay_unassessed():
    mono = study(window(700000), key="ariane-1n")
    row = mono.candidates[0].catalogue_pressure_checks[0]
    assert row.status == "unassessed" and row.published.pressure_reference == "unspecified"
    bi = study(window(1400000, port="fuel_inlet"), key="ariane-s10-13")
    checks = bi.candidates[0].catalogue_pressure_checks
    assert [(item.port, item.status) for item in checks] == [
        ("fuel_inlet", "unassessed"), ("oxidizer_inlet", "unassessed"),
    ]
    assert checks[1].reason == "Inlet window not supplied."
    assert bi.scores[0].score is None


def test_declared_port_limits_cannot_substitute_for_companion_port():
    record = sass.thruster("ariane-s10-13")
    reference = record.field("inlet_pressure_range")
    # The real source does not establish these port-specific bounds.
    fuel = replace(reference, name="synthetic_fuel_pressure", pressure_port="fuel_inlet",
                   pressure_reference="absolute")
    oxidizer = replace(reference, name="synthetic_oxidizer_pressure", pressure_port="oxidizer_inlet",
                       pressure_reference="absolute", lower=replace(reference.lower, value=1_500_000))
    edited = replace(record, fields=(*record.fields, fuel, oxidizer))
    with patch.object(catalog, "THRUSTERS", tuple(edited if r is record else r for r in catalog.THRUSTERS)):
        missing = study(window(1400000, port="fuel_inlet"), key=record.key)
        off = study(window(1400000, port="fuel_inlet"), window(1400000, port="oxidizer_inlet"), key=record.key)
        good = study(window(1400000, port="fuel_inlet"), window(1500000, port="oxidizer_inlet"), key=record.key)
    assert [row.status for row in missing.candidates[0].catalogue_pressure_checks] == [
        "within_published_range", "unassessed",
    ]
    assert [row.status for row in off.candidates[0].catalogue_pressure_checks] == [
        "within_published_range", "outside_reference",
    ]
    assert missing.scores[0].score is None and off.scores[0].score is None
    assert good.scores[0].eligibility == "eligible"


def test_nominal_pressure_at_unpaired_point_is_not_reused():
    record = sass.thruster("moog-058e151")
    source = record.field("inlet_pressure")
    for changed in (
        replace(source, nominal_point_id="synthetic different point"),
        replace(source, lower=replace(source.lower, citations=(sass.Citation("ARIANE-1N", "synthetic"),))),
    ):
        edited = replace(record, fields=tuple(changed if item is source else item for item in record.fields))
        with patch.object(catalog, "THRUSTERS", tuple(edited if r is record else r for r in catalog.THRUSTERS)):
            result = study(window(690000))
        assert result.candidates[0].catalogue_pressure_checks[0].status == "unassessed"
        assert result.scores[0].score is None


def test_json_roundtrip_preserves_pressure_provenance_and_rejects_missing_datum():
    result = study(window(580000, reference="gauge", datum=(110000, 111000)))
    saved = result.to_dict()
    assert sass.project_from_json(dumps(saved["inputs"])).to_dict() == saved
    pressure = saved["candidates"][0]["catalogue_pressure_checks"][0]
    assert pressure["published"]["pressure_reference"] == "absolute"
    assert pressure["absolute_window_pa"] == [690000, 691000]
    assert "Thruster inlet pressure comparison" in sass.build_project_html(result)
    request = result.input_data
    del request["candidates"][0]["inlet_pressures"][0]["reference"]
    with pytest.raises(DomainError, match="missing required"):
        sass.project_from_json(dumps(request))
    request["candidates"][0]["inlet_pressures"] = None
    with pytest.raises(DomainError, match="JSON array"):
        sass.project_from_json(dumps(request))


def test_malformed_units_datums_and_wrong_ports_fail_closed():
    good = window(690000)
    for changes in (
        {"port": "tank"}, {"reference": ""}, {"lower": q(-1, "Pa")},
        {"upper": q(690000, "bar")}, {"upper": q(0, "Pa")},
        {"gauge_reference_record": "invented"},
    ):
        with pytest.raises(DomainError):
            replace(good, **changes)
    for args in ((-150000, -140000), (0, 20)):
        ref = (100000, 110000) if args[0] < 0 else (110000, 100000)
        with pytest.raises(DomainError):
            window(*args, reference="gauge", datum=ref)
    with pytest.raises(DomainError, match="unique ports"):
        study(good, good)
    with pytest.raises(DomainError, match="unique ports"):
        study(good, key="ariane-s10-13")
    with pytest.raises(DomainError, match="same operating condition"):
        study(window(1000000, port="fuel_inlet"), window(1100000, port="oxidizer_inlet", condition="another"),
              key="ariane-s10-13")
    assert asdict(good)["reference"] == "absolute"

