"""Keep catalogue operating, proof and burst pressure roles separate."""

from dataclasses import replace
from decimal import Decimal, localcontext
from math import nextafter
from unittest.mock import patch

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies import sass
from rocket_propulsion.studies.sass import catalog
from tests.test_sass_pressure_interface import study, window


def with_moog_fields(record, fields):
    updated = replace(record, fields=fields)
    return patch.object(catalog, "THRUSTERS", tuple(updated if entry is record else entry
                                                      for entry in catalog.THRUSTERS))


def test_moog_pressure_ratings_keep_their_published_roles_and_scope():
    item = sass.thruster("moog-058e151")
    for name, kind, limit in (
        ("inlet_pressure", "nominal", 690000),
        ("maximum_expected_operating_pressure", "upper_bound", 2760000),
        ("proof_pressure", "nominal", 7000000),
        ("burst_pressure", "nominal", 11140000),
    ):
        field = item.field(name)
        assert (field.kind, field.upper.value, field.upper.unit) == (kind, limit, "Pa")
        assert field.pressure_reference == "absolute" and field.pressure_port == "inlet"
        assert {citation.source_id for citation in field.upper.citations} == {"MOOG-COLD"}
    rating = item.field("maximum_expected_operating_pressure")
    assert (rating.reported_psia.value, rating.reported_psia.unit) == (400, "psia")
    assert rating.reported_psia.citations == rating.upper.citations
    valid = study(window(690000))
    pressure = valid.candidates[0].catalogue_pressure_checks[0]
    assert pressure.status == "matches_nominal"
    assert pressure.meop.status == "within_meop_reference"
    assert valid.scores[0].eligibility == "eligible"


def test_outside_nominal_but_below_meop_has_two_distinct_results():
    result = study(window(700000, 2000000))
    pressure = result.candidates[0].catalogue_pressure_checks[0]
    assert pressure.status == "outside_reference"
    assert pressure.meop.status == "within_meop_reference"
    assert result.candidates[0].sizing is not None
    assert result.scores[0].score is None


def test_meop_classifies_exact_boundary_crossing_and_exceedance():
    above = nextafter(2760000, float("inf"))
    for low, high, expected in (
        (2758000, 2760000, "meop_unit_disagreement"),
        (2760000, 2760000, "meop_unit_disagreement"),
        (2757000, 2758000, "crosses_meop_reference"),
        (2760000, above, "crosses_meop_reference"),
        (above, above, "exceeds_meop_reference"),
        (5000000, 6000000, "exceeds_meop_reference"),
    ):
        result = study(window(low, high))
        pressure = result.candidates[0].catalogue_pressure_checks[0]
        assert pressure.meop.status == expected
        assert pressure.meop.published_limit.name == "maximum_expected_operating_pressure"
        assert result.scores[0].score is None


def test_meop_keeps_exact_conversion_between_adjacent_binary_values():
    # NIST's exact lbf divided by the exact square-inch area, independently
    # evaluated in Decimal. Gauge addition can land between two Pa floats.
    with localcontext() as context:
        context.prec = 70
        exact = Decimal(400) * Decimal("4.4482216152605") / Decimal("0.00064516")
        lower = float(exact)
        upper = nextafter(lower, float("inf"))
        assert Decimal.from_float(lower) < exact < Decimal.from_float(upper)
        close_below, close_above = 0.91726734450, 0.91726734455
        assert Decimal.from_float(lower) < 2757902 + Decimal.from_float(close_below) < exact
        assert exact < 2757902 + Decimal.from_float(close_above) < Decimal.from_float(upper)

    for inlet, expected in (
        (window(lower), "within_meop_reference"),
        (window(upper), "meop_unit_disagreement"),
        (window(2757902, reference="gauge", datum=(close_below, close_below)), "within_meop_reference"),
        (window(2757902, reference="gauge", datum=(close_above, close_above)), "meop_unit_disagreement"),
    ):
        result = study(inlet)
        assert result.candidates[0].catalogue_pressure_checks[0].meop.status == expected


@pytest.mark.parametrize(("psia", "at_metric", "above_metric"), [
    (None, "within_meop_reference", "exceeds_meop_reference"),
    (401, "within_meop_reference", "meop_unit_disagreement"),
])
def test_metric_boundary_is_inclusive_for_a_single_limit_or_a_higher_psia_limit(psia, at_metric, above_metric):
    record = sass.thruster("moog-058e151")
    rating = record.field("maximum_expected_operating_pressure")
    alternate = None if psia is None else replace(rating.reported_psia, value=psia)
    fields = tuple(replace(field, reported_psia=alternate) if field is rating else field for field in record.fields)
    with with_moog_fields(record, fields):
        for inlet, expected in ((2760000, at_metric), (nextafter(2760000, float("inf")), above_metric)):
            result = study(window(inlet))
            assert result.candidates[0].catalogue_pressure_checks[0].meop.status == expected


def test_alternate_pressure_requires_same_row_and_absolute_inlet_scope():
    record = sass.thruster("moog-058e151")
    rating = record.field("maximum_expected_operating_pressure")
    for change in (
        {"pressure_reference": "unspecified"},
        {"pressure_port": None},
        {"reported_psia": replace(rating.reported_psia, citations=record.field("proof_pressure").upper.citations)},
    ):
        reason = None
        try:
            replace(rating, **change)
        except DomainError as error:
            reason = str(error)
        assert reason == "The alternate psia value requires a published absolute inlet ceiling from the same row."
    with pytest.raises(DomainError, match="psia"):
        replace(rating, reported_psia=replace(rating.reported_psia, unit="psig"))


def test_gauge_datum_is_used_for_meop_and_missing_reference_stays_unassessed():
    unreferenced = study(window(2800000, reference="gauge"))
    pressure = unreferenced.candidates[0].catalogue_pressure_checks[0]
    assert pressure.meop.status == "unassessed"
    converted = study(window(2670000, reference="gauge", datum=(100000, 110000)))
    pressure = converted.candidates[0].catalogue_pressure_checks[0]
    assert pressure.absolute_window_pa == (2770000, 2780000)
    assert pressure.meop.status == "exceeds_meop_reference"


def test_proof_and_burst_cannot_be_used_as_performance_or_meop_fallback():
    record = sass.thruster("moog-058e151")
    fields = tuple(field for field in record.fields
                   if field.name not in {"inlet_pressure", "maximum_expected_operating_pressure"})
    with with_moog_fields(record, fields):
        report = study(window(7000000))
    pressure = report.candidates[0].catalogue_pressure_checks[0]
    assert pressure.status == "unassessed" and pressure.meop is None and pressure.published is None
    assert report.scores[0].score is None


@pytest.mark.parametrize(("inlet", "expected"), [
    (3000000, "exceeds_meop_reference"),
    (2759000, "meop_unit_disagreement"),
])
def test_matching_nominal_does_not_hide_an_operating_rating_conflict(inlet, expected):
    record = sass.thruster("moog-058e151")
    performance = record.field("inlet_pressure")
    # Simulate a conflicting catalogue in which the nominal performance
    # reference itself exceeds the pressure rating of the selected component.
    endpoint = replace(performance.upper, value=inlet)
    fields = tuple(replace(field, lower=endpoint, upper=endpoint) if field is performance else field
                   for field in record.fields)
    with with_moog_fields(record, fields):
        result = study(window(inlet))
    pressure = result.candidates[0].catalogue_pressure_checks[0]
    assert pressure.status == "matches_nominal"
    assert pressure.meop.status == expected
    assert result.scores[0].score is None
    assert result.candidates[0].sizing is not None
    assert "Operating pressure review required" in sass.build_project_html(result)


def test_report_and_json_retain_the_operating_comparison_and_test_ratings():
    result = study(window(5000000, 6000000))
    row = result.to_dict()["candidates"][0]
    review = row["catalogue_pressure_checks"][0]["meop"]
    assert review["status"] == "exceeds_meop_reference"
    assert review["published_limit"]["reported_psia"]["value"] == 400
    assert review["published_limit"]["reported_psia"]["unit"] == "psia"
    assert {item["name"] for item in row["catalog_record"]["fields"]} >= {"proof_pressure", "burst_pressure"}
    report = sass.build_project_html(result)
    assert "Published component pressure ratings" in report
    assert "2760000 Pa (400 psia)" in report
    assert "entire inlet window exceeds published MEOP" in report

