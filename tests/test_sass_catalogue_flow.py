"""Attack catalogue flow-reference provenance, comparison and trade eligibility."""

from dataclasses import replace
from fractions import Fraction
from math import nextafter
from unittest.mock import patch

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies import sass
from rocket_propulsion.studies.sass import catalog as catalogue
from tests.test_sass_studies import q


def study(*keys):
    return sass.analyze_project(
        "Synthetic catalogue comparison", q(500, "kg"),
        (sass.DeltaVItem("short burn", q(5, "m/s"), q(0), 1),),
        tuple(sass.Candidate(key, key, q(3, "kg"), 1, 1, q(1), q(1)) for key in keys),
        residual_fraction=q(0), fixed_reserve=q(0, "kg"), requirements=(),
        criteria=(sass.Criterion("wet_mass", q(500, "kg"), q(600, "kg"), q(1)),),
    )


def with_field(record, field):
    changed = replace(record, fields=tuple(field if item.name == field.name else item
                                           for item in record.fields))
    return patch.object(catalogue, "THRUSTERS", tuple(changed if item is record else item
                                                      for item in catalogue.THRUSTERS))


def test_manufacturer_1n_conflict_is_retained_but_cannot_win_trade():
    report = study("ariane-1n", "ariane-s10-13")
    one, ten = report.candidates
    assert one.catalogue_flow_check.status == "above_published_range"
    assert ten.catalogue_flow_check.status == "within_published_range"
    assert one.sizing is not None
    assert one.catalogue_flow_check.reported_range.upper.value == .000447
    assert one.operating_point.get("reported_mass_flow").value == .00044
    assert one.operating_point.get("closed_mass_flow").value == float(
        Fraction(1) / (Fraction(9.80665) * 220)
    )
    assert [(row.eligibility, row.score is None) for row in report.scores] == [
        ("incomplete", True), ("eligible", False),
    ]
    entry = report.to_dict()["candidates"][0]["catalogue_flow_check"]
    assert entry["status"] == "above_published_range"
    assert "ARIANE-1N" in {row["source_id"] for row in entry["reported_range"]["upper"]["citations"]}
    rendered = sass.build_project_html(report)
    assert "Published mass-flow range" in rendered
    assert "Above published range" in rendered
    assert "0.000447" in rendered


def test_flow_comparison_does_not_assume_endpoint_from_different_sources():
    record = sass.thruster("ariane-1n")
    bounds = record.field("reported_mass_flow_range")
    for field in (
        replace(bounds, lower=replace(bounds.lower, citations=(
            sass.Citation("ARIANE-10N", "Synthetic incompatible lower endpoint"),))),
        replace(bounds, upper=replace(bounds.upper, citations=(
            sass.Citation("ARIANE-10N", "Synthetic incompatible upper endpoint"),))),
        replace(bounds, kind="upper_bound", lower=None),
    ):
        with with_field(record, field):
            result = study(record.key)
        point = result.candidates[0]
        assert point.catalogue_flow_check.status == "unassessed"
        assert point.operating_point.get("reported_mass_flow").value == .00044
        assert result.scores[0].eligibility == "eligible"


def test_below_range_also_blocks_ranking_without_removing_budget():
    record = sass.thruster("ariane-1n")
    bounds = record.field("reported_mass_flow_range")
    synthetic = replace(bounds, lower=replace(bounds.lower, value=.0005),
                         upper=replace(bounds.upper, value=.0006))
    with with_field(record, synthetic):
        report = study(record.key)
    assert report.candidates[0].catalogue_flow_check.status == "below_published_range"
    assert report.candidates[0].sizing is not None
    assert report.scores[0].score is None


def test_exact_flow_outside_range_cannot_round_itself_inside():
    record = sass.thruster("ariane-1n")
    bounds = record.field("reported_mass_flow_range")
    flow = Fraction(1) / (Fraction(9.80665) * 220)
    rounded = float(flow)
    if Fraction(rounded) >= flow:
        endpoint, direction = "lower", "below_published_range"
        assert Fraction(rounded) > flow
    else:
        endpoint, direction = "upper", "above_published_range"
    other = "upper" if endpoint == "lower" else "lower"
    safe = 1.0 if endpoint == "lower" else .0001
    for value, expected in ((rounded, direction),
                            (nextafter(rounded, -1 if endpoint == "lower" else 1),
                             "within_published_range")):
        field = replace(bounds, **{endpoint: replace(getattr(bounds, endpoint), value=value),
                                   other: replace(getattr(bounds, other), value=safe)})
        with with_field(record, field):
            report = study(record.key)
        assert report.candidates[0].catalogue_flow_check.status == expected


def test_same_source_range_rejects_wrong_units_and_nonpositive_values():
    record = sass.thruster("ariane-1n")
    bounds = record.field("reported_mass_flow_range")
    for changes in (
        {"lower": replace(bounds.lower, unit="g/s"), "upper": replace(bounds.upper, unit="g/s")},
        {"lower": replace(bounds.lower, value=0)},
    ):
        with with_field(record, replace(bounds, **changes)), pytest.raises(DomainError):
            study(record.key)

