"""Normal-shock oracle built from conservation of mass, momentum and enthalpy."""

import json
import math
from decimal import ROUND_DOWN, Inexact, localcontext

import mpmath as mp
import pytest

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.compressible.normal_shock import calculate_normal_shock
from rocket_propulsion.core.errors import DomainError


def conservation_reference(mach, gamma):
    with mp.workdps(520):
        m, g = mp.mpf(mach), mp.mpf(gamma)
        # Let volume be rho1/rho2, which is also the velocity ratio by mass
        # conservation. Substituting the momentum and ideal-gas relations
        # in enthalpy conservation gives a quadratic. Factoring out the
        # trivial no-shock root yields the compressive volume ratio below.
        volume = ((g - 1) + 2 / m**2) / (g + 1)
        pressure = 1 + g * m**2 * (1 - volume)
        temperature = pressure * volume
        # Ideal-gas entropy difference, divided by R, with constant total T.
        entropy = g / (g - 1) * mp.log(temperature) - mp.log(pressure)
        return {
            "downstream_mach": float(m * volume / mp.sqrt(temperature)),
            "pressure_ratio": float(pressure),
            "density_ratio": float(1 / volume),
            "temperature_ratio": float(temperature),
            "stagnation_pressure_ratio": float(mp.exp(-entropy)),
        }


@pytest.mark.parametrize(("mach", "gamma"), [
    (math.nextafter(1.0, math.inf), math.nextafter(1.0, math.inf)),
    (math.nextafter(1.0, math.inf), 1.4),
    (1.001, 1.4),
    (2.0, math.nextafter(1.0, math.inf)),
    (10.0, 1.001),
    (2.0, 1.4),
    (50.0, 1.4),
    (2.0, 1e308),
    (1e100, 1e308),
    (1e60, 1.4),
    (1e64, 1.4),
    (1e153, 10.0),
])
def test_normal_shock_complete_outputs_against_conservation(mach, gamma):
    expected = conservation_reference(mach, gamma)
    assert all(math.isfinite(value) and value > 0 for value in expected.values())
    try:
        result = calculate_normal_shock(mach, gamma)
    except (DomainError, ArithmeticError) as exc:
        pytest.fail(f"All independent outputs are representable for {mach=}, {gamma=}: {exc!r}")
    assert 0 < result.stagnation_pressure_ratio <= 1
    for field, value in expected.items():
        actual = getattr(result, field)
        assert math.isfinite(actual)
        assert abs(actual - value) <= 2 * math.ulp(value), field


@pytest.mark.parametrize(("mach", "gamma", "field"), [
    (1e155, 10.0, "pressure_ratio"),
    (50.0, 1.001, "stagnation_pressure_ratio"),
    (1e66, 1.4, "stagnation_pressure_ratio"),
])
def test_normal_shock_unrepresentable_outputs_have_typed_failures(mach, gamma, field):
    independent = conservation_reference(mach, gamma)
    assert independent[field] == 0 or math.isinf(independent[field])
    try:
        dispatch_calculation("/api/v1/normal-shock", {"upstream_mach": mach, "gamma": gamma})
    except DomainError as error:
        assert error.code == "numeric_range"
        assert error.field == field
    else:
        raise AssertionError(f"Unrepresentable {field} escaped")


@pytest.mark.parametrize(("field", "value"), [
    ("gamma", True), ("gamma", 1.0), ("gamma", math.inf),
    ("upstream_mach", "2.0"), ("upstream_mach", 1.0), ("upstream_mach", math.nan),
])
def test_normal_shock_direct_input_failure_identifies_the_field(field, value):
    inputs = {"upstream_mach": 2.0, "gamma": 1.4} | {field: value}
    try:
        calculate_normal_shock(**inputs)
    except DomainError as error:
        assert error.field == field
    else:
        raise AssertionError(f"Invalid input accepted for {field}")


def test_normal_shock_route_retains_context_and_strict_json():
    mach, gamma = 2.0, math.nextafter(1.0, math.inf)
    expected = conservation_reference(mach, gamma)
    with localcontext() as ctx:
        ctx.prec = 5
        ctx.rounding = ROUND_DOWN
        ctx.Emin, ctx.Emax = -9, 9
        ctx.traps[Inexact] = True
        flags = dict(ctx.flags)
        try:
            answer = dispatch_calculation("/api/v1/normal-shock", {"upstream_mach": mach, "gamma": gamma})
        except (DomainError, ArithmeticError) as exc:
            raise AssertionError(f"Representable shock failed under caller context: {exc!r}") from exc
        assert dict(ctx.flags) == flags
        for key, value in expected.items():
            assert abs(answer[key] - value) <= 2 * math.ulp(value)
        json.dumps(answer, allow_nan=False)
