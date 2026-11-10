"""Candidate red tests for cheap rejection before an expensive equilibrium solve."""

from unittest.mock import patch

import pytest

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.errors import RocketPropulsionError


@pytest.mark.parametrize(("field", "value"), [
    ("oxidizer_drop_fraction", 0.0),
    ("oxidizer_drop_fraction", -1.0),
    ("fuel_drop_fraction", 0.0),
    ("fuel_drop_fraction", -1.0),
    ("combustion_time_lag_s", -1.0),
    ("oxidizer_drop_fraction", 1e308),
    ("fuel_drop_fraction", 1e308),
])
def test_invalid_stability_request_refused_before_equilibrium(field, value):
    with (
        patch("rocket_propulsion.api.combustion_dynamics_routes._chamber_state") as state,
        patch("rocket_propulsion.api.combustion_dynamics_routes.cea_characteristic_velocity_function") as cstar,
    ):
        state.side_effect = AssertionError(f"Expensive equilibrium started for {field}={value}")
        cstar.side_effect = AssertionError(f"Expensive performance started for {field}={value}")
        with pytest.raises(RocketPropulsionError) as error:
            dispatch_calculation("/api/v1/combustion-stability", {field: value})
        assert error.value.field == field
        state.assert_not_called()
        cstar.assert_not_called()
