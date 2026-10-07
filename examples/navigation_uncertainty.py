"""Turn propulsion uncertainty into navigation inputs.

1. A prior over the propagator-bridge drivers from engineering sigmas.
2. The Gates execution-error covariance of a 25 m/s burn, with the
   proportional magnitude term derived from that prior.
3. A Gauss-Markov thrust-scale process for a long low-thrust arc.

All numbers are illustrative.
"""

import random

from rocket_propulsion.propulsion.burns import (
    FirstOrderGaussMarkov,
    GatesExecutionErrorModel,
    PropulsionParameterPrior,
)

prior = PropulsionParameterPrior.from_standard_deviations(
    {"thrust_scale": 0.015, "mass_flow_scale": 0.010, "ignition_time_bias_s": 0.05},
    correlation={("thrust_scale", "mass_flow_scale"): 0.6},
    source="illustrative acceptance-test sigmas",
    rationale="example only",
)
print("prior sigmas:", dict(zip(prior.parameter_names, prior.standard_deviations, strict=True)))
print("prior hash:", prior.prior_hash)

gates = GatesExecutionErrorModel.from_parameter_prior(
    prior,
    consumed_mass_fraction=0.02,
    fixed_magnitude_sigma_m_s=0.002,
    fixed_pointing_sigma_m_s=0.001,
    proportional_pointing_sigma_rad=0.003,
)
delta_v = (0.0, 25.0, 0.0)
print("proportional magnitude sigma:", round(gates.proportional_magnitude_sigma, 6))
for row in gates.covariance(delta_v):
    print("Q row [m^2/s^2]:", tuple(round(value, 8) for value in row))
print("one executed delta-v:", gates.sample(delta_v, random.Random(1)))

process = FirstOrderGaussMarkov("thrust_scale", steady_state_sigma=0.01, time_constant_s=6 * 3600)
print("phi over 1 h:", round(process.transition(3600.0), 6))
print("discrete noise variance over 1 h:", process.discrete_noise_variance(3600.0))
