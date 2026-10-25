# Calibration fields over the operating-condition surface

## Purpose

`rocket_propulsion.propulsion.burns.calibration_field` turns per-condition
hot-fire calibrations ([`HOT_FIRE_CALIBRATION.md`](HOT_FIRE_CALIBRATION.md))
into a field over the rectilinear operating-condition grid used by
`RectilinearPerformanceSurface`
([`OPERATING_CONDITION_SURFACES.md`](OPERATING_CONDITION_SURFACES.md)). Each
grid node (for example one chamber pressure and one mixture ratio) is one
`condition_id` stratum with its own Mandel–Paule random-effects calibration.
Between nodes the field interpolates; outside the grid it refuses.

## Building a field

```python
field = calibrate_hot_fire_field(
    "thrust_scale",
    axes,                                   # PerformanceAxis tuple
    {"pc1-mr55": (1.0e6, 5.5), ...},        # condition_id -> grid node
    measurements,                           # HotFireScaleMeasurement tuple
    source="qualification campaign Q3",
    rationale="per-condition thrust/model ratios",
)
```

Or pass ready-made calibrations as `CalibrationFieldNode(coordinates,
calibration)` to `HotFireCalibrationField(...)`. The constructor fails closed
on a missing node, a duplicate node, a node outside an axis range or between
grid coordinates, a `condition_id` used at two nodes, or a node that calibrates
another scale. The grouping helper also refuses measurements whose condition
has no coordinates and mapped conditions without measurements.

## Interpolation rule

With the multilinear weights `w_c >= 0`, `sum_c w_c = 1`, of the corner nodes
`c` of the active cell:

```text
mu(x)        = sum_c w_c mu_c                    consensus scale
u(x)         = sum_c w_c u_c                     consensus standard uncertainty
u_indep(x)   = sqrt(sum_c w_c^2 u_c^2)           reported, uncorrelated nodes
sigma_w^2(x) = sum_c w_c sigma_w,c^2             within-engine process variance
sigma_b^2(x) = sum_c w_c sigma_b,c^2             between-engine process variance
sigma_p^2(x) = sigma_w^2(x) + sigma_b^2(x)       population process variance
```

Justification:

- `mu(x)` is a linear function of the node estimates with sensitivity
  coefficients `w_c`. JCGM 100:2008 §5.2.2 gives its combined variance for
  correlation coefficients `-1 <= r <= +1` (equations 13–16); NOTE 1 states
  that for `r = +1` the combined standard uncertainty is the linear sum
  `sum_c w_c u_c`. Because all `w_c >= 0` and `|r| <= 1`, that linear sum is
  an upper bound for every correlation between node estimates, so `u(x)` is
  a conservative envelope. The uncorrelated value (equation 10) is reported
  as `independent_nodes_standard_uncertainty` and never exceeds `u(x)`.
- Variance components are parameters, not estimates combined by a
  measurement function, so no propagation law applies. Interpolating them
  linearly in variance is this library's stated engineering convention: it
  reproduces every node, stays non-negative, keeps
  `sigma_p^2 = sigma_w^2 + sigma_b^2` exactly, and is never smaller than
  linear interpolation of standard deviations (Cauchy–Schwarz).

At a node the weight of that node is one, so the field reproduces the node
calibration exactly. Values are summed in exact rational arithmetic and
rounded once, as in the performance surface.

## Domain and events

`evaluate()` uses the performance-surface domain rule: keys must equal the
axis names, a value within `max(1e-12, 1e-12 max(|min|, |max|))` of an
endpoint snaps to it, and anything farther outside is refused.
`validity_surfaces()` returns stopping `PerformanceValiditySurface` events
`calibration-field:<scale>:<axis>:minimum|maximum` (positive inside), and
`grid_line_surfaces()` returns the non-stopping kink events of the
multilinear interpolant.

## Composition

`apply_calibration_fields(surface, thrust_field=..., specific_impulse_field=...)`
requires fields on exactly the surface axes and returns a
`CalibratedPerformanceSurface`:

- node force scale `s_F mu_F` and node flow scale `s_mdot mu_F / mu_I`
  (mass flow is `F / (g0 Isp)`), interpolated with the base surface's own
  rule; the source hash binds the base surface hash and the field hashes;
- `evaluate()` adds the field evaluations and relative standard uncertainties
  `u_F/mu_F`, `u_I/mu_I`, and the first-order flow bound
  `u_F/mu_F + u_I/mu_I` (JCGM 100:2008 §5.1.6 and §5.2.2 NOTE 4 with
  `|r| <= 1`, since the two fields usually come from the same runs).

For the ensemble, `field.ensemble_parameters(conditions)` returns two
`UncertainParameter` factors centred at one, to be used after the
interpolated consensus has been applied deterministically:

- `<scale>::aleatory:hot-fire-field`, relative standard deviation
  `sigma_p(x)/mu(x)` (or `sigma_w(x)/mu(x)` for `CalibrationScope.SAME_ENGINE`);
- `<scale>::epistemic:hot-fire-field-consensus`, relative standard deviation
  `u(x)/mu(x)`.

Both use the existing truncation `[max(1e-12, 1 - k r), 1 + k r]`, `k = 3` by
default, and are accepted by `ScenarioConditionalUncertainty`.
`residual_uncertain_parameter()` alone, with its plain scale name, matches
`HotFireScaleCalibration.residual_uncertain_parameter()` at every node.

## Persistence

`rocket_propulsion_hot_fire_calibration_field_v1` stores the axes, every node
calibration with its measurements and hash, the references, the numerical
policies `interpolation_numerics` and `uncertainty_numerics`, and
`field_hash`. Reopening reruns every node calibration from its measurements
and rejects hash changes, falsified derived results and changed policies.

## Credibility boundary

- No interpolation (model-form) error between nodes is estimated; dense
  enough grids or separate discrepancy evidence
  ([`MODEL_DISCREPANCY_TRANSFER.md`](MODEL_DISCREPANCY_TRANSFER.md)) are
  required to cover it.
- Variance-component sampling uncertainty is not included, and the
  fewer-than-six-engines warning of each node is carried over.
- Covariance between nodes is bounded, not estimated; the thrust and Isp
  ensemble factors are independent unless a correlation model is supplied.
- Only multilinear interpolation is offered for the field, because the
  envelope bound needs non-negative weights.

## Primary references

- JCGM 100:2008, *Guide to the Expression of Uncertainty in Measurement*,
  §5.1.2, §5.1.6, §5.2.2:
  <https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf>
- NIST/SEMATECH e-Handbook §7.4.4, *Variance Components*:
  <https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm>
- NIST Dataplot, *Consensus Mean*:
  <https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm>
- NASA Glenn, *Specific Impulse*:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/>
- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
