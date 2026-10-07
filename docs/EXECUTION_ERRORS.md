# Execution errors, priors and process noise

`propulsion/burns/execution_errors.py` turns propulsion uncertainty into the
compact second-moment models that orbit determination, linear covariance
analysis and maneuver design use. The Monte Carlo ensembles
([`PROPULSION_UNCERTAINTY.md`](PROPULSION_UNCERTAINTY.md)) remain the tool for
full non-Gaussian dispersion; this module is its linear counterpart.

## Parameter prior

`PropulsionParameterPrior` holds a mean vector and a symmetric
positive-semidefinite covariance over named propagator-bridge drivers
(`thrust_scale`, `mass_flow_scale`, `ignition_time_bias_s`). Parameters not
listed are treated as known. It is hash-verified like every other artifact.

Two constructors:

- `from_standard_deviations({...}, correlation={(a, b): rho})` for engineering
  judgement or vendor data;
- `from_multivariate_calibration(calibration, scope=..., ...)` for a paired
  hot-fire campaign ([`HOT_FIRE_CALIBRATION.md`](HOT_FIRE_CALIBRATION.md)).

### Mapping from hot-fire calibration

The calibration estimates thrust, specific-impulse and optional axial-efficiency
scales \(T, I, A\). The bridge drivers are

\[
  s_F = T A,\qquad s_{\dot m} = T / I .
\]

The prior uses the JCGM 100 law of propagation of uncertainty,
\(P = J C J^{\mathsf T}\), with

\[
  J=\begin{pmatrix}A & 0 & T\\ 1/I & -T/I^2 & 0\end{pmatrix},
\]

and \(C\) built from the marginal process variances of the chosen scope
(`POPULATION` or `SAME_ENGINE`) and the calibration's effective correlation
matrix. With `artifact_includes_consensus=True` (default) the tabulated
artifact already carries the consensus correction, so the prior is the relative
residual about 1:

\[
  \sigma^2_{s_F}=r_T^2,\quad
  \sigma^2_{s_{\dot m}}=r_T^2+r_I^2-2\rho r_T r_I,\quad
  \operatorname{cov}=r_T^2-\rho r_T r_I,\qquad r_X=\sigma_X/\bar X .
\]

A 40 000-sample nonlinear Monte Carlo of \(T\) and \(T/I\) agrees with these
values to within 5 % in the tests. `time_scale` stretches burn duration and is
not a bridge driver. Ignition-time uncertainty is not observable from scale
calibration and must be supplied explicitly.

`prior.check_against(bridge.parameters)` fails on unknown names or means outside
the driver bounds, and warns when mean ± 3σ crosses a bound.

## Gates execution-error model

For a commanded impulsive \(\Delta\mathbf v = V\hat e\):

\[
  Q=(\sigma_r^2+V^2\sigma_s^2)\,\hat e\hat e^{\mathsf T}
   +(\sigma_a^2+V^2\sigma_p^2)\,(I-\hat e\hat e^{\mathsf T}),
\]

with fixed magnitude \(\sigma_r\) [m/s], proportional magnitude \(\sigma_s\),
fixed transverse \(\sigma_a\) [m/s] and proportional pointing \(\sigma_p\)
[rad]. `covariance()` returns \(Q\) and `sample()` draws an executed
\(\Delta\mathbf v\); 40 000 samples reproduce \(Q\) within 3 % in the tests.

### Proportional magnitude from propulsion uncertainty

For a duration-terminated burn with constant exhaust velocity and nominal
consumed mass fraction \(\mu\),

\[
  \frac{\partial\ln\Delta v}{\partial s_F}=1,\qquad
  \frac{\partial\ln\Delta v}{\partial s_{\dot m}}
  =-1+\frac{\mu}{(1-\mu)\,(-\ln(1-\mu))}\approx\frac{\mu}{2}.
\]

A flow-scale error largely cancels because it changes both exhaust velocity and
propellant consumed. `GatesExecutionErrorModel.from_parameter_prior` sets
\(\sigma_s=\sqrt{g^{\mathsf T}P_\mathrm{rel}g}\). The fixed and pointing terms
come from cutoff logic and attitude control and must be supplied. An
accelerometer (Δv) cutoff largely removes the proportional magnitude term;
model its resolution in \(\sigma_r\) instead.

## First-order Gauss-Markov process noise

For long low-thrust arcs, a thrust- or flow-scale deviation \(x\) can be
modelled as

\[
  \dot x=-x/\tau+w,\qquad q=2\sigma^2/\tau,
\]

with exact discrete form \(x_{k+1}=\phi x_k+w_k\), \(\phi=e^{-\Delta t/\tau}\),
\(\operatorname{Var}(w_k)=\sigma^2(1-\phi^2)\). `FirstOrderGaussMarkov` gives
the transition, discrete noise, moment propagation and exact sample paths. The
consumer adds \(x\) to its state, multiplies the bridge's `thrust_scale` by
\(1+x\), and uses the bridge's `thrust_scale` partial for \(\partial a/\partial x\).

## References

- [Navigation Filter Best Practices, NASA/TP-2018-219822, NTRS 20180003657](https://ntrs.nasa.gov/citations/20180003657)
- [JCGM 100:2008, Guide to the Expression of Uncertainty in Measurement](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf)
- C. R. Gates, *A Simplified Model of Midcourse Maneuver Execution Errors*,
  JPL Technical Report 32-504, 1963 (original model; cited bibliographically).
