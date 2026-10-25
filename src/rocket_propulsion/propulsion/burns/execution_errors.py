"""Navigation-facing propulsion uncertainty: priors, Gates errors, process noise.

Monte Carlo ensembles (``uncertainty.py``, ``ensemble.py``) carry complete
realizations.  Orbit determination, linear covariance analysis and maneuver
design instead need compact second-moment models:

``PropulsionParameterPrior``
    Mean and covariance of the propagator-bridge drivers (``thrust_scale``,
    ``mass_flow_scale``, ``ignition_time_bias_s``), built from explicit
    standard deviations or linearized from a paired hot-fire calibration.
    It is the a-priori / consider covariance for those drivers.
``GatesExecutionErrorModel``
    The classical four-parameter impulsive maneuver execution-error model:
    fixed and proportional magnitude errors along the burn, fixed and
    proportional pointing errors across it.  The proportional magnitude term
    can be derived from a parameter prior.
``FirstOrderGaussMarkov``
    Exponentially correlated (first-order Gauss-Markov) process noise for a
    slowly wandering thrust or flow scale during long low-thrust arcs.

Pointing statistics come from attitude control and remain the consumer's
input; this module only supplies the propulsion-owned terms. Prior/consider
covariance, Gauss-Markov process noise and execution-error models follow [1]_;
the execution-error model originates in [2]_, and covariance mapping uses the
first-order law of propagation of uncertainty [3]_ [4]_.

References
----------
.. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
       https://ntrs.nasa.gov/citations/20180003657
.. [2] C. R. Gates, "A Simplified Model of Midcourse Maneuver Execution
       Errors", JPL Technical Report 32-504, 1963.
.. [3] JCGM 100:2008, "Evaluation of Measurement Data - Guide to the Expression
       of Uncertainty in Measurement".
       https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
.. [4] NIST-TN1297-COVARIANCE. NIST TN 1297 Appendix A, first-order
       measurement uncertainty with sensitivity coefficients and covariance.
"""

from __future__ import annotations

import random
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise
from math import exp, isfinite, log1p, sqrt

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    real,
    real_matrix,
    real_tuple,
    records,
    strings,
    text,
    vector3,
)

from .calibration import CalibrationScope
from .multivariate_calibration import HotFireMultivariateCalibration
from .propagator import EstimableParameter
from .serialization import canonical_sha256
from .vehicle_geometry import _cross

PROPULSION_PRIOR_SCHEMA = "rocket_propulsion_parameter_prior_v1"
EXECUTION_ERROR_REFERENCE_IDS = ("NASA-20180003657", "JCGM-100-2008")
_SYMMETRY_TOLERANCE = 1e-12

Matrix = tuple[tuple[float, ...], ...]




def _cholesky(matrix: Matrix) -> Matrix | None:
    """Return a lower Cholesky factor, or ``None`` if not positive semidefinite.

    A diagonal residual down to ``-tol`` is treated as zero, with
    ``tol = max(1e-300, 1e-12 max_i |a_ii|)``; a zero pivot then requires the
    rest of its column to vanish within ``tol``. The factor of a singular
    positive semidefinite matrix therefore has zero columns.
    """

    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    scale = max((abs(matrix[i][i]) for i in range(size)), default=0.0)
    tolerance = max(1e-300, scale * 1e-12)
    for row in range(size):
        for column in range(row + 1):
            total = matrix[row][column] - sum(
                lower[row][k] * lower[column][k] for k in range(column)
            )
            if row == column:
                if total < -tolerance:
                    return None
                lower[row][column] = sqrt(max(total, 0.0))
            elif lower[column][column] > 0.0:
                lower[row][column] = total / lower[column][column]
            elif abs(total) > tolerance:
                return None
    return tuple(tuple(row) for row in lower)


@dataclass(frozen=True, slots=True)
class PropulsionParameterPrior:
    """Gaussian prior (mean and covariance) over named bridge drivers.

    The prior serves as the a-priori or consider covariance of the drivers in
    orbit determination and linear covariance analysis [1]_. Parameters not
    listed are treated as known.

    Parameters
    ----------
    parameter_names : tuple of str
        Unique, nonempty driver names.
    mean : tuple of float
        Finite prior means in the drivers' own units (dimensionless scales,
        seconds for timing bias).
    covariance : tuple of tuple of float
        Finite, symmetric (relative tolerance ``1e-12``) positive semidefinite
        covariance in the same units [2]_.
    source, rationale : str
        Nonempty provenance and justification.
    reference_ids : tuple of str, optional
        Nonempty, unique source identifiers.
    schema : str, optional
        Fixed schema identifier.

    Raises
    ------
    DomainError
        If any constraint above is violated.

    References
    ----------
    .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
           https://ntrs.nasa.gov/citations/20180003657
    .. [2] JCGM 100:2008, section 5.2, "Evaluation of Measurement Data -
           Guide to the Expression of Uncertainty in Measurement".
           https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    """

    parameter_names: tuple[str, ...]
    mean: tuple[float, ...]
    covariance: Matrix
    source: str
    rationale: str
    reference_ids: tuple[str, ...] = EXECUTION_ERROR_REFERENCE_IDS
    schema: str = PROPULSION_PRIOR_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != PROPULSION_PRIOR_SCHEMA:
            raise DomainError("Unsupported propulsion-prior schema.")
        for field_name in ("parameter_names", "reference_ids"):
            object.__setattr__(
                self, field_name, strings(getattr(self, field_name), f"Prior {field_name}")
            )
        object.__setattr__(self, "mean", real_tuple(self.mean, "Prior mean"))
        object.__setattr__(self, "covariance", real_matrix(self.covariance, "Prior covariance"))
        if not isinstance(self.source, str) or not isinstance(self.rationale, str):
            raise DomainError("Prior source and rationale must be strings.")
        size = len(self.parameter_names)
        if size == 0:
            raise DomainError("A parameter prior requires at least one parameter.")
        if len(set(self.parameter_names)) != size or any(
            not name.strip() for name in self.parameter_names
        ):
            raise DomainError("Prior parameter names must be unique and non-empty.")
        if len(self.mean) != size or len(self.covariance) != size or any(
            len(row) != size for row in self.covariance
        ):
            raise DomainError("Prior mean and covariance dimensions must match the names.")
        values = (*self.mean, *(value for row in self.covariance for value in row))
        if any(not isfinite(value) for value in values):
            raise DomainError("Prior mean and covariance must be finite.")
        for i in range(size):
            for j in range(i):
                tolerance = _SYMMETRY_TOLERANCE * max(
                    1.0, abs(self.covariance[i][j]), abs(self.covariance[j][i])
                )
                if abs(self.covariance[i][j] - self.covariance[j][i]) > tolerance:
                    raise DomainError("Prior covariance must be symmetric.")
        if _cholesky(self.covariance) is None:
            raise DomainError("Prior covariance must be positive semidefinite.")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Prior source and rationale must not be empty.")
        if not self.reference_ids or len(set(self.reference_ids)) != len(self.reference_ids):
            raise DomainError("Prior references must be non-empty and unique.")

    @classmethod
    def from_standard_deviations(
        cls,
        standard_deviations: Mapping[str, float],
        *,
        source: str,
        rationale: str,
        mean: Mapping[str, float] | None = None,
        correlation: Mapping[tuple[str, str], float] | None = None,
    ) -> PropulsionParameterPrior:
        """Build a prior from standard deviations and pairwise correlations.

        The covariance is ``C_ij = r_ij s_i s_j`` [1]_ with ``r_ii = 1`` and
        unspecified pairs uncorrelated.

        Parameters
        ----------
        standard_deviations : mapping of str to float
            Finite, nonnegative sigma of each driver; defines the order.
        source, rationale : str
            Nonempty provenance and justification.
        mean : mapping of str to float, optional
            Means of listed drivers. Missing means default to the bridge
            reference values (1 for ``thrust_scale`` and ``mass_flow_scale``,
            0 for ``ignition_time_bias_s``); any other name must be given.
        correlation : mapping of (str, str) to float, optional
            Coefficients in ``[-1, 1]`` for pairs of distinct listed drivers;
            applied symmetrically.

        Returns
        -------
        PropulsionParameterPrior
            The constructed prior.

        Raises
        ------
        DomainError
            If an argument has the wrong type or value, or the resulting
            covariance is not positive semidefinite.

        References
        ----------
        .. [1] JCGM 100:2008, section 5.2, "Evaluation of Measurement Data -
               Guide to the Expression of Uncertainty in Measurement".
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        """

        if not isinstance(standard_deviations, Mapping):
            raise DomainError("Standard deviations must be a mapping.")
        if mean is not None and not isinstance(mean, Mapping):
            raise DomainError("Means must be a mapping.")
        if correlation is not None and not isinstance(correlation, Mapping):
            raise DomainError("Correlations must be a mapping of name pairs.")
        names = tuple(text(name, "Prior parameter name") for name in standard_deviations)
        defaults = {"thrust_scale": 1.0, "mass_flow_scale": 1.0, "ignition_time_bias_s": 0.0}
        supplied_mean = dict(mean or {})
        if set(supplied_mean) - set(names):
            raise DomainError("Means were supplied for parameters without a sigma.")
        means = []
        for name in names:
            if name in supplied_mean:
                means.append(real(supplied_mean[name], f"Prior mean {name!r}"))
            elif name in defaults:
                means.append(defaults[name])
            else:
                raise DomainError(f"Parameter {name!r} needs an explicit mean.")
        sigmas = []
        for name in names:
            sigma = real(standard_deviations[name], f"Prior sigma {name!r}")
            if sigma < 0.0:
                raise DomainError("Prior standard deviations must be finite and non-negative.")
            sigmas.append(sigma)
        rho = {name: {other: (1.0 if name == other else 0.0) for other in names} for name in names}
        for pair, raw in (correlation or {}).items():
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise DomainError("Correlation keys must be (name, name) pairs.")
            first, second = pair
            if first not in rho or second not in rho or first == second:
                raise DomainError("Correlation pairs must name two distinct prior parameters.")
            value = real(raw, "Correlation coefficient")
            if not -1.0 <= value <= 1.0:
                raise DomainError("Correlation coefficients must lie in [-1, 1].")
            rho[first][second] = rho[second][first] = float(value)
        covariance = tuple(
            tuple(rho[a][b] * sigmas[i] * sigmas[j] for j, b in enumerate(names))
            for i, a in enumerate(names)
        )
        return cls(names, tuple(means), covariance, source, rationale)

    @classmethod
    def from_multivariate_calibration(
        cls,
        calibration: HotFireMultivariateCalibration,
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        artifact_includes_consensus: bool = True,
        ignition_time_bias_sigma_s: float | None = None,
    ) -> PropulsionParameterPrior:
        """Linearize a paired hot-fire calibration onto the bridge drivers.

        With calibrated scales ``T`` (thrust), ``I`` (specific impulse) and
        optional ``A`` (axial efficiency), the bridge drivers are
        ``thrust_scale = T A`` and ``mass_flow_scale = T / I``.  First-order
        propagation (law of propagation of uncertainty [1]_ [2]_) gives
        ``P = J C J^T`` with Jacobian ``J`` of the mapping at the consensus
        means.  ``C`` uses the marginal process variances of the
        selected ``scope`` and the calibration's effective (possibly
        regularized) correlation, exactly as ``residual_model`` does.

        If ``artifact_includes_consensus`` (default) the artifact was built
        with the consensus corrections already applied, so the prior is the
        *relative* residual about 1.  Otherwise the mean is the consensus
        product and the covariance is absolute.

        ``time_scale`` stretches burn duration and is not a bridge driver; it
        is ignored.  The ignition-time sigma is not observable from scale
        calibration and must be supplied explicitly to be included.

        Parameters
        ----------
        calibration : HotFireMultivariateCalibration
            Paired calibration containing ``thrust_scale`` and
            ``specific_impulse_scale``.
        scope : CalibrationScope, optional
            Population or same-engine variances and correlation.
        artifact_includes_consensus : bool, optional
            Return the relative residual about one (default) or the absolute
            consensus mean and covariance.
        ignition_time_bias_sigma_s : float or None, optional
            Nonnegative ignition-time sigma [s]; adds an uncorrelated
            ``ignition_time_bias_s`` driver.

        Returns
        -------
        PropulsionParameterPrior
            Prior over ``thrust_scale``, ``mass_flow_scale`` and, optionally,
            ``ignition_time_bias_s``, for use as an a-priori or consider
            covariance [3]_.

        Raises
        ------
        DomainError
            If the calibration or scope is invalid, a required scale is
            missing, the same-engine correlation is undefined, or the sigma is
            negative.

        Notes
        -----
        The mapping is first order [2]_; it holds only while the relative
        variances are small enough for the linearization of ``T A`` and
        ``T / I`` about the consensus means to be adequate.

        References
        ----------
        .. [1] JCGM 100:2008, section 5.1.2, "Evaluation of Measurement Data -
               Guide to the Expression of Uncertainty in Measurement".
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        .. [2] NIST-TN1297-COVARIANCE. NIST TN 1297 Appendix A, first-order
               measurement uncertainty with sensitivity coefficients and covariance.
        .. [3] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        if not isinstance(calibration, HotFireMultivariateCalibration):
            raise DomainError("A HotFireMultivariateCalibration is required.")
        names = calibration.parameter_names
        if "thrust_scale" not in names or "specific_impulse_scale" not in names:
            raise DomainError(
                "Calibration must contain thrust_scale and specific_impulse_scale."
            )
        if not isinstance(scope, CalibrationScope):
            raise DomainError("Calibration scope must be a CalibrationScope value.")
        evidence = (
            calibration.population_correlation
            if scope is CalibrationScope.POPULATION
            else calibration.same_engine_correlation
        )
        if evidence is None:
            raise DomainError("Same-engine correlation is undefined for this calibration.")
        index = {name: position for position, name in enumerate(names)}
        means = [marginal.consensus_scale for marginal in calibration.marginals]
        variances = [
            marginal.population_process_variance
            if scope is CalibrationScope.POPULATION
            else marginal.within_engine_process_variance
            for marginal in calibration.marginals
        ]
        sigmas = [sqrt(max(value, 0.0)) for value in variances]
        correlation = evidence.effective_matrix
        selected = [index["thrust_scale"], index["specific_impulse_scale"]]
        if "axial_efficiency_scale" in index:
            selected.append(index["axial_efficiency_scale"])
        covariance_in = [
            [correlation[a][b] * sigmas[a] * sigmas[b] for b in selected] for a in selected
        ]
        thrust = means[selected[0]]
        isp = means[selected[1]]
        axial = means[selected[2]] if len(selected) == 3 else 1.0
        jacobian = [[axial, 0.0], [1.0 / isp, -thrust / isp**2]]
        if len(selected) == 3:
            jacobian[0].append(thrust)
            jacobian[1].append(0.0)
        mapped_mean = [thrust * axial, thrust / isp]
        mapped = [
            [
                sum(
                    jacobian[r][a] * covariance_in[a][b] * jacobian[c][b]
                    for a in range(len(selected))
                    for b in range(len(selected))
                )
                for c in range(2)
            ]
            for r in range(2)
        ]
        if artifact_includes_consensus:
            mapped = [
                [mapped[r][c] / (mapped_mean[r] * mapped_mean[c]) for c in range(2)]
                for r in range(2)
            ]
            mapped_mean = [1.0, 1.0]
        names_out: list[str] = ["thrust_scale", "mass_flow_scale"]
        if ignition_time_bias_sigma_s is not None:
            ignition_time_bias_sigma_s = real(ignition_time_bias_sigma_s, "Ignition-time sigma")
            if ignition_time_bias_sigma_s < 0.0:
                raise DomainError("Ignition-time sigma must be finite and non-negative.")
            names_out.append("ignition_time_bias_s")
            mapped_mean.append(0.0)
            mapped = [[*row, 0.0] for row in mapped]
            mapped.append([0.0, 0.0, ignition_time_bias_sigma_s**2])
        symmetric = tuple(
            tuple(0.5 * (mapped[r][c] + mapped[c][r]) for c in range(len(mapped)))
            for r in range(len(mapped))
        )
        return cls(
            tuple(names_out),
            tuple(mapped_mean),
            symmetric,
            source=(
                f"{calibration.source}; paired hot-fire calibration "
                f"sha256={calibration.calibration_hash}"
            ),
            rationale=(
                f"{scope.value} scope, first-order mapping T*A and T/I"
                + (", relative to consensus" if artifact_includes_consensus else ", absolute")
            ),
        )

    @property
    def standard_deviations(self) -> tuple[float, ...]:
        """Return the prior standard deviations ``sqrt(C_ii)`` in name order."""

        return tuple(sqrt(max(self.covariance[i][i], 0.0)) for i in range(len(self.mean)))

    def relative_covariance(self) -> Matrix:
        """Return the relative covariance ``C_ij / (mean_i mean_j)``.

        Returns
        -------
        tuple of tuple of float
            Dimensionless covariance of the relative deviations.

        Raises
        ------
        DomainError
            If any mean is zero.
        """

        if any(value == 0.0 for value in self.mean):
            raise DomainError("Relative covariance requires non-zero means.")
        return tuple(
            tuple(
                self.covariance[i][j] / (self.mean[i] * self.mean[j])
                for j in range(len(self.mean))
            )
            for i in range(len(self.mean))
        )

    def variance_of(self, name: str) -> float:
        """Return the variance of ``name``, or zero for a driver not in the prior.

        A driver absent from the prior is treated as known.
        """

        if name not in self.parameter_names:
            return 0.0
        position = self.parameter_names.index(name)
        return self.covariance[position][position]

    def check_against(self, parameters: tuple[EstimableParameter, ...]) -> tuple[str, ...]:
        """Validate the prior against bridge drivers and flag 3-sigma bound overlap.

        Parameters
        ----------
        parameters : tuple of EstimableParameter
            Drivers of the propagator bridge.

        Returns
        -------
        tuple of str
            One warning per driver whose ``mean -/+ 3 sigma`` interval leaves
            the driver bounds; such a prior can push an estimator to a bound
            [1]_.

        Raises
        ------
        DomainError
            If a prior name is not a driver or a mean violates its driver
            bounds.

        References
        ----------
        .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        parameters = records(parameters, EstimableParameter, "Bridge parameters")
        drivers = {parameter.name: parameter for parameter in parameters}
        warnings: list[str] = []
        for name, value, sigma in zip(
            self.parameter_names, self.mean, self.standard_deviations, strict=True
        ):
            if name not in drivers:
                raise DomainError(f"Prior parameter {name!r} is not a driver of this model.")
            driver = drivers[name]
            driver.validate(value)
            if value - 3.0 * sigma < driver.minimum or value + 3.0 * sigma > driver.maximum:
                warnings.append(
                    f"{name}: mean +/- 3 sigma exceeds driver bounds "
                    f"[{driver.minimum}, {driver.maximum}]."
                )
        return tuple(warnings)

    @property
    def prior_hash(self) -> str:
        """SHA-256 hex digest of the canonical prior content without its hash.

        A reproducible content identity of names, moments and provenance, not
        an authentication.
        """

        return canonical_sha256(self.to_dict(include_hash=False))

    def to_dict(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic JSON-compatible prior document.

        Parameters
        ----------
        include_hash : bool, optional
            Add ``prior_hash`` computed over the document without it.

        Returns
        -------
        dict
            Schema, names, mean, covariance, provenance and references.
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "parameter_names": list(self.parameter_names),
            "mean": list(self.mean),
            "covariance": [list(row) for row in self.covariance],
            "source": self.source,
            "rationale": self.rationale,
            "reference_ids": list(self.reference_ids),
        }
        if include_hash:
            payload["prior_hash"] = canonical_sha256(payload)
        return payload

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PropulsionParameterPrior:
        """Reopen a prior document and verify its content hash.

        Parameters
        ----------
        payload : mapping
            Decoded JSON document [1]_ [2]_ produced by :meth:`to_dict`.

        Returns
        -------
        PropulsionParameterPrior
            The validated prior.

        Raises
        ------
        InputError
            If the document is malformed, violates the prior constraints, or
            its ``prior_hash`` does not match the content.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
               duplicate-member and nonstandard numeric token behavior.
        """

        if not isinstance(payload, Mapping):
            raise InputError("Parameter-prior artifact root must be an object.")
        try:
            prior = cls(
                parameter_names=payload["parameter_names"],
                mean=payload["mean"],
                covariance=payload["covariance"],
                source=payload["source"],
                rationale=payload["rationale"],
                reference_ids=payload["reference_ids"],
                schema=payload["schema"],
            )
        except (KeyError, TypeError, ValueError, DomainError) as error:
            raise InputError("Malformed parameter-prior artifact.") from error
        if payload.get("prior_hash") != prior.prior_hash:
            raise InputError("Parameter-prior hash does not match its contents.")
        return prior


def delta_v_log_sensitivities(consumed_mass_fraction: float) -> tuple[float, float]:
    """Return ``(dlnDv/ds_T, dlnDv/ds_m)`` for a duration-terminated burn.

    From the ideal rocket equation [1]_, for constant effective exhaust
    velocity and a fixed burn duration,
    ``Dv = (s_T / s_m) c ln(1 / (1 - s_m mu))`` where ``mu`` is the nominal
    consumed fraction of initial mass.  At ``s_T = s_m = 1``::

        dlnDv/ds_T = 1
        dlnDv/ds_m = -1 + mu / ((1 - mu) (-ln(1 - mu)))   (~ mu/2 for small mu)

    A flow-scale error mostly cancels: it changes both the exhaust velocity
    and the propellant consumed. The sensitivities feed the proportional
    magnitude term of the execution-error model [2]_.

    Parameters
    ----------
    consumed_mass_fraction : float
        Nominal consumed fraction ``mu`` of initial mass, in ``(0, 1)``.

    Returns
    -------
    tuple of float
        Logarithmic sensitivities to the thrust and mass-flow scales.

    Raises
    ------
    DomainError
        If ``mu`` is not in ``(0, 1)``.

    References
    ----------
    .. [1] NASA Glenn, "Ideal Rocket Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] NASA/TP-2018-219822, "Navigation Filter Best Practices".
           https://ntrs.nasa.gov/citations/20180003657
    """

    mu = real(consumed_mass_fraction, "Consumed mass fraction")
    if not 0.0 < mu < 1.0:
        raise DomainError("Consumed mass fraction must lie in (0, 1).")
    return 1.0, -1.0 + mu / ((1.0 - mu) * -log1p(-mu))


@dataclass(frozen=True, slots=True)
class GatesExecutionErrorModel:
    """Four-parameter impulsive maneuver execution-error covariance.

    For a commanded velocity change ``Dv = V e`` (unit ``e``)::

        Q = (sr^2 + V^2 ss^2) e e^T + (sa^2 + V^2 sp^2) (I - e e^T)

    ``sr`` fixed magnitude (resolution/shutoff) error [m/s], ``ss``
    proportional magnitude error [-], ``sa`` fixed transverse error [m/s] and
    ``sp`` proportional pointing error [rad], all one-sigma.  The transverse
    term applies independently to each of the two axes perpendicular to
    ``e``.  The model originates in [1]_ and is used for maneuver execution
    errors in navigation filters [2]_.

    Parameters
    ----------
    fixed_magnitude_sigma_m_s : float
        ``sr`` [m/s], finite and nonnegative.
    proportional_magnitude_sigma : float
        ``ss`` [-], finite and nonnegative.
    fixed_pointing_sigma_m_s : float
        ``sa`` [m/s], finite and nonnegative.
    proportional_pointing_sigma_rad : float
        ``sp`` [rad], finite and nonnegative.

    Raises
    ------
    DomainError
        If a sigma is negative or not finite.

    References
    ----------
    .. [1] C. R. Gates, "A Simplified Model of Midcourse Maneuver Execution
           Errors", JPL Technical Report 32-504, 1963.
    .. [2] NASA/TP-2018-219822, "Navigation Filter Best Practices".
           https://ntrs.nasa.gov/citations/20180003657
    """

    fixed_magnitude_sigma_m_s: float
    proportional_magnitude_sigma: float
    fixed_pointing_sigma_m_s: float
    proportional_pointing_sigma_rad: float

    def __post_init__(self) -> None:
        for field_name in (
            "fixed_magnitude_sigma_m_s",
            "proportional_magnitude_sigma",
            "fixed_pointing_sigma_m_s",
            "proportional_pointing_sigma_rad",
        ):
            value = real(getattr(self, field_name), f"Gates {field_name}")
            if value < 0.0:
                raise DomainError("Gates error sigmas must be finite and non-negative.")
            object.__setattr__(self, field_name, value)

    @classmethod
    def from_parameter_prior(
        cls,
        prior: PropulsionParameterPrior,
        *,
        consumed_mass_fraction: float,
        fixed_magnitude_sigma_m_s: float,
        fixed_pointing_sigma_m_s: float,
        proportional_pointing_sigma_rad: float,
    ) -> GatesExecutionErrorModel:
        """Derive the proportional magnitude term from propulsion uncertainty.

        ``ss = sqrt(g^T P_rel g)`` with ``g`` from
        :func:`delta_v_log_sensitivities` and ``P_rel`` the relative
        covariance of ``thrust_scale`` and ``mass_flow_scale`` (a scale not in
        the prior contributes nothing).  Valid for duration-terminated burns;
        an accelerometer (delta-v) cutoff removes most of this term and leaves
        its own resolution error in the fixed magnitude sigma instead [1]_.
        The quadratic form is the first-order propagation of the correlated
        relative scale uncertainty [2]_.

        Parameters
        ----------
        prior : PropulsionParameterPrior
            Prior with nonzero means for the listed scales.
        consumed_mass_fraction : float
            Nominal consumed mass fraction in ``(0, 1)``.
        fixed_magnitude_sigma_m_s, fixed_pointing_sigma_m_s : float
            ``sr`` and ``sa`` [m/s].
        proportional_pointing_sigma_rad : float
            ``sp`` [rad], supplied by attitude control.

        Returns
        -------
        GatesExecutionErrorModel
            Model with the derived ``ss``.

        Raises
        ------
        DomainError
            If the prior type, a mean, the mass fraction or a sigma is invalid.

        References
        ----------
        .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        .. [2] JCGM 100:2008, section 5.2, "Evaluation of Measurement Data -
               Guide to the Expression of Uncertainty in Measurement".
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        """

        if not isinstance(prior, PropulsionParameterPrior):
            raise DomainError("A PropulsionParameterPrior is required.")
        gradient = dict(
            zip(
                ("thrust_scale", "mass_flow_scale"),
                delta_v_log_sensitivities(consumed_mass_fraction),
                strict=True,
            )
        )
        relevant = [name for name in prior.parameter_names if name in gradient]
        positions = {name: prior.parameter_names.index(name) for name in relevant}
        if any(prior.mean[positions[name]] == 0.0 for name in relevant):
            raise DomainError("Thrust and flow scale means must be non-zero.")
        variance = 0.0
        for a in relevant:
            for b in relevant:
                i, j = positions[a], positions[b]
                relative = prior.covariance[i][j] / (prior.mean[i] * prior.mean[j])
                variance += gradient[a] * relative * gradient[b]
        return cls(
            fixed_magnitude_sigma_m_s,
            sqrt(max(variance, 0.0)),
            fixed_pointing_sigma_m_s,
            proportional_pointing_sigma_rad,
        )

    @staticmethod
    def _unit(delta_v: tuple[float, float, float]) -> tuple[tuple[float, float, float], float]:
        """Return the unit direction and magnitude of a nonzero delta-v.

        Raises
        ------
        DomainError
            If ``delta_v`` is not a finite 3-vector or is zero.
        """
        delta_v = vector3(delta_v, "Delta-v")
        magnitude = sqrt(sum(value * value for value in delta_v))
        if magnitude <= 0.0:
            raise DomainError("Delta-v must be non-zero.")
        return tuple(value / magnitude for value in delta_v), magnitude  # type: ignore[return-value]

    def along_and_transverse_variance(self, delta_v_magnitude_m_s: float) -> tuple[float, float]:
        """Return the along-burn and per-transverse-axis variances [m^2/s^2].

        ``(sr**2 + V**2 ss**2, sa**2 + V**2 sp**2)`` for magnitude ``V`` [1]_.

        Raises
        ------
        DomainError
            If the magnitude is negative or not finite.

        References
        ----------
        .. [1] C. R. Gates, "A Simplified Model of Midcourse Maneuver Execution
               Errors", JPL Technical Report 32-504, 1963.
        """

        delta_v_magnitude_m_s = real(delta_v_magnitude_m_s, "Delta-v magnitude")
        if delta_v_magnitude_m_s < 0.0:
            raise DomainError("Delta-v magnitude must be finite and non-negative.")
        v2 = delta_v_magnitude_m_s**2
        return (
            self.fixed_magnitude_sigma_m_s**2 + v2 * self.proportional_magnitude_sigma**2,
            self.fixed_pointing_sigma_m_s**2 + v2 * self.proportional_pointing_sigma_rad**2,
        )

    def covariance(self, delta_v: tuple[float, float, float]) -> Matrix:
        """Return the 3x3 execution-error covariance in the frame of ``delta_v``.

        ``Q = s_along**2 e e^T + s_trans**2 (I - e e^T)`` [1]_ [2]_.

        Parameters
        ----------
        delta_v : tuple of float
            Nonzero commanded velocity change [m/s].

        Returns
        -------
        tuple of tuple of float
            Symmetric covariance [m^2/s^2].

        Raises
        ------
        DomainError
            If ``delta_v`` is invalid or zero.

        References
        ----------
        .. [1] C. R. Gates, "A Simplified Model of Midcourse Maneuver Execution
               Errors", JPL Technical Report 32-504, 1963.
        .. [2] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        unit, magnitude = self._unit(delta_v)
        along, transverse = self.along_and_transverse_variance(magnitude)
        return tuple(
            tuple(
                along * unit[i] * unit[j]
                + transverse * ((1.0 if i == j else 0.0) - unit[i] * unit[j])
                for j in range(3)
            )
            for i in range(3)
        )

    def sample(
        self, delta_v: tuple[float, float, float], generator: random.Random
    ) -> tuple[float, float, float]:
        """Return one executed delta-v drawn from the Gaussian error model.

        Independent normal errors are drawn along ``e`` and along two unit
        axes perpendicular to it, with the variances of
        :meth:`along_and_transverse_variance` [1]_ [2]_.

        Parameters
        ----------
        delta_v : tuple of float
            Nonzero commanded velocity change [m/s].
        generator : random.Random
            Random stream; the draw is reproducible for a given state.

        Returns
        -------
        tuple of float
            Executed velocity change [m/s].

        Raises
        ------
        DomainError
            If the generator type or ``delta_v`` is invalid.

        References
        ----------
        .. [1] C. R. Gates, "A Simplified Model of Midcourse Maneuver Execution
               Errors", JPL Technical Report 32-504, 1963.
        .. [2] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        if not isinstance(generator, random.Random):
            raise DomainError("A random.Random generator is required.")
        unit, magnitude = self._unit(delta_v)
        delta_v = vector3(delta_v, "Delta-v")
        along, transverse = self.along_and_transverse_variance(magnitude)
        helper = (1.0, 0.0, 0.0) if abs(unit[0]) < 0.9 else (0.0, 1.0, 0.0)
        first = _cross(unit, helper)
        norm = sqrt(sum(value * value for value in first))
        first = tuple(value / norm for value in first)
        second = _cross(unit, first)
        a = generator.gauss(0.0, sqrt(along))
        b = generator.gauss(0.0, sqrt(transverse))
        c = generator.gauss(0.0, sqrt(transverse))
        return tuple(
            delta_v[i] + a * unit[i] + b * first[i] + c * second[i] for i in range(3)
        )  # type: ignore[return-value]



@dataclass(frozen=True, slots=True)
class FirstOrderGaussMarkov:
    """Exponentially correlated process noise for one driver deviation.

    The deviation ``x`` (for example ``thrust_scale - 1``) obeys
    ``dx/dt = -x / tau + w`` with white noise of spectral density
    ``q = 2 sigma^2 / tau``, so the steady-state standard deviation is
    ``sigma``.  Over a step ``dt`` the exact discrete form is
    ``x_{k+1} = phi x_k + w_k`` with ``phi = exp(-dt / tau)`` and
    ``Var(w_k) = sigma^2 (1 - phi^2)`` [1]_.

    Parameters
    ----------
    parameter_name : str
        Name of the driver whose deviation is modeled.
    steady_state_sigma : float
        ``sigma``, finite and nonnegative, in the driver's units.
    time_constant_s : float
        Correlation time ``tau`` [s], finite and positive.

    Raises
    ------
    DomainError
        If a value is invalid.

    References
    ----------
    .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
           https://ntrs.nasa.gov/citations/20180003657
    """

    parameter_name: str
    steady_state_sigma: float
    time_constant_s: float

    def __post_init__(self) -> None:
        text(self.parameter_name, "Gauss-Markov parameter name")
        object.__setattr__(
            self, "steady_state_sigma", real(self.steady_state_sigma, "Gauss-Markov sigma")
        )
        object.__setattr__(
            self, "time_constant_s", real(self.time_constant_s, "Gauss-Markov time constant")
        )
        if self.steady_state_sigma < 0.0:
            raise DomainError("Gauss-Markov sigma must be finite and non-negative.")
        if self.time_constant_s <= 0.0:
            raise DomainError("Gauss-Markov time constant must be finite and positive.")

    @property
    def power_spectral_density(self) -> float:
        """Return the white-noise spectral density ``q = 2 sigma^2 / tau`` [unit^2/s]."""

        return 2.0 * self.steady_state_sigma**2 / self.time_constant_s

    def transition(self, dt_s: float) -> float:
        """Return the state transition ``phi = exp(-|dt| / tau)`` over ``dt_s``.

        Raises
        ------
        DomainError
            If ``dt_s`` is not finite.
        """

        dt_s = real(dt_s, "Gauss-Markov step")
        return exp(-abs(dt_s) / self.time_constant_s)

    def discrete_noise_variance(self, dt_s: float) -> float:
        """Return the exact discrete noise variance ``sigma^2 (1 - phi^2)`` over ``dt_s``."""

        phi = self.transition(dt_s)
        return self.steady_state_sigma**2 * (1.0 - phi * phi)

    def propagate(self, mean: float, variance: float, dt_s: float) -> tuple[float, float]:
        """Propagate the deviation's mean and variance over ``dt_s``.

        ``m' = phi m`` and ``P' = phi**2 P + sigma**2 (1 - phi**2)`` [1]_.

        Parameters
        ----------
        mean, variance : float
            Current mean and nonnegative variance.
        dt_s : float
            Step [s]; its sign is ignored.

        Returns
        -------
        tuple of float
            Propagated ``(mean, variance)``.

        Raises
        ------
        DomainError
            If an input is not finite or the variance is negative.

        References
        ----------
        .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        mean = real(mean, "Gauss-Markov mean")
        variance = real(variance, "Gauss-Markov variance")
        if variance < 0.0:
            raise DomainError("Gauss-Markov mean and variance must be finite and valid.")
        phi = self.transition(dt_s)
        return phi * mean, phi * phi * variance + self.discrete_noise_variance(dt_s)

    def sample_path(
        self, times_s: tuple[float, ...], generator: random.Random, *, initial: float | None = None
    ) -> tuple[float, ...]:
        """Draw one exact discrete realization at increasing ``times_s``.

        Without ``initial`` the first value is drawn from the steady state
        ``N(0, sigma**2)``. Each later value applies the exact discrete
        transition and noise variance for its step [1]_.

        Parameters
        ----------
        times_s : tuple of float
            Finite, strictly increasing sample times [s].
        generator : random.Random
            Random stream.
        initial : float or None, optional
            Deviation at the first time.

        Returns
        -------
        tuple of float
            Deviation at each time; empty for no times.

        Raises
        ------
        DomainError
            If the generator type, times or initial value are invalid.

        References
        ----------
        .. [1] NASA/TP-2018-219822, "Navigation Filter Best Practices".
               https://ntrs.nasa.gov/citations/20180003657
        """

        if not isinstance(generator, random.Random):
            raise DomainError("A random.Random generator is required.")
        times_s = real_tuple(times_s, "Gauss-Markov sample times")
        if not times_s:
            return ()
        if any(
            later <= earlier for earlier, later in pairwise(times_s)
        ):
            raise DomainError("Gauss-Markov sample times must be finite and increasing.")
        value = (
            generator.gauss(0.0, self.steady_state_sigma)
            if initial is None
            else real(initial, "Gauss-Markov initial value")
        )
        values = [value]
        for earlier, later in pairwise(times_s):
            dt = later - earlier
            value = self.transition(dt) * value + generator.gauss(
                0.0, sqrt(self.discrete_noise_variance(dt))
            )
            values.append(value)
        return tuple(values)
