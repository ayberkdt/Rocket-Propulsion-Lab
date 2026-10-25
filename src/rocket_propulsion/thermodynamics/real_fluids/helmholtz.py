"""Reduced Helmholtz-energy equations of state for pure fluids.

A fundamental equation of state gives the molar Helmholtz energy ``a`` as a
function of temperature ``T`` and molar density ``rho``. In reduced form

    alpha(tau, delta) = a / (R T) = alpha0(tau, delta) + alphar(tau, delta),

with ``tau = T_r / T`` and ``delta = rho / rho_r`` [1]_ [2]_. The ideal-gas
part ``alpha0`` follows from an ideal-gas heat capacity; the residual part
``alphar`` is an empirical sum of polynomial terms ``n delta^d tau^t``,
exponential terms ``n delta^d tau^t exp(-c delta^l)`` and Gaussian bell-shaped
terms ``n delta^d tau^t exp(-eta (delta - epsilon)^2 - beta (tau - gamma)^2)``
[1]_ (Eq. 32), [2]_ (Eq. 55), [3]_ (Eq. 5.3), [4]_ (Eq. 4).

All thermodynamic properties follow from ``alpha`` and its first and second
derivatives, using the relations printed as Eqs. (56)-(64) and (68)-(71) of
Span et al. [2]_. Saturation states are found from the Maxwell criterion,
equal pressure and equal Gibbs energy of the coexisting phases at one
temperature [2]_ (p. 1365); a vapor-pressure ancillary equation of the
source paper, when one is given, supplies the starting pressure.

The module is a framework only: the coefficient sets of the reference fluids
live in :mod:`rocket_propulsion.thermodynamics.real_fluids.reference_fluids`
and the cubic equations in
:mod:`rocket_propulsion.thermodynamics.real_fluids.cubic`. Every calculator
refuses states outside the temperature and pressure range declared by the
fluid's source, and refuses two-phase ``(T, rho)`` states.

References
----------
.. [1] LEACHMAN-2009-H2. J. W. Leachman, R. T Jacobsen, S. G. Penoncello and
       E. W. Lemmon, "Fundamental Equations of State for Parahydrogen, Normal
       Hydrogen, and Orthohydrogen", J. Phys. Chem. Ref. Data 38(3), 721-748,
       2009, Eqs. (31)-(32), p. 730,
       https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
.. [2] SPAN-2000-N2. R. Span, E. W. Lemmon, R. T Jacobsen, W. Wagner and
       A. Yokozeki, "A Reference Equation of State for the Thermodynamic
       Properties of Nitrogen for Temperatures from 63.151 to 1000 K and
       Pressures to 2200 MPa", J. Phys. Chem. Ref. Data 29(6), 1361-1433,
       2000, Eqs. (55)-(71), pp. 1384-1386,
       https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
.. [3] SETZMANN-WAGNER-1991-CH4. U. Setzmann and W. Wagner, "A New Equation
       of State and Tables of Thermodynamic Properties for Methane Covering
       the Range from the Melting Line to 625 K at Pressures up to 1000 MPa",
       J. Phys. Chem. Ref. Data 20(6), 1061-1155, 1991, Eq. (5.3), p. 1091,
       https://srd.nist.gov/jpcrdreprint/1.555898.pdf
.. [4] ORTIZVEGA-2023-HE. D. O. Ortiz Vega, K. R. Hall, J. C. Holste,
       A. H. Harvey and E. W. Lemmon, "An Equation of State for the
       Thermodynamic Properties of Helium", NIST IR 8474, 2023, Eq. (4), p. 3,
       https://doi.org/10.6028/NIST.IR.8474
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, isfinite, log, log1p
from typing import Literal, Protocol

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

Phase = Literal["liquid", "vapor", "supercritical"]
TermKind = Literal["power", "exponential", "gaussian"]
ResidualDerivatives = tuple[float, float, float, float, float, float]

_TERM_KINDS = ("power", "exponential", "gaussian")
_DENSITY_REL_TOL = 1e-14
_SAT_STEP_TOL = 1e-12
_SAT_RESIDUAL_TOL = 1e-13
_SCAN_LOG_FACTOR = 1.15
_SCAN_LINEAR_START = 0.05
_SCAN_LINEAR_STEP = 0.01
_SCAN_FIRST_DELTA = 1e-8
_MAX_SCAN_POINTS = 20000


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise InputError(f"{label} must be positive.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class ResidualTerm:
    """One term of an empirical residual Helmholtz-energy sum.

    Parameters
    ----------
    kind : {"power", "exponential", "gaussian"}
        ``"power"`` is ``n delta^d tau^t``; ``"exponential"`` multiplies it
        by ``exp(-c delta^l)``; ``"gaussian"`` multiplies it by
        ``exp(-eta (delta - epsilon)^2 - beta (tau - gamma)^2)`` [1]_ [2]_.
    n : float
        Coefficient.
    t : float
        Temperature exponent.
    d : float
        Density exponent, non-negative.
    l : float
        Density exponent inside the exponential; positive for
        ``"exponential"`` terms, ignored otherwise.
    c : float
        Factor of ``delta^l`` inside the exponential (1 in every reference
        equation of this package), positive.
    eta, beta, gamma, epsilon : float
        Gaussian parameters; ``eta`` and ``beta`` are non-negative. The sign
        convention is the one above; a source that prints the exponent as
        ``+phi (delta - D)^2 + beta (tau - gamma)^2`` with negative ``phi`` and
        ``beta`` [1]_ is converted with ``eta = -phi`` and ``beta = -beta``.

    Raises
    ------
    InputError
        If the kind is unknown or a parameter violates its sign constraint.
    DomainError
        If a parameter is not a finite real number.

    References
    ----------
    .. [1] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), Eq. (32), p. 730,
           https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    .. [2] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eq. (55), p. 1384,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    kind: TermKind
    n: float
    t: float
    d: float
    l: float = 0.0  # noqa: E741 - exponent name used by the sources
    c: float = 1.0
    eta: float = 0.0
    beta: float = 0.0
    gamma: float = 0.0
    epsilon: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in _TERM_KINDS:
            raise InputError(f"Unknown residual term kind {self.kind!r}.", field="kind")
        for name in ("n", "t", "d", "l", "c", "eta", "beta", "gamma", "epsilon"):
            object.__setattr__(self, name, real(getattr(self, name), name))
        if self.d < 0.0:
            raise InputError("Density exponent d must be non-negative.", field="d")
        if self.kind == "exponential" and (self.l <= 0.0 or self.c <= 0.0):
            raise InputError("Exponential terms need positive l and c.", field="l")
        if self.kind == "gaussian" and (self.eta < 0.0 or self.beta < 0.0):
            raise InputError("Gaussian terms need non-negative eta and beta.", field="eta")


@dataclass(frozen=True, slots=True)
class ResidualHelmholtz:
    """Residual reduced Helmholtz energy as a sum of :class:`ResidualTerm`.

    Parameters
    ----------
    terms : tuple of ResidualTerm
        Terms in source order.

    Attributes
    ----------
    reduced_density_limit : float
        Upper bound of ``delta`` used by the root searches; infinite for an
        empirical sum.

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eqs. (55) and (81)-(86), pp. 1384-1386,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    terms: tuple[ResidualTerm, ...]
    reduced_density_limit: float = field(default=float("inf"), init=False)
    _packed: tuple[tuple[float, ...], ...] = field(
        default=(), init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        terms = tuple(self.terms)
        if not terms or not all(isinstance(term, ResidualTerm) for term in terms):
            raise InputError("terms must be a non-empty tuple of ResidualTerm.", field="terms")
        object.__setattr__(self, "terms", terms)
        packed = tuple(
            (
                float(_TERM_KINDS.index(term.kind)), term.n, term.t, term.d, term.l, term.c,
                term.eta, term.beta, term.gamma, term.epsilon,
            )
            for term in terms
        )
        object.__setattr__(self, "_packed", packed)

    def derivatives(self, tau: float, delta: float) -> ResidualDerivatives:
        """Return ``alphar`` and its derivatives at ``(tau, delta)``.

        Parameters
        ----------
        tau, delta : float
            Inverse reduced temperature and reduced density, both positive.

        Returns
        -------
        tuple of float
            ``(alphar, alphar_delta, alphar_delta_delta, alphar_tau,
            alphar_tau_tau, alphar_delta_tau)``.

        Notes
        -----
        Every term is a product ``n f(delta) g(tau)``; with
        ``A = d ln f / d delta`` and ``B = d2 ln f / d delta2`` the
        derivatives are ``f' = f A`` and ``f'' = f (A^2 + B)``, and likewise
        for ``g``. This reproduces the closed forms of Span et al. [1]_
        (Eqs. 81-86) for all three term kinds.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eqs. (81)-(86), p. 1386,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        inv_d = 1.0 / delta
        inv_t = 1.0 / tau
        log_d = log(delta)
        log_t = log(tau)
        ar = ar_d = ar_dd = ar_t = ar_tt = ar_dt = 0.0
        for kind, n, t, d, el, c, eta, beta, gam, eps in self._packed:
            f = exp(d * log_d)
            fa = d * inv_d
            fb = -d * inv_d * inv_d
            g = exp(t * log_t)
            ga = t * inv_t
            gb = -t * inv_t * inv_t
            if kind == 1.0:
                e = c * exp(el * log_d)
                f *= exp(-e)
                fa -= el * e * inv_d
                fb -= el * (el - 1.0) * e * inv_d * inv_d
            elif kind == 2.0:
                x = delta - eps
                y = tau - gam
                f *= exp(-eta * x * x)
                fa -= 2.0 * eta * x
                fb -= 2.0 * eta
                g *= exp(-beta * y * y)
                ga -= 2.0 * beta * y
                gb -= 2.0 * beta
            f1 = f * fa
            g1 = g * ga
            ar += n * f * g
            ar_d += n * f1 * g
            ar_dd += n * f * (fa * fa + fb) * g
            ar_t += n * f * g1
            ar_tt += n * f * g * (ga * ga + gb)
            ar_dt += n * f1 * g1
        return ar, ar_d, ar_dd, ar_t, ar_tt, ar_dt


class ResidualModel(Protocol):
    """Interface of a residual Helmholtz model (empirical sum or cubic).

    Attributes
    ----------
    reduced_density_limit : float
        Largest admissible ``delta`` (covolume limit for a cubic equation).
    """

    reduced_density_limit: float

    def derivatives(self, tau: float, delta: float) -> ResidualDerivatives:
        """Return ``alphar`` and its five first and second derivatives."""


@dataclass(frozen=True, slots=True)
class IdealHelmholtz:
    """Ideal-gas reduced Helmholtz energy ``alpha0(tau, delta)``.

    ``alpha0 = ln(delta) + a1 + a2 tau + log_tau ln(tau) + sum n tau^t
    + sum n ln(1 - exp(-theta tau)) + sum n ln(1 + g exp(-theta tau))``.

    Parameters
    ----------
    a1, a2 : float
        Integration constants fixing the enthalpy and entropy reference.
    log_tau : float
        Coefficient of ``ln(tau)``; equals ``c_p0/R - 1`` for a constant
        ideal-gas heat capacity.
    power : tuple of (n, t)
        Terms ``n tau^t``. A heat-capacity term ``c_p0/R = c T^k``
        (``k`` not 0 or -1) contributes ``-c T_r^k / (k (k + 1)) tau^-k``.
    planck_einstein : tuple of (n, theta)
        Terms ``n ln(1 - exp(-theta tau))``, the Helmholtz form of the
        Planck-Einstein heat capacity ``n u^2 e^u / (e^u - 1)^2`` with
        ``u = theta tau`` [1]_ (Eqs. 25 and 31).
    two_level : tuple of (n, g, theta)
        Terms ``n ln(1 + g exp(-theta tau))``. With ``n = -N`` this is the
        Helmholtz form of the electronic heat-capacity term
        ``N g u^2 e^-u / (1 + g e^-u)^2`` of Stewart et al. [2]_ (Eq. 11).

    References
    ----------
    .. [1] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), Eqs. (25)-(31), pp. 728-730,
           https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    .. [2] STEWART-1991-O2. R. B. Stewart, R. T Jacobsen and W. Wagner,
           "Thermodynamic Properties of Oxygen from the Triple Point to 300 K
           with Pressures to 80 MPa", J. Phys. Chem. Ref. Data 20(5),
           917-1021, 1991, Eqs. (11)-(14), pp. 926-927,
           https://srd.nist.gov/jpcrdreprint/1.555897.pdf
    """

    a1: float
    a2: float
    log_tau: float
    power: tuple[tuple[float, float], ...] = ()
    planck_einstein: tuple[tuple[float, float], ...] = ()
    two_level: tuple[tuple[float, float, float], ...] = ()

    def __post_init__(self) -> None:
        for name in ("a1", "a2", "log_tau"):
            object.__setattr__(self, name, real(getattr(self, name), name))
        object.__setattr__(self, "power", tuple(
            (real(n, "power n"), real(t, "power t")) for n, t in self.power))
        object.__setattr__(self, "planck_einstein", tuple(
            (real(n, "planck_einstein n"), _positive(th, "planck_einstein theta"))
            for n, th in self.planck_einstein))
        object.__setattr__(self, "two_level", tuple(
            (real(n, "two_level n"), _positive(g, "two_level g"),
             _positive(th, "two_level theta")) for n, g, th in self.two_level))

    def evaluate(self, tau: float, delta: float) -> tuple[float, float, float]:
        """Return ``(alpha0, alpha0_tau, alpha0_tau_tau)``.

        Parameters
        ----------
        tau, delta : float
            Inverse reduced temperature and reduced density, both positive.

        Returns
        -------
        tuple of float
            The function and its first and second ``tau`` derivatives; the
            ``delta`` derivatives are ``1/delta`` and ``-1/delta^2``.

        References
        ----------
        .. [1] SETZMANN-WAGNER-1991-CH4. Setzmann and Wagner, J. Phys. Chem.
               Ref. Data 20, 1061 (1991), Table 34, p. 1092,
               https://srd.nist.gov/jpcrdreprint/1.555898.pdf
        """

        a = log(delta) + self.a1 + self.a2 * tau + self.log_tau * log(tau)
        a_t = self.a2 + self.log_tau / tau
        a_tt = -self.log_tau / (tau * tau)
        for n, t in self.power:
            a += n * tau ** t
            a_t += n * t * tau ** (t - 1.0)
            a_tt += n * t * (t - 1.0) * tau ** (t - 2.0)
        for n, theta in self.planck_einstein:
            em = exp(-theta * tau)
            one = 1.0 - em
            a += n * log1p(-em)
            a_t += n * theta * em / one
            a_tt -= n * theta * theta * em / (one * one)
        for n, g, theta in self.two_level:
            em = g * exp(-theta * tau)
            one = 1.0 + em
            a += n * log1p(em)
            a_t -= n * theta * em / one
            a_tt += n * theta * theta * em / (one * one)
        return a, a_t, a_tt


@dataclass(frozen=True, slots=True)
class VaporPressureAncillary:
    """Ancillary vapor-pressure equation ``ln(p/pc) = (Tc/T) sum n theta^k``.

    ``theta = 1 - T/Tc``. It only seeds the Maxwell solve; reported
    saturation states always come from the equation of state, as the sources
    recommend [1]_ (p. 924) [2]_ (Eq. 33, p. 733).

    Parameters
    ----------
    t_c_k, p_c_pa : float
        Critical temperature and pressure used by the ancillary equation.
    terms : tuple of (n, k)
        Coefficients and exponents.
    source : str
        Citation key and locator of the coefficients.

    References
    ----------
    .. [1] STEWART-1991-O2. Stewart et al., J. Phys. Chem. Ref. Data 20, 917
           (1991), Eq. (7), pp. 924-925,
           https://srd.nist.gov/jpcrdreprint/1.555897.pdf
    .. [2] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), Eq. (33) and Table 8, pp. 732-733,
           https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    """

    t_c_k: float
    p_c_pa: float
    terms: tuple[tuple[float, float], ...]
    source: str

    def pressure(self, temperature_k: float) -> float:
        """Return the ancillary vapor pressure in Pa at ``temperature_k``.

        Parameters
        ----------
        temperature_k : float
            Temperature below ``t_c_k``.

        Returns
        -------
        float
            Approximate saturation pressure in Pa.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eq. (4), p. 1368,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        theta = 1.0 - temperature_k / self.t_c_k
        total = sum(n * theta ** k for n, k in self.terms)
        return self.p_c_pa * exp(self.t_c_k / temperature_k * total)


@dataclass(frozen=True, slots=True)
class LiquidDensityAncillary:
    """Ancillary saturated-liquid density equation of a source paper.

    ``form="linear"``: ``rho'/rho_c = 1 + sum n theta^k``;
    ``form="log"``: ``ln(rho'/rho_c) = sum n theta^k``; ``theta = 1 - T/Tc``.
    It only seeds the liquid root of the Maxwell solve [1]_ (p. 924).

    Parameters
    ----------
    t_c_k, rho_c_mol_m3 : float
        Critical temperature and molar density used by the equation.
    terms : tuple of (n, k)
        Coefficients and exponents.
    form : {"linear", "log"}
        Functional form, see above.
    source : str
        Citation key and locator.

    References
    ----------
    .. [1] STEWART-1991-O2. Stewart et al., J. Phys. Chem. Ref. Data 20, 917
           (1991), Eq. (10), p. 926, and p. 924,
           https://srd.nist.gov/jpcrdreprint/1.555897.pdf
    .. [2] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eq. (5), p. 1368,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    t_c_k: float
    rho_c_mol_m3: float
    terms: tuple[tuple[float, float], ...]
    form: Literal["linear", "log"]
    source: str

    def density(self, temperature_k: float) -> float:
        """Return the ancillary saturated-liquid molar density in mol/m3.

        Parameters
        ----------
        temperature_k : float
            Temperature below ``t_c_k``.

        Returns
        -------
        float
            Approximate saturated-liquid molar density.

        References
        ----------
        .. [1] SETZMANN-WAGNER-1991-CH4. Setzmann and Wagner, J. Phys. Chem.
               Ref. Data 20, 1061 (1991), Eq. (3.4), p. 1082,
               https://srd.nist.gov/jpcrdreprint/1.555898.pdf
        """

        theta = 1.0 - temperature_k / self.t_c_k
        total = sum(n * theta ** k for n, k in self.terms)
        ratio = 1.0 + total if self.form == "linear" else exp(total)
        return self.rho_c_mol_m3 * ratio


@dataclass(frozen=True, slots=True)
class FluidModel:
    """Complete pure-fluid equation of state with its declared validity.

    Parameters
    ----------
    name : str
        Display name.
    key : str
        Registry key.
    model_kind : {"reference", "cubic"}
        ``"cubic"`` models are labelled lower accuracy in every result.
    molar_mass_kg_mol, gas_constant_j_mol_k : float
        Constants used by the source (the source's own gas constant).
    t_reduce_k, rho_reduce_mol_m3 : float
        Reducing temperature and molar density of ``tau`` and ``delta``.
    t_crit_k, p_crit_pa : float
        Critical temperature and pressure (phase labels, saturation limit).
    t_min_k, t_max_k, p_max_pa : float
        Declared validity range; states outside it are refused.
    ideal : IdealHelmholtz
        Ideal-gas part.
    residual : ResidualModel
        Residual part.
    ancillary : VaporPressureAncillary or None
        Optional saturation starting pressure.
    liquid_ancillary : LiquidDensityAncillary or None
        Optional saturated-liquid starting density; without it the liquid
        root is taken from a scan of the isotherm.
    references : tuple of str
        Citation keys with locators.
    notes : tuple of str
        Provenance and limitation notes.

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), physical constants table, p. 1364,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    name: str
    key: str
    model_kind: Literal["reference", "cubic"]
    molar_mass_kg_mol: float
    gas_constant_j_mol_k: float
    t_reduce_k: float
    rho_reduce_mol_m3: float
    t_crit_k: float
    p_crit_pa: float
    t_min_k: float
    t_max_k: float
    p_max_pa: float
    ideal: IdealHelmholtz
    residual: ResidualModel
    ancillary: VaporPressureAncillary | None = None
    liquid_ancillary: LiquidDensityAncillary | None = None
    references: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("molar_mass_kg_mol", "gas_constant_j_mol_k", "t_reduce_k",
                     "rho_reduce_mol_m3", "t_crit_k", "p_crit_pa", "t_min_k", "t_max_k",
                     "p_max_pa"):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if not self.t_min_k < self.t_crit_k < self.t_max_k:
            raise InputError("Need t_min_k < t_crit_k < t_max_k.", field="t_crit_k")
        if self.model_kind not in ("reference", "cubic"):
            raise InputError("model_kind must be 'reference' or 'cubic'.", field="model_kind")


@dataclass(frozen=True, slots=True)
class FluidState:
    """Single-phase thermodynamic state from a Helmholtz-energy model.

    Attributes
    ----------
    fluid, model_kind : str
        Fluid name and ``"reference"`` or ``"cubic"``.
    phase : {"liquid", "vapor", "supercritical"}
        Label in the NIST WebBook convention: ``"liquid"`` below the
        critical temperature above the vapor pressure, ``"vapor"`` below the
        vapor pressure or above ``T_c`` below ``p_c``, ``"supercritical"``
        above both critical values.
    temperature_k, pressure_pa : float
        State temperature and pressure.
    density_kg_m3, molar_density_mol_m3 : float
        Mass and molar density.
    compressibility_factor : float
        ``Z = p / (rho R T)``.
    internal_energy_j_kg, enthalpy_j_kg, entropy_j_kg_k, gibbs_energy_j_kg : float
        Caloric properties on the model's reference state.
    cv_j_kg_k, cp_j_kg_k : float
        Isochoric and isobaric heat capacities.
    speed_of_sound_m_s : float
        Thermodynamic speed of sound.
    joule_thomson_k_pa : float
        Joule-Thomson coefficient ``(dT/dp)_h``.
    dp_drho_t_pa_m3_kg, dp_dt_rho_pa_k : float
        ``(dp/drho)_T`` (mass density) and ``(dp/dT)_rho``.
    warnings : tuple of WarningMessage
        Accuracy and limitation notes.
    """

    fluid: str
    model_kind: str
    phase: Phase
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    molar_density_mol_m3: float
    compressibility_factor: float
    internal_energy_j_kg: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float
    gibbs_energy_j_kg: float
    cv_j_kg_k: float
    cp_j_kg_k: float
    speed_of_sound_m_s: float
    joule_thomson_k_pa: float
    dp_drho_t_pa_m3_kg: float
    dp_dt_rho_pa_k: float
    warnings: tuple[WarningMessage, ...] = ()


@dataclass(frozen=True, slots=True)
class SaturationState:
    """Vapor-liquid equilibrium from the Maxwell criterion.

    Attributes
    ----------
    fluid : str
        Fluid name.
    temperature_k, pressure_pa : float
        Saturation temperature and pressure (pressure from the vapor side).
    liquid, vapor : FluidState
        Coexisting phases.
    enthalpy_of_vaporization_j_kg : float
        ``h_vapor - h_liquid``.
    ancillary_pressure_pa : float or None
        Starting pressure from the source's ancillary equation, if used.
    pressure_mismatch, gibbs_mismatch : float
        Final relative pressure difference and Gibbs-energy difference
        divided by ``R T``.
    iterations : int
        Newton iterations of the Maxwell solve.
    warnings : tuple of WarningMessage
        Accuracy notes.
    """

    fluid: str
    temperature_k: float
    pressure_pa: float
    liquid: FluidState
    vapor: FluidState
    enthalpy_of_vaporization_j_kg: float
    ancillary_pressure_pa: float | None
    pressure_mismatch: float
    gibbs_mismatch: float
    iterations: int
    warnings: tuple[WarningMessage, ...] = ()


@dataclass(frozen=True, slots=True)
class DepartureFunctions:
    """Real-fluid minus ideal-gas properties at the same ``T`` and ``p``.

    Attributes
    ----------
    temperature_k, pressure_pa : float
        State.
    compressibility_factor : float
        ``Z``.
    enthalpy_j_mol : float
        ``h - h_ig(T)``.
    entropy_j_mol_k : float
        ``s - s_ig(T, p)``.
    gibbs_energy_j_mol : float
        ``g - g_ig(T, p) = R T ln(phi)``.
    ln_fugacity_coefficient : float
        ``ln(phi) = alphar + delta alphar_delta - ln Z`` [1]_ (Eq. 65).
    cv_j_mol_k : float
        ``cv - cv_ig``.

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eq. (65), p. 1385,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    temperature_k: float
    pressure_pa: float
    compressibility_factor: float
    enthalpy_j_mol: float
    entropy_j_mol_k: float
    gibbs_energy_j_mol: float
    ln_fugacity_coefficient: float
    cv_j_mol_k: float


class RealFluid:
    """Property calculator over one :class:`FluidModel`.

    Parameters
    ----------
    model : FluidModel
        Equation of state and validity range.

    Notes
    -----
    Saturation solutions are cached per temperature on the instance. The
    calculator is otherwise stateless.

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eqs. (56)-(64), p. 1385, and Maxwell criterion, p. 1365,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    __slots__ = ("_sat_cache", "model")

    def __init__(self, model: FluidModel) -> None:
        if not isinstance(model, FluidModel):
            raise InputError("model must be a FluidModel.", field="model")
        self.model = model
        self._sat_cache: dict[float, tuple[float, float, float | None, int, float, float]] = {}

    # ------------------------------------------------------------------ validation
    def _temperature(self, value: object, *, saturation: bool = False) -> float:
        m = self.model
        temperature = real(value, "temperature_k")
        if not m.t_min_k <= temperature <= m.t_max_k:
            raise DomainError(
                f"{m.name}: T = {temperature} K is outside the validity range "
                f"[{m.t_min_k}, {m.t_max_k}] K.", field="temperature_k",
                details={"t_min_k": m.t_min_k, "t_max_k": m.t_max_k})
        if saturation and temperature >= m.t_crit_k:
            raise DomainError(
                f"{m.name}: saturation needs T < T_c = {m.t_crit_k} K.", field="temperature_k")
        return temperature

    def _pressure_value(self, value: object) -> float:
        m = self.model
        pressure = real(value, "pressure_pa")
        if not 0.0 < pressure <= m.p_max_pa:
            raise DomainError(
                f"{m.name}: p = {pressure} Pa is outside the validity range "
                f"(0, {m.p_max_pa}] Pa.", field="pressure_pa", details={"p_max_pa": m.p_max_pa})
        return pressure

    # ------------------------------------------------------------------ kernels
    def _jkd(self, tau: float, delta: float) -> tuple[float, float, float]:
        ar, ar_d, ar_dd, *_ = self.model.residual.derivatives(tau, delta)
        j = delta * (1.0 + delta * ar_d)
        k = log(delta) + ar + delta * ar_d
        dd = 1.0 + 2.0 * delta * ar_d + delta * delta * ar_dd
        return j, k, dd

    def _p_of_delta(self, temperature: float, delta: float) -> tuple[float, float]:
        m = self.model
        tau = m.t_reduce_k / temperature
        _ar, ar_d, ar_dd, *_ = m.residual.derivatives(tau, delta)
        rt_rho = m.gas_constant_j_mol_k * temperature * m.rho_reduce_mol_m3
        p = rt_rho * delta * (1.0 + delta * ar_d)
        dp = rt_rho * (1.0 + 2.0 * delta * ar_d + delta * delta * ar_dd)
        return p, dp

    def _solve_branch(self, temperature: float, pressure: float, lo: float, hi: float,
                      start: float | None = None) -> float:
        """Safeguarded Newton for ``p(delta) = pressure`` on an increasing bracket."""

        x = start if start is not None and lo < start < hi else (lo * hi) ** 0.5
        for _ in range(300):
            p, dp = self._p_of_delta(temperature, x)
            f = p - pressure
            if f > 0.0:
                hi = x
            else:
                lo = x
            if f == 0.0 or abs(f) <= 1e-15 * pressure:
                return x
            step_ok = dp > 0.0 and isfinite(dp)
            xn = x - f / dp if step_ok else -1.0
            if not lo < xn < hi:
                xn = (lo * hi) ** 0.5 if hi > 4.0 * lo else 0.5 * (lo + hi)
            if abs(xn - x) <= _DENSITY_REL_TOL * x or hi - lo <= _DENSITY_REL_TOL * hi:
                return xn
            x = xn
        raise ConvergenceError(
            f"{self.model.name}: density solve at T = {temperature} K, p = {pressure} Pa "
            "did not converge.")

    def _upper_bracket(self, temperature: float, pressure: float, lo: float) -> float:
        limit = self.model.residual.reduced_density_limit
        hi = max(lo, 1e-12)
        if hi >= limit:
            hi = 0.5 * limit
        for _ in range(2000):
            nxt = hi * 1.05
            if nxt >= limit:
                nxt = 0.5 * (hi + limit)
            hi = nxt
            p, dp = self._p_of_delta(temperature, hi)
            if p > pressure and dp > 0.0:
                return hi
            if hi >= 1e3 or (isfinite(limit) and limit - hi <= 1e-12 * limit):
                break
        raise ConvergenceError(
            f"{self.model.name}: no density reaches p = {pressure} Pa at T = {temperature} K.")

    def _local_root(self, temperature: float, pressure: float, start: float) -> float | None:
        """Bracket ``p(delta) = pressure`` on a stable branch next to ``start``."""

        p0, dp0 = self._p_of_delta(temperature, start)
        if dp0 <= 0.0:
            return None
        factor = 1.002 if p0 < pressure else 1.0 / 1.002
        x = start
        for _ in range(500):
            nxt = x * factor
            p, dp = self._p_of_delta(temperature, nxt)
            if dp <= 0.0:
                return None
            if (p - pressure) * (p0 - pressure) <= 0.0:
                lo, hi = (x, nxt) if nxt > x else (nxt, x)
                return self._solve_branch(temperature, pressure, lo, hi)
            x = nxt
        return None

    def _scan(self, temperature: float, p_stop: float) -> list[tuple[float, float]]:
        limit = self.model.residual.reduced_density_limit
        points: list[tuple[float, float]] = []
        delta = _SCAN_FIRST_DELTA
        passed_max = False
        while len(points) < _MAX_SCAN_POINTS:
            p, dp = self._p_of_delta(temperature, delta)
            points.append((delta, p))
            if dp <= 0.0:
                passed_max = True
            if passed_max and dp > 0.0 and p > p_stop:
                break
            if delta < _SCAN_LINEAR_START:
                delta *= _SCAN_LOG_FACTOR
            else:
                delta += _SCAN_LINEAR_STEP
            if delta >= limit:
                delta = 0.5 * (points[-1][0] + limit)
                if limit - delta <= 1e-9 * limit:
                    break
        return points

    def _saturation_deltas(self, temperature: float) -> tuple[float, float, float | None, int,
                                                               float, float]:
        cached = self._sat_cache.get(temperature)
        if cached is not None:
            return cached
        m = self.model
        tau = m.t_reduce_k / temperature
        anc_p = None
        if m.ancillary is not None and temperature < m.ancillary.t_c_k:
            anc_p = m.ancillary.pressure(temperature)
        scan = self._scan(temperature, min(m.p_max_pa, max(100.0 * m.p_crit_pa,
                                                         2.0 * (anc_p or 0.0))))
        pressures = [p for _, p in scan]
        imax = next((i for i in range(1, len(scan) - 1)
                     if pressures[i] >= pressures[i - 1] and pressures[i] > pressures[i + 1]), None)
        if imax is None:
            raise ConvergenceError(
                f"{m.name}: no van der Waals loop found at T = {temperature} K; the state is "
                "too close to the critical point for the Maxwell solve.", field="temperature_k")
        minima = [i for i in range(imax + 1, len(scan) - 1)
                  if pressures[i] <= pressures[i - 1] and pressures[i] < pressures[i + 1]]
        imin = minima[-1] if minima else min(range(imax, len(scan)), key=pressures.__getitem__)
        p_spin_v, p_spin_l = pressures[imax], pressures[imin]
        guess = anc_p
        if guess is None or not max(p_spin_l, 0.0) < guess < p_spin_v:
            guess = 0.5 * (max(p_spin_l, 0.0) + p_spin_v)
        iv = next(i for i in range(imax) if pressures[i] < guess <= pressures[i + 1])
        delta_v = self._solve_branch(temperature, guess, scan[iv][0], scan[iv + 1][0])
        delta_l = None
        if m.liquid_ancillary is not None and temperature < m.liquid_ancillary.t_c_k:
            start = m.liquid_ancillary.density(temperature) / m.rho_reduce_mol_m3
            delta_l = self._local_root(temperature, guess, start)
        if delta_l is None:
            il = next(i for i in range(len(scan) - 2, imin - 1, -1)
                      if pressures[i] < guess <= pressures[i + 1])
            delta_l = self._solve_branch(temperature, guess, scan[il][0], scan[il + 1][0])
        result = self._maxwell(temperature, tau, delta_l, delta_v)
        if not (result[1] <= scan[imax + 1][0] and result[0] >= scan[imin - 1][0]):
            raise ConvergenceError(
                f"{m.name}: Maxwell solve at T = {temperature} K converged to a metastable or "
                "unstable root pair.", field="temperature_k")
        out = (result[0], result[1], anc_p, result[2], result[3], result[4])
        self._sat_cache[temperature] = out
        return out

    def _maxwell(self, temperature: float, tau: float, delta_l: float,
                 delta_v: float) -> tuple[float, float, int, float, float]:
        ln_l, ln_v = log(delta_l), log(delta_v)
        for iteration in range(1, 101):
            jl, kl, dl = self._jkd(tau, delta_l)
            jv, kv, dv = self._jkd(tau, delta_v)
            if dl <= 0.0 or dv <= 0.0:
                raise ConvergenceError(
                    f"{self.model.name}: Maxwell iteration left the stable branches at "
                    f"T = {temperature} K.", field="temperature_k")
            f1 = jl - jv
            f2 = kl - kv
            if iteration > 1 and abs(f1) <= _SAT_RESIDUAL_TOL * jv and abs(f2) <= _SAT_RESIDUAL_TOL:
                break
            a11, a12, a21, a22 = delta_l * dl, -delta_v * dv, dl, -dv
            det = a11 * a22 - a12 * a21
            x1 = (-f1 * a22 + a12 * f2) / det
            x2 = (-a11 * f2 + f1 * a21) / det
            scale = max(abs(x1), abs(x2))
            if scale > 0.25:
                x1 *= 0.25 / scale
                x2 *= 0.25 / scale
            ln_l += x1
            ln_v += x2
            delta_l, delta_v = exp(ln_l), exp(ln_v)
            if scale <= _SAT_STEP_TOL:
                break
        else:
            raise ConvergenceError(
                f"{self.model.name}: Maxwell solve did not converge at T = {temperature} K.",
                field="temperature_k")
        if delta_l <= delta_v * (1.0 + 1e-7):
            raise ConvergenceError(
                f"{self.model.name}: Maxwell solve collapsed to one phase at T = "
                f"{temperature} K (too close to the critical point).", field="temperature_k")
        jl, kl, _ = self._jkd(tau, delta_l)
        jv, kv, _ = self._jkd(tau, delta_v)
        return delta_l, delta_v, iteration, (jl - jv) / jv, kl - kv

    def _phase(self, temperature: float, pressure: float, p_sat: float | None) -> Phase:
        m = self.model
        if temperature < m.t_crit_k and p_sat is not None:
            return "liquid" if pressure > p_sat else "vapor"
        return "supercritical" if pressure >= m.p_crit_pa else "vapor"

    def _warnings(self, phase: Phase) -> tuple[WarningMessage, ...]:
        found: list[WarningMessage] = []
        if self.model.model_kind == "cubic":
            found.append(WarningMessage(
                "cubic_eos_lower_accuracy",
                f"{self.model.name} is a cubic equation of state: a fast fallback with "
                "percent-level (or worse, near saturation and in the liquid) errors; use a "
                "reference Helmholtz equation where one exists."))
        if phase == "liquid":
            found.append(WarningMessage(
                "melting_line_not_checked",
                "The melting line is not modelled; compressed-liquid states near the triple "
                "point are not checked for solidification.", severity="info"))
        return tuple(found)

    def _build(self, temperature: float, delta: float, phase: Phase) -> FluidState:
        m = self.model
        r = m.gas_constant_j_mol_k
        mm = m.molar_mass_kg_mol
        tau = m.t_reduce_k / temperature
        rho = delta * m.rho_reduce_mol_m3
        a0, a0_t, a0_tt = m.ideal.evaluate(tau, delta)
        ar, ar_d, ar_dd, ar_t, ar_tt, ar_dt = m.residual.derivatives(tau, delta)
        z = 1.0 + delta * ar_d
        rt = r * temperature
        u_rt = tau * (a0_t + ar_t)
        cv_r = -tau * tau * (a0_tt + ar_tt)
        num = 1.0 + delta * ar_d - delta * tau * ar_dt
        den = 1.0 + 2.0 * delta * ar_d + delta * delta * ar_dd
        if not (den > 0.0 and cv_r > 0.0):
            raise DomainError(
                f"{m.name}: state T = {temperature} K, rho = {rho} mol/m3 is mechanically or "
                "thermally unstable.", field="density_kg_m3")
        cp_r = cv_r + num * num / den
        w2 = rt / mm * (den + num * num / cv_r)
        dp_drho = rt * den
        dp_dt = rho * r * num
        beta = dp_dt / (rho * dp_drho)
        cp_molar = cp_r * r
        jt = (temperature * beta - 1.0) / (rho * cp_molar)
        return FluidState(
            fluid=m.name, model_kind=m.model_kind, phase=phase,
            temperature_k=temperature, pressure_pa=rho * rt * z,
            density_kg_m3=rho * mm, molar_density_mol_m3=rho, compressibility_factor=z,
            internal_energy_j_kg=u_rt * rt / mm, enthalpy_j_kg=(u_rt + z) * rt / mm,
            entropy_j_kg_k=(tau * (a0_t + ar_t) - a0 - ar) * r / mm,
            gibbs_energy_j_kg=(1.0 + a0 + ar + delta * ar_d) * rt / mm,
            cv_j_kg_k=cv_r * r / mm, cp_j_kg_k=cp_molar / mm, speed_of_sound_m_s=w2 ** 0.5,
            joule_thomson_k_pa=jt, dp_drho_t_pa_m3_kg=dp_drho / mm, dp_dt_rho_pa_k=dp_dt,
            warnings=self._warnings(phase))

    # ------------------------------------------------------------------ public API
    def saturation_t(self, temperature_k: float) -> SaturationState:
        """Return the saturated liquid and vapor at ``temperature_k``.

        Parameters
        ----------
        temperature_k : float
            Temperature in ``[t_min_k, t_crit_k)``.

        Returns
        -------
        SaturationState
            Coexisting states satisfying equal pressure and Gibbs energy.

        Raises
        ------
        DomainError
            If the temperature is outside the range or not below ``T_c``.
        ConvergenceError
            If the Maxwell solve fails (typically within about 1e-5 ``T_c``
            of the critical point).

        Notes
        -----
        The unknowns are ``ln(delta')`` and ``ln(delta'')``; the residuals
        are the differences of ``J = delta (1 + delta alphar_delta)`` (reduced
        pressure, Eq. 56 of [1]_) and of ``K = ln(delta) + alphar +
        delta alphar_delta`` (Gibbs energy over ``RT`` minus terms common to
        both phases, Eq. 61 of [1]_). Starting densities are the stable roots
        of ``p(delta) = p_anc`` where ``p_anc`` is the source's ancillary
        vapor pressure, or the mid pressure of the spinodal loop when no
        ancillary is available.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eqs. (56) and (61), p. 1385, Maxwell criterion,
               p. 1365, https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        temperature = self._temperature(temperature_k, saturation=True)
        delta_l, delta_v, anc_p, iterations, p_mis, g_mis = self._saturation_deltas(temperature)
        liquid = self._build(temperature, delta_l, "liquid")
        vapor = self._build(temperature, delta_v, "vapor")
        return SaturationState(
            fluid=self.model.name, temperature_k=temperature, pressure_pa=vapor.pressure_pa,
            liquid=liquid, vapor=vapor,
            enthalpy_of_vaporization_j_kg=vapor.enthalpy_j_kg - liquid.enthalpy_j_kg,
            ancillary_pressure_pa=anc_p, pressure_mismatch=p_mis, gibbs_mismatch=g_mis,
            iterations=iterations, warnings=self._warnings("vapor"))

    def saturation_p(self, pressure_pa: float) -> SaturationState:
        """Return the saturation state whose vapor pressure is ``pressure_pa``.

        Parameters
        ----------
        pressure_pa : float
            Pressure between the vapor pressure at ``t_min_k`` and the vapor
            pressure at ``(1 - 1e-4) T_c``.

        Returns
        -------
        SaturationState
            Saturation state at the matching temperature.

        Raises
        ------
        DomainError
            If no saturation temperature in the solvable range matches.
        ConvergenceError
            If the temperature iteration fails.

        Notes
        -----
        Regula falsi (Illinois) iteration on ``ln p_sat`` against ``1/T``.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Maxwell criterion, p. 1365,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        pressure = self._pressure_value(pressure_pa)
        m = self.model
        t_lo, t_hi = m.t_min_k, m.t_crit_k * (1.0 - 1e-4)
        target = log(pressure)

        def f(temperature: float) -> float:
            return log(self.saturation_t(temperature).pressure_pa) - target

        f_lo, f_hi = f(t_lo), f(t_hi)
        if not f_lo <= 0.0 <= f_hi:
            raise DomainError(
                f"{m.name}: p = {pressure} Pa has no solvable saturation temperature in "
                f"[{t_lo}, {t_hi}] K.", field="pressure_pa")
        x_lo, x_hi = 1.0 / t_hi, 1.0 / t_lo  # increasing x = 1/T decreases ln p
        g_lo, g_hi = f_hi, f_lo
        side = 0
        x = x_lo
        for _ in range(200):
            x = (x_lo * g_hi - x_hi * g_lo) / (g_hi - g_lo)
            gx = f(1.0 / x)
            if abs(gx) <= 1e-13 or (x_hi - x_lo) <= 1e-15 * x:
                break
            if gx > 0.0:
                x_lo, g_lo = x, gx
                if side == 1:
                    g_hi *= 0.5
                side = 1
            else:
                x_hi, g_hi = x, gx
                if side == -1:
                    g_lo *= 0.5
                side = -1
        else:
            raise ConvergenceError(f"{m.name}: saturation temperature iteration failed.")
        return self.saturation_t(1.0 / x)

    def state_tp(self, temperature_k: float, pressure_pa: float) -> FluidState:
        """Return the stable single-phase state at ``(T, p)``.

        Parameters
        ----------
        temperature_k : float
            Temperature in the model's validity range.
        pressure_pa : float
            Pressure in ``(0, p_max]``.

        Returns
        -------
        FluidState
            Liquid, vapor or supercritical state.

        Raises
        ------
        DomainError
            Outside the validity range, or exactly on the saturation curve.
        ConvergenceError
            If the density solve fails, or below ``T_c`` within the range
            where the saturation solve cannot classify the phase.

        Notes
        -----
        Below ``T_c`` the vapor pressure decides the branch and the saturated
        densities bound the root, so metastable roots are never returned.
        Above ``T_c`` the isotherm is searched upward from the ideal-gas
        density.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eq. (56), p. 1385,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        temperature = self._temperature(temperature_k)
        pressure = self._pressure_value(pressure_pa)
        m = self.model
        rho_ig = pressure / (m.gas_constant_j_mol_k * temperature) / m.rho_reduce_mol_m3
        if temperature < m.t_crit_k:
            delta_l, delta_v, *_ = self._saturation_deltas(temperature)
            p_sat, _ = self._p_of_delta(temperature, delta_v)
            if abs(pressure - p_sat) <= 1e-12 * p_sat:
                raise DomainError(
                    f"{m.name}: (T, p) lies on the saturation curve; use saturation_t.",
                    field="pressure_pa")
            if pressure > p_sat:
                hi = self._upper_bracket(temperature, pressure, delta_l)
                delta = self._solve_branch(temperature, pressure, delta_l, hi)
                return self._build(temperature, delta, "liquid")
            lo = min(rho_ig * 1e-3, delta_v * 1e-3)
            delta = self._solve_branch(temperature, pressure, lo, delta_v, start=rho_ig)
            return self._build(temperature, delta, "vapor")
        lo = rho_ig * 1e-3
        hi = self._upper_bracket(temperature, pressure, rho_ig)
        delta = self._solve_branch(temperature, pressure, lo, hi, start=rho_ig)
        return self._build(temperature, delta, self._phase(temperature, pressure, None))

    def state_rho_t(self, density_kg_m3: float, temperature_k: float) -> FluidState:
        """Return the state at mass density ``density_kg_m3`` and temperature.

        Parameters
        ----------
        density_kg_m3 : float
            Positive mass density.
        temperature_k : float
            Temperature in the validity range.

        Returns
        -------
        FluidState
            Single-phase state.

        Raises
        ------
        DomainError
            Outside the validity range (including the resulting pressure),
            or inside the two-phase region below ``T_c``.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eqs. (56)-(64), p. 1385,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        temperature = self._temperature(temperature_k)
        density = _positive(density_kg_m3, "density_kg_m3")
        m = self.model
        delta = density / m.molar_mass_kg_mol / m.rho_reduce_mol_m3
        if delta >= m.residual.reduced_density_limit:
            raise DomainError(f"{m.name}: density exceeds the model's limit.",
                              field="density_kg_m3")
        p_sat = None
        if temperature < m.t_crit_k:
            delta_l, delta_v, *_ = self._saturation_deltas(temperature)
            if delta_v < delta < delta_l:
                raise DomainError(
                    f"{m.name}: rho = {density} kg/m3 at T = {temperature} K is inside the "
                    "two-phase region; use saturation_t.", field="density_kg_m3")
            p_sat, _ = self._p_of_delta(temperature, delta_v)
        pressure, _ = self._p_of_delta(temperature, delta)
        self._pressure_value(pressure)
        phase = self._phase(temperature, pressure, p_sat)
        return self._build(temperature, delta, phase)

    def departures(self, temperature_k: float, pressure_pa: float) -> DepartureFunctions:
        """Return departure (residual) functions at ``(T, p)``.

        Parameters
        ----------
        temperature_k, pressure_pa : float
            State inside the validity range.

        Returns
        -------
        DepartureFunctions
            ``h - h_ig``, ``s - s_ig``, ``ln(phi)`` and ``cv - cv_ig``.

        Notes
        -----
        ``h - h_ig = R T (tau alphar_tau + delta alphar_delta)``,
        ``s - s_ig(T, p) = R (tau alphar_tau - alphar + ln Z)``,
        ``ln(phi) = alphar + delta alphar_delta - ln Z`` [1]_ and
        ``cv - cv_ig = -R tau^2 alphar_tau_tau``; they follow from Eqs.
        (57)-(62) and (65) of [1]_ by removing the ideal-gas part at equal
        ``T`` and ``p``.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eqs. (57)-(65), p. 1385,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        state = self.state_tp(temperature_k, pressure_pa)
        m = self.model
        r = m.gas_constant_j_mol_k
        tau = m.t_reduce_k / state.temperature_k
        delta = state.molar_density_mol_m3 / m.rho_reduce_mol_m3
        ar, ar_d, _ar_dd, ar_t, ar_tt, _ar_dt = m.residual.derivatives(tau, delta)
        z = state.compressibility_factor
        ln_phi = ar + delta * ar_d - log(z)
        rt = r * state.temperature_k
        return DepartureFunctions(
            temperature_k=state.temperature_k, pressure_pa=state.pressure_pa,
            compressibility_factor=z, enthalpy_j_mol=rt * (tau * ar_t + delta * ar_d),
            entropy_j_mol_k=r * (tau * ar_t - ar + log(z)), gibbs_energy_j_mol=rt * ln_phi,
            ln_fugacity_coefficient=ln_phi, cv_j_mol_k=-r * tau * tau * ar_tt)

    def pressure_unchecked(self, density_kg_m3: float, temperature_k: float) -> float:
        """Return ``p(rho, T)`` without phase or range checks (diagnostics).

        Parameters
        ----------
        density_kg_m3, temperature_k : float
            Positive density and temperature.

        Returns
        -------
        float
            Pressure in Pa from Eq. (56) of [1]_, possibly metastable or
            unstable.

        References
        ----------
        .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
               (2000), Eq. (56), p. 1385,
               https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
        """

        m = self.model
        delta = _positive(density_kg_m3, "density_kg_m3") / m.molar_mass_kg_mol
        return self._p_of_delta(_positive(temperature_k, "temperature_k"),
                                delta / m.rho_reduce_mol_m3)[0]


def with_reference_state(ideal: IdealHelmholtz, *, gas_constant_j_mol_k: float,
                         t_reduce_k: float, rho_reduce_mol_m3: float, temperature_k: float,
                         pressure_pa: float, enthalpy_j_mol: float,
                         entropy_j_mol_k: float) -> IdealHelmholtz:
    """Return ``ideal`` with ``a1, a2`` set to an ideal-gas reference state.

    Parameters
    ----------
    ideal : IdealHelmholtz
        Ideal part whose ``a1`` and ``a2`` are replaced.
    gas_constant_j_mol_k, t_reduce_k, rho_reduce_mol_m3 : float
        Constants of the model the ideal part belongs to.
    temperature_k, pressure_pa : float
        Reference temperature and pressure of the ideal gas.
    enthalpy_j_mol, entropy_j_mol_k : float
        Ideal-gas enthalpy and entropy prescribed there.

    Returns
    -------
    IdealHelmholtz
        Shifted copy. ``a2`` moves enthalpy by ``R T_r a2`` and ``a1``
        moves entropy by ``-R a1`` (Eqs. 58-60 of [1]_).

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), Eqs. (58)-(60), p. 1385,
           https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    base = IdealHelmholtz(0.0, 0.0, ideal.log_tau, ideal.power, ideal.planck_einstein,
                          ideal.two_level)
    tau = t_reduce_k / temperature_k
    delta = pressure_pa / (gas_constant_j_mol_k * temperature_k) / rho_reduce_mol_m3
    a0, a0_t, _ = base.evaluate(tau, delta)
    h0 = gas_constant_j_mol_k * temperature_k * (tau * a0_t + 1.0)
    s0 = gas_constant_j_mol_k * (tau * a0_t - a0)
    a2 = (enthalpy_j_mol - h0) / (gas_constant_j_mol_k * t_reduce_k)
    a1 = -(entropy_j_mol_k - s0) / gas_constant_j_mol_k
    return IdealHelmholtz(a1, a2, ideal.log_tau, ideal.power, ideal.planck_einstein,
                          ideal.two_level)
