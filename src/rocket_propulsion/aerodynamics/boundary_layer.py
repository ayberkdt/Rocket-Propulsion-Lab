"""Flat-plate boundary layers: Blasius, reference temperature and Van Driest II.

Laminar incompressible results come from a numerical shooting solution of the
Blasius equation ``2 f''' + f f'' = 0`` [1]_ (the textbook quotes the rounded
wall value ``f''(0) = 0.332``). Compressible laminar and turbulent skin friction
use the reference-temperature method in the two forms printed by Anderson
(Eckert/Rubesin-Johnson Eq. (18.53) and the Meador-Smart laminar/turbulent
fits) [1]_. Turbulent incompressible friction uses Anderson Eqs. (19.1)-(19.2)
and the local forms quoted by Tauber [2]_; compressible turbulent friction can
also be computed with the Van Driest (II) transformation exactly as written
by Hopkins [3]_ (Karman-Schoenherr incompressible law, Keyes viscosity).
Heating follows from Reynolds analogy (Anderson Eq. (18.50) = Tauber's
Colburn form, Eq. (52)).

There is no transition correlation: the transition Reynolds number is always
an explicit input. Calorically perfect gas, constant Prandtl number and zero
pressure gradient are assumed throughout; no real-gas, catalytic, roughness
or pressure-gradient effects are modelled.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson, Fundamentals of Aerodynamics, 6th
       ed., McGraw-Hill, 2017, Chs. 15, 18 and 19 (printed pp. 935-937,
       1015-1037, 1052-1057).
.. [2] TAUBER-1989-TP2914. M. E. Tauber, NASA TP-2914, 1989, printed pp. 12-13,
       Eqs. (52)-(55), https://ntrs.nasa.gov/citations/19890017745
.. [3] HOPKINS-1972-TND6945. E. J. Hopkins, NASA TN D-6945, 1972, printed
       pp. 2-3, Eqs. (1)-(18), https://ntrs.nasa.gov/citations/19730001588
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import asin, isfinite, log10, sqrt

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import numeric_range_guard, real

__all__ = [
    "ANDERSON_AIR_SUTHERLAND",
    "ANDERSON_RULE_OF_THUMB_TRANSITION_REYNOLDS",
    "BlasiusSolution",
    "CompressibleFlatPlateResult",
    "IncompressibleLaminarFlatPlate",
    "IncompressibleTurbulentFlatPlate",
    "SutherlandViscosity",
    "VanDriestIIResult",
    "adiabatic_wall_temperature",
    "compressible_flat_plate_reference_temperature",
    "flat_plate_regime",
    "laminar_flat_plate_incompressible",
    "laminar_recovery_factor",
    "reference_temperature_ratio",
    "reynolds_analogy_stanton",
    "solve_blasius",
    "stagnation_temperature_ratio",
    "turbulent_flat_plate_incompressible",
    "turbulent_recovery_factor",
    "van_driest_ii",
]


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _nonnegative(value: object, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must not be negative.", field=label)
    return number


# ---------------------------------------------------------------------------
# Transport property and recovery relations
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SutherlandViscosity:
    """Sutherland viscosity law ``mu/mu0 = (T/T0)^1.5 (T0 + S)/(T + S)``.

    Attributes
    ----------
    reference_viscosity_pa_s : float
        ``mu0`` at ``reference_temperature_k``.
    reference_temperature_k : float
        ``T0`` (a reference temperature, not the total temperature).
    sutherland_constant_k : float
        ``S``; Anderson prints 110 K for air.

    Notes
    -----
    Anderson Eq. (15.3) with the sea-level reference ``mu0 = 1.7894e-5`` kg/(m s)
    at ``T0 = 288.16`` K [1]_; the book states the law is approximate and does
    not hold at high temperatures. The worked examples of Chs. 18-19 use
    ``T0 = 288`` K, which callers may pass explicitly.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (15.3) and reference values, printed
           pp. 936-937.
    """

    reference_viscosity_pa_s: float = 1.7894e-5
    reference_temperature_k: float = 288.16
    sutherland_constant_k: float = 110.0

    def __post_init__(self) -> None:
        _positive(self.reference_viscosity_pa_s, "reference_viscosity_pa_s")
        _positive(self.reference_temperature_k, "reference_temperature_k")
        _nonnegative(self.sutherland_constant_k, "sutherland_constant_k")

    @numeric_range_guard("temperature_k")
    def viscosity(self, temperature_k: float) -> float:
        """Return the dynamic viscosity in Pa s at ``temperature_k``.

        Parameters
        ----------
        temperature_k : float
            Static temperature in K, positive.

        Returns
        -------
        float
            Dynamic viscosity in Pa s.

        References
        ----------
        .. [1] ANDERSON-FOA-2017. Eq. (15.3), printed p. 936.
        """

        temperature = _positive(temperature_k, "temperature_k")
        t0 = self.reference_temperature_k
        s = self.sutherland_constant_k
        return self.reference_viscosity_pa_s * (temperature / t0) ** 1.5 * (t0 + s) / (temperature + s)


ANDERSON_AIR_SUTHERLAND = SutherlandViscosity()
ANDERSON_RULE_OF_THUMB_TRANSITION_REYNOLDS = 5.0e5
"""Anderson's "rule of thumb" critical Reynolds number (printed p. 935).

Supplied only as a documented value a caller may choose; no function uses it
as a default.
"""


def laminar_recovery_factor(prandtl: float) -> float:
    """Return the laminar flat-plate recovery factor ``r = Pr^(1/2)``.

    Parameters
    ----------
    prandtl : float
        Constant Prandtl number, positive.

    Returns
    -------
    float
        Recovery factor.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (18.47), printed p. 1030.
    """

    return sqrt(_positive(prandtl, "prandtl"))


def turbulent_recovery_factor(prandtl: float) -> float:
    """Return the approximate turbulent recovery factor ``r = Pr^(1/3)``.

    Parameters
    ----------
    prandtl : float
        Constant Prandtl number, positive.

    Returns
    -------
    float
        Recovery factor.

    Notes
    -----
    Anderson states this as an approximation in Example 19.3 [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Example 19.3, printed p. 1056.
    """

    return _positive(prandtl, "prandtl") ** (1.0 / 3.0)


def stagnation_temperature_ratio(mach: float, gamma: float = 1.4) -> float:
    """Return ``T0/T = 1 + (gamma - 1) M^2 / 2`` for a calorically perfect gas.

    Parameters
    ----------
    mach : float
        Mach number, non-negative.
    gamma : float, optional
        Ratio of specific heats, greater than 1.

    Returns
    -------
    float
        Total-to-static temperature ratio.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (8.40) as used on printed p. 991.
    """

    m = _nonnegative(mach, "mach")
    g = real(gamma, "gamma")
    if g <= 1.0:
        raise DomainError("gamma must exceed 1.", field="gamma")
    return 1.0 + 0.5 * (g - 1.0) * m * m


def adiabatic_wall_temperature(
    edge_temperature_k: float, mach: float, recovery_factor: float, gamma: float = 1.4
) -> float:
    """Return ``Taw = Te [1 + r (T0/Te - 1)]``.

    Parameters
    ----------
    edge_temperature_k : float
        Boundary-layer edge static temperature in K.
    mach : float
        Edge Mach number.
    recovery_factor : float
        Recovery factor ``r`` in (0, 1.5].
    gamma : float, optional
        Ratio of specific heats.

    Returns
    -------
    float
        Adiabatic wall temperature in K.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.46b), (16.49) and (18.56), printed
           pp. 976 and 1034.
    """

    te = _positive(edge_temperature_k, "edge_temperature_k")
    r = _positive(recovery_factor, "recovery_factor")
    if r > 1.5:
        raise DomainError("recovery_factor above 1.5 is outside this model.", field="recovery_factor")
    return te * (1.0 + r * (stagnation_temperature_ratio(mach, gamma) - 1.0))


def reynolds_analogy_stanton(skin_friction: float, prandtl: float, form: str = "flat_plate") -> float:
    """Return the Stanton number from a skin-friction coefficient.

    Parameters
    ----------
    skin_friction : float
        Local (or matching average) skin-friction coefficient, positive.
    prandtl : float
        Constant Prandtl number.
    form : {"flat_plate", "couette"}
        ``"flat_plate"``: ``CH = (cf/2) Pr^(-2/3)``, Anderson Eq. (18.50), the
        same as Tauber's Colburn form, Eq. (52). ``"couette"``:
        ``CH = (cf/2) Pr^(-1)``, Anderson Eqs. (16.59)/(16.98).

    Returns
    -------
    float
        Stanton number ``CH = qw / (rho_e u_e (h_aw - h_w))``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.59), (16.98), (18.50), printed pp. 978,
           994, 1030.
    .. [2] TAUBER-1989-TP2914. Eq. (52), printed p. 12.
    """

    cf = _positive(skin_friction, "skin_friction")
    pr = _positive(prandtl, "prandtl")
    if form == "flat_plate":
        return 0.5 * cf * pr ** (-2.0 / 3.0)
    if form == "couette":
        return 0.5 * cf / pr
    raise InputError("form must be 'flat_plate' or 'couette'.", field="form")


# ---------------------------------------------------------------------------
# Blasius solution
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BlasiusSolution:
    """Numerical solution of ``2 f''' + f f'' = 0``, ``f(0)=f'(0)=0``, ``f'(inf)=1``.

    Attributes
    ----------
    wall_curvature : float
        ``f''(0)``; Anderson quotes 0.332.
    eta_99 : float
        Similarity coordinate where ``f' = 0.99`` (Anderson rounds to 5.0).
    displacement_coefficient : float
        ``eta - f(eta)`` at the outer edge, ``delta* sqrt(Re_x)/x`` (book 1.72).
    momentum_coefficient : float
        ``integral f'(1 - f') d eta``, ``theta sqrt(Re_x)/x`` (book 0.664).
    eta_max : float
        Outer integration limit used to impose ``f'(eta_max) = 1``.
    steps : int
        Number of RK4 steps.
    edge_residual : float
        ``|f'(eta_max) - 1|`` of the accepted shot.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (18.15), shooting method and values,
           printed pp. 1016-1020.
    """

    wall_curvature: float
    eta_99: float
    displacement_coefficient: float
    momentum_coefficient: float
    eta_max: float
    steps: int
    edge_residual: float


def _blasius_shot(curvature: float, eta_max: float, steps: int, record: bool = False):
    h = eta_max / steps

    def rhs(state):
        f, fp, fpp = state
        return (fp, fpp, -0.5 * f * fpp)

    state = (0.0, 0.0, curvature)
    history: list[tuple[float, tuple[float, ...]]] = [(0.0, state)] if record else []
    for index in range(steps):
        k1 = rhs(state)
        k2 = rhs(tuple(s + 0.5 * h * k for s, k in zip(state, k1)))
        k3 = rhs(tuple(s + 0.5 * h * k for s, k in zip(state, k2)))
        k4 = rhs(tuple(s + h * k for s, k in zip(state, k3)))
        state = tuple(s + h / 6.0 * (a + 2.0 * b + 2.0 * c + d) for s, a, b, c, d in zip(state, k1, k2, k3, k4))
        if record:
            history.append(((index + 1) * h, state))
    return state, history


@lru_cache(maxsize=8)
def _solve_blasius_cached(eta_max: float, steps: int) -> BlasiusSolution:
    lower, upper = 0.1, 1.0
    f_lower = _blasius_shot(lower, eta_max, steps)[0][1] - 1.0
    f_upper = _blasius_shot(upper, eta_max, steps)[0][1] - 1.0
    if f_lower * f_upper > 0.0:
        raise ConvergenceError("Blasius shooting bracket does not enclose f'(inf) = 1.")
    # Secant iterations safeguarded by the bracket (regula falsi, Illinois variant).
    side = 0
    curvature = lower
    residual = f_lower
    for _ in range(100):
        curvature = upper - f_upper * (upper - lower) / (f_upper - f_lower)
        residual = _blasius_shot(curvature, eta_max, steps)[0][1] - 1.0
        if abs(residual) < 1.0e-13:
            break
        if residual * f_upper > 0.0:
            upper, f_upper = curvature, residual
            if side == -1:
                f_lower *= 0.5
            side = -1
        else:
            lower, f_lower = curvature, residual
            if side == 1:
                f_upper *= 0.5
            side = 1
    else:
        raise ConvergenceError("Blasius shooting did not converge.")
    final, history = _blasius_shot(curvature, eta_max, steps, record=True)
    eta_99 = None
    momentum = 0.0
    previous_eta, previous_state = history[0]
    for eta, state in history[1:]:
        a = previous_state[1] * (1.0 - previous_state[1])
        b = state[1] * (1.0 - state[1])
        momentum += 0.5 * (a + b) * (eta - previous_eta)
        if eta_99 is None and previous_state[1] < 0.99 <= state[1]:
            fraction = (0.99 - previous_state[1]) / (state[1] - previous_state[1])
            eta_99 = previous_eta + fraction * (eta - previous_eta)
        previous_eta, previous_state = eta, state
    if eta_99 is None:
        raise ConvergenceError("Blasius profile never reached f' = 0.99.")
    return BlasiusSolution(
        wall_curvature=curvature,
        eta_99=eta_99,
        displacement_coefficient=eta_max - final[0],
        momentum_coefficient=momentum,
        eta_max=eta_max,
        steps=steps,
        edge_residual=abs(final[1] - 1.0),
    )


def solve_blasius(eta_max: float = 10.0, steps: int = 4000) -> BlasiusSolution:
    """Solve the Blasius equation by RK4 shooting on ``f''(0)``.

    Parameters
    ----------
    eta_max : float, optional
        Outer limit (between 6 and 20) where ``f'(eta_max) = 1`` is imposed.
    steps : int, optional
        Number of fixed RK4 steps (at least 200).

    Returns
    -------
    BlasiusSolution
        Wall curvature, thickness coefficients and the shot residual.

    Raises
    ------
    ConvergenceError
        If the safeguarded secant iteration fails.

    Notes
    -----
    This is the shooting procedure Anderson describes for Eq. (18.15) [1]_;
    the momentum integral is evaluated with the trapezoidal rule on the RK4
    grid. The result is cached per ``(eta_max, steps)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (18.15) and printed pp. 1016-1020.
    """

    limit = real(eta_max, "eta_max")
    if not 6.0 <= limit <= 20.0:
        raise DomainError("eta_max must be between 6 and 20.", field="eta_max")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 200:
        raise DomainError("steps must be an integer of at least 200.", field="steps")
    return _solve_blasius_cached(limit, steps)


@dataclass(frozen=True, slots=True)
class IncompressibleLaminarFlatPlate:
    """Blasius flat-plate quantities at one station.

    Attributes
    ----------
    reynolds_x : float
        Local Reynolds number ``rho V x / mu``.
    local_skin_friction : float
        ``cf = 2 f''(0) / sqrt(Re_x)`` (Anderson Eq. (18.20): 0.664/sqrt(Re_x)).
    average_skin_friction : float
        ``Cf = 4 f''(0) / sqrt(Re_x)`` for a plate of length ``x``
        (Anderson Eq. (18.22): 1.328/sqrt(Re_c)).
    thickness_over_x, displacement_thickness_over_x, momentum_thickness_over_x : float
        ``delta/x``, ``delta*/x`` and ``theta/x`` (Eqs. (18.23), (18.26), (18.28)).
    local_stanton : float or None
        ``CH = 0.332 Pr^(-2/3)/sqrt(Re_x)`` form (Eq. (18.54)) when a Prandtl
        number was supplied.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (18.20)-(18.28), (18.54), printed
           pp. 1019-1021, 1034.
    """

    reynolds_x: float
    local_skin_friction: float
    average_skin_friction: float
    thickness_over_x: float
    displacement_thickness_over_x: float
    momentum_thickness_over_x: float
    local_stanton: float | None


def laminar_flat_plate_incompressible(reynolds_x: float, prandtl: float | None = None) -> IncompressibleLaminarFlatPlate:
    """Return incompressible laminar flat-plate friction and thicknesses.

    Parameters
    ----------
    reynolds_x : float
        Reynolds number based on distance from the leading edge, positive.
    prandtl : float, optional
        Prandtl number for the Reynolds-analogy Stanton number.

    Returns
    -------
    IncompressibleLaminarFlatPlate
        Coefficients built from the numerical Blasius solution.

    Notes
    -----
    Coefficients use the computed ``f''(0)``, ``eta_99`` and integrals instead
    of the rounded textbook constants [1]_. No laminar-stability check is made;
    use :func:`flat_plate_regime` with an explicit transition Reynolds number.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (18.19)-(18.28), (18.50), (18.54), printed
           pp. 1019-1021, 1030, 1034.
    """

    re = _positive(reynolds_x, "reynolds_x")
    solution = solve_blasius()
    root = sqrt(re)
    cf = 2.0 * solution.wall_curvature / root
    stanton = None
    if prandtl is not None:
        stanton = reynolds_analogy_stanton(cf, prandtl)
    return IncompressibleLaminarFlatPlate(
        reynolds_x=re,
        local_skin_friction=cf,
        average_skin_friction=2.0 * cf,
        thickness_over_x=solution.eta_99 / root,
        displacement_thickness_over_x=solution.displacement_coefficient / root,
        momentum_thickness_over_x=solution.momentum_coefficient / root,
        local_stanton=stanton,
    )


@dataclass(frozen=True, slots=True)
class IncompressibleTurbulentFlatPlate:
    """Turbulent incompressible flat-plate correlations at one station.

    Attributes
    ----------
    reynolds_x : float
        Reynolds number based on the plate length or station distance.
    law : str
        ``"anderson_power_law"`` or ``"meador_smart"``.
    local_skin_friction : float
        Local ``cf``: ``0.0592/Re^0.2`` (Tauber Eq. (54), ``cf/2 = 0.0296/Re^0.2``)
        or ``0.02296/Re^0.139`` (Meador-Smart).
    average_skin_friction : float
        ``Cf``: ``0.074/Re^0.2`` (Anderson Eq. (19.2)) or ``0.02667/Re^0.139``.
    thickness_over_x : float or None
        ``0.37/Re^0.2`` (Anderson Eq. (19.1)); ``None`` for Meador-Smart.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (19.1)-(19.2) and Sec. 19.2.2, printed
           pp. 1052-1053, 1056.
    .. [2] TAUBER-1989-TP2914. Eq. (54), printed p. 13.
    """

    reynolds_x: float
    law: str
    local_skin_friction: float
    average_skin_friction: float
    thickness_over_x: float | None


def turbulent_flat_plate_incompressible(reynolds_x: float, law: str = "anderson_power_law") -> IncompressibleTurbulentFlatPlate:
    """Return incompressible turbulent flat-plate skin friction.

    Parameters
    ----------
    reynolds_x : float
        Reynolds number, positive. Tauber states the 1/5-power local law was
        used up to ``Re = 1e7`` and the Schultz-Grunow form above it [2]_; a
        warning-free result is returned regardless, the caller owns the range.
    law : {"anderson_power_law", "meador_smart"}
        Correlation family.

    Returns
    -------
    IncompressibleTurbulentFlatPlate
        Local and average skin friction (and thickness for the power law).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (19.1)-(19.2) and Meador-Smart turbulent
           relations, printed pp. 1052-1053, 1056.
    .. [2] TAUBER-1989-TP2914. Eqs. (54)-(55), printed p. 13.
    """

    re = _positive(reynolds_x, "reynolds_x")
    if law == "anderson_power_law":
        return IncompressibleTurbulentFlatPlate(re, law, 2.0 * 0.0296 / re**0.2, 0.074 / re**0.2, 0.37 / re**0.2)
    if law == "meador_smart":
        return IncompressibleTurbulentFlatPlate(re, law, 0.02296 / re**0.139, 0.02667 / re**0.139, None)
    raise InputError("law must be 'anderson_power_law' or 'meador_smart'.", field="law")


def flat_plate_regime(reynolds_x: float, transition_reynolds: float) -> str:
    """Classify a flat-plate station as laminar or turbulent.

    Parameters
    ----------
    reynolds_x : float
        Local Reynolds number.
    transition_reynolds : float
        Caller-chosen transition Reynolds number (no default: Anderson calls
        ``5e5`` only a rule of thumb, and Tauber's review shows measured
        transition Reynolds numbers spanning more than an order of magnitude).

    Returns
    -------
    str
        ``"laminar"`` if ``reynolds_x < transition_reynolds`` else ``"turbulent"``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Definition of Re_cr, printed p. 935.
    .. [2] TAUBER-1989-TP2914. Sec. 3 and Fig. 9, printed p. 10.
    """

    re = _positive(reynolds_x, "reynolds_x")
    tr = _positive(transition_reynolds, "transition_reynolds")
    return "laminar" if re < tr else "turbulent"


# ---------------------------------------------------------------------------
# Reference temperature method
# ---------------------------------------------------------------------------


def reference_temperature_ratio(
    mach: float,
    wall_to_edge_temperature: float,
    method: str = "eckert",
    recovery_factor: float | None = None,
    gamma: float = 1.4,
) -> float:
    """Return ``T*/Te`` for the reference-temperature method.

    Parameters
    ----------
    mach : float
        Edge Mach number.
    wall_to_edge_temperature : float
        ``Tw/Te``.
    method : {"eckert", "meador_smart_laminar", "meador_smart_turbulent"}
        ``"eckert"``: ``1 + 0.032 Me^2 + 0.58 (Tw/Te - 1)``, Anderson
        Eq. (18.53). Meador-Smart laminar:
        ``0.45 + 0.55 Tw/Te + 0.16 r (gamma-1)/2 Me^2``; turbulent:
        ``0.5 (1 + Tw/Te) + 0.16 r (gamma-1)/2 Me^2`` (both as printed by
        Anderson).
    recovery_factor : float, optional
        Required for the Meador-Smart forms.
    gamma : float, optional
        Ratio of specific heats.

    Returns
    -------
    float
        Reference-to-edge temperature ratio.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (18.53), Secs. 18.4.1 and 19.2.2, printed
           pp. 1033, 1036, 1056.
    """

    m = _nonnegative(mach, "mach")
    ratio = _positive(wall_to_edge_temperature, "wall_to_edge_temperature")
    if method == "eckert":
        return 1.0 + 0.032 * m * m + 0.58 * (ratio - 1.0)
    if method in {"meador_smart_laminar", "meador_smart_turbulent"}:
        if recovery_factor is None:
            raise InputError("Meador-Smart reference temperature needs recovery_factor.", field="recovery_factor")
        r = _positive(recovery_factor, "recovery_factor")
        g = real(gamma, "gamma")
        if g <= 1.0:
            raise DomainError("gamma must exceed 1.", field="gamma")
        compressibility = 0.16 * r * 0.5 * (g - 1.0) * m * m
        if method == "meador_smart_laminar":
            return 0.45 + 0.55 * ratio + compressibility
        return 0.5 * (1.0 + ratio) + compressibility
    raise InputError(
        "method must be 'eckert', 'meador_smart_laminar' or 'meador_smart_turbulent'.", field="method"
    )


@dataclass(frozen=True, slots=True)
class CompressibleFlatPlateResult:
    """Reference-temperature flat-plate result at distance ``x`` (one side).

    Attributes
    ----------
    regime : str
        ``"laminar"`` or ``"turbulent"`` (chosen by the caller).
    method : str
        Reference-temperature variant.
    mach : float
        Edge Mach number.
    recovery_factor : float
        Recovery factor used for ``Taw``.
    adiabatic_wall_temperature_k, wall_temperature_k, reference_temperature_k : float
        ``Taw``, ``Tw`` and ``T*`` in K.
    reference_density_kg_m3, reference_viscosity_pa_s : float
        ``rho* = p/(R T*)`` and ``mu*`` from the viscosity law.
    reference_reynolds : float
        ``Re*_x = rho* u_e x / mu*``.
    local_skin_friction, average_skin_friction : float
        ``cf*`` and ``Cf*`` based on ``rho*`` (Anderson Eqs. (18.52), (18.57),
        (19.3)/(19.4)).
    wall_shear_pa : float
        ``tau_w = cf* rho* u_e^2 / 2`` at ``x``.
    friction_drag_per_area_pa : float
        ``Cf* rho* u_e^2 / 2``: mean shear over ``[0, x]`` (one side).
    reference_stanton : float
        ``CH*`` from Reynolds analogy, Eq. (18.50)/(18.55).
    wall_heat_flux_w_m2 : float
        ``CH* rho* u_e cp (Taw - Tw)``; positive into the wall.
    warnings : tuple of WarningMessage
        Model-applicability notes.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (18.50)-(18.57), (19.3)-(19.4), printed
           pp. 1030-1036, 1054.
    """

    regime: str
    method: str
    mach: float
    recovery_factor: float
    adiabatic_wall_temperature_k: float
    wall_temperature_k: float
    reference_temperature_k: float
    reference_density_kg_m3: float
    reference_viscosity_pa_s: float
    reference_reynolds: float
    local_skin_friction: float
    average_skin_friction: float
    wall_shear_pa: float
    friction_drag_per_area_pa: float
    reference_stanton: float
    wall_heat_flux_w_m2: float
    warnings: tuple[WarningMessage, ...]


def compressible_flat_plate_reference_temperature(
    *,
    edge_temperature_k: float,
    edge_pressure_pa: float,
    edge_velocity_m_s: float,
    x_m: float,
    gas_constant_j_kg_k: float,
    gamma: float = 1.4,
    prandtl: float = 0.71,
    wall_temperature_k: float | None = None,
    regime: str = "laminar",
    method: str = "eckert",
    recovery_factor: float | None = None,
    viscosity: SutherlandViscosity = ANDERSON_AIR_SUTHERLAND,
) -> CompressibleFlatPlateResult:
    """Return compressible flat-plate friction and heating by reference temperature.

    Parameters
    ----------
    edge_temperature_k, edge_pressure_pa, edge_velocity_m_s : float
        Edge (free-stream for a flat plate at zero incidence) state.
    x_m : float
        Distance from the leading edge (or plate length for averages), m.
    gas_constant_j_kg_k : float
        Specific gas constant ``R``.
    gamma, prandtl : float, optional
        Ratio of specific heats and constant Prandtl number.
    wall_temperature_k : float, optional
        Wall temperature; ``None`` selects an adiabatic wall (``Tw = Taw``).
    regime : {"laminar", "turbulent"}
        Incompressible law used at ``T*``. Laminar: Blasius ``cf* = 2 f''(0)
        / sqrt(Re*)`` (book 0.664); turbulent: ``Cf* = 0.074/Re*^0.2``
        (Eq. (19.3)) with ``method="eckert"`` or ``0.02667/Re*^0.139`` with
        ``method="meador_smart"``.
    method : {"eckert", "meador_smart"}
        Reference-temperature variant (see :func:`reference_temperature_ratio`).
    recovery_factor : float, optional
        Defaults to ``sqrt(Pr)`` (laminar) or ``Pr^(1/3)`` (turbulent).
        Anderson's Example 19.2 deliberately reuses the laminar value.
    viscosity : SutherlandViscosity, optional
        Viscosity law for ``mu*``.

    Returns
    -------
    CompressibleFlatPlateResult
        Reference state, friction, heating and warnings.

    Notes
    -----
    Reynolds analogy at the reference state, Eq. (18.55), is applied for both
    regimes; for turbulent flow this is the Colburn form that Tauber reports
    gave the better agreement with data (TP-2914 Eq. (52)). Calorically
    perfect gas, ``cp = gamma R/(gamma - 1)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (18.47)-(18.57), (19.3)-(19.4), Secs.
           18.4.1 and 19.2.2, Examples 18.2, 18.3, 19.2, 19.3, printed
           pp. 1030-1037, 1054-1057.
    .. [2] TAUBER-1989-TP2914. Eq. (52), printed p. 12.
    """

    te = _positive(edge_temperature_k, "edge_temperature_k")
    pe = _positive(edge_pressure_pa, "edge_pressure_pa")
    ue = _positive(edge_velocity_m_s, "edge_velocity_m_s")
    x = _positive(x_m, "x_m")
    gas = _positive(gas_constant_j_kg_k, "gas_constant_j_kg_k")
    g = real(gamma, "gamma")
    if g <= 1.0:
        raise DomainError("gamma must exceed 1.", field="gamma")
    pr = _positive(prandtl, "prandtl")
    if regime not in {"laminar", "turbulent"}:
        raise InputError("regime must be 'laminar' or 'turbulent'.", field="regime")
    if method not in {"eckert", "meador_smart"}:
        raise InputError("method must be 'eckert' or 'meador_smart'.", field="method")
    mach = ue / sqrt(g * gas * te)
    if recovery_factor is None:
        r = laminar_recovery_factor(pr) if regime == "laminar" else turbulent_recovery_factor(pr)
    else:
        r = _positive(recovery_factor, "recovery_factor")
    taw = adiabatic_wall_temperature(te, mach, r, g)
    tw = taw if wall_temperature_k is None else _positive(wall_temperature_k, "wall_temperature_k")
    variant = "eckert" if method == "eckert" else f"meador_smart_{regime}"
    t_star = te * reference_temperature_ratio(mach, tw / te, variant, r, g)
    rho_star = pe / (gas * t_star)
    mu_star = viscosity.viscosity(t_star)
    re_star = rho_star * ue * x / mu_star
    if regime == "laminar":
        curvature = solve_blasius().wall_curvature
        cf = 2.0 * curvature / sqrt(re_star)
        cf_average = 2.0 * cf
    elif method == "eckert":
        cf = 2.0 * 0.0296 / re_star**0.2
        cf_average = 0.074 / re_star**0.2
    else:
        cf = 0.02296 / re_star**0.139
        cf_average = 0.02667 / re_star**0.139
    dynamic = 0.5 * rho_star * ue * ue
    stanton = reynolds_analogy_stanton(cf, pr)
    cp = g * gas / (g - 1.0)
    warnings: list[WarningMessage] = []
    if t_star > 1000.0 or taw > 1000.0:
        warnings.append(WarningMessage(
            "calorically_perfect_limit",
            "Temperatures above about 1000 K: Anderson notes constant cp and Pr are justified only up to "
            "roughly 1000 K; real-gas effects are not modelled.",
        ))
    if regime == "turbulent":
        warnings.append(WarningMessage(
            "turbulent_reference_temperature_uncertain",
            "Anderson's Example 19.2 shows a 20 percent spread between turbulent methods; treat as an estimate.",
            severity="info",
        ))
    return CompressibleFlatPlateResult(
        regime=regime,
        method=method,
        mach=mach,
        recovery_factor=r,
        adiabatic_wall_temperature_k=taw,
        wall_temperature_k=tw,
        reference_temperature_k=t_star,
        reference_density_kg_m3=rho_star,
        reference_viscosity_pa_s=mu_star,
        reference_reynolds=re_star,
        local_skin_friction=cf,
        average_skin_friction=cf_average,
        wall_shear_pa=cf * dynamic,
        friction_drag_per_area_pa=cf_average * dynamic,
        reference_stanton=stanton,
        wall_heat_flux_w_m2=stanton * rho_star * ue * cp * (taw - tw),
        warnings=tuple(warnings),
    )


# ---------------------------------------------------------------------------
# Van Driest II
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VanDriestIIResult:
    """Van Driest (II) turbulent flat-plate skin friction.

    Attributes
    ----------
    mach, edge_temperature_k, wall_temperature_k, adiabatic_wall_temperature_k : float
        Inputs and ``Taw = Te (1 + r m)``, ``m = 0.2 Me^2`` (Eq. (17)/(18)).
    reynolds_x : float
        Length Reynolds number at edge conditions.
    f_c, f_theta, f_x : float
        Transformation functions ``F_c``, ``F_theta = mu_e/mu_w`` and
        ``F_x = F_theta/F_c`` (Eqs. (8)-(11)).
    local_skin_friction, average_skin_friction : float
        Compressible ``cf`` and ``CF`` (Eqs. (4), (5)).
    momentum_reynolds : float
        ``Re_theta = CF Re_x / 2`` (Eq. (3) inverted on the compressible side).
    incompressible_average_skin_friction : float
        ``CF_bar`` from the Karman-Schoenherr law at ``Re_x_bar = F_x Re_x``.

    References
    ----------
    .. [1] HOPKINS-1972-TND6945. Eqs. (1)-(18), printed pp. 2-3.
    """

    mach: float
    edge_temperature_k: float
    wall_temperature_k: float
    adiabatic_wall_temperature_k: float
    reynolds_x: float
    f_c: float
    f_theta: float
    f_x: float
    local_skin_friction: float
    average_skin_friction: float
    momentum_reynolds: float
    incompressible_average_skin_friction: float


def _karman_schoenherr_average(reynolds: float) -> float:
    # 0.242/sqrt(C) = log10(Re C); the left side minus the right side decreases in C.
    lower, upper = 1.0e-7, 0.2

    def residual(c: float) -> float:
        return 0.242 / sqrt(c) - log10(reynolds * c)

    if residual(lower) < 0.0 or residual(upper) > 0.0:
        raise DomainError("Reynolds number outside the Karman-Schoenherr solution bracket.", field="reynolds_x")
    for _ in range(200):
        middle = sqrt(lower * upper)
        if residual(middle) > 0.0:
            lower = middle
        else:
            upper = middle
        if upper / lower - 1.0 < 1.0e-14:
            break
    return sqrt(lower * upper)


def _keyes_viscosity_ratio(edge_temperature_k: float, wall_temperature_k: float) -> float:
    te, tw = edge_temperature_k, wall_temperature_k
    return sqrt(te / tw) * (1.0 + 122.0 / tw * 10.0 ** (-5.0 / tw)) / (1.0 + 122.0 / te * 10.0 ** (-5.0 / te))


@numeric_range_guard("mach")
def van_driest_ii(
    *,
    mach: float,
    edge_temperature_k: float,
    reynolds_x: float,
    wall_temperature_k: float | None = None,
    wall_to_adiabatic_temperature: float | None = None,
    recovery_factor: float = 0.88,
    viscosity_law: str = "keyes",
    viscosity_exponent: float = 0.76,
) -> VanDriestIIResult:
    """Return Van Driest (II) turbulent flat-plate skin friction.

    Parameters
    ----------
    mach : float
        Edge Mach number, ``0 <= Me <= 10`` (the chart range of [1]_).
    edge_temperature_k : float
        Edge static temperature in K.
    reynolds_x : float
        Length Reynolds number at edge conditions; Hopkins notes the
        Karman-Schoenherr law is corroborated by data for about
        ``3e5 <= Re <= 4.5e8`` and the charts span ``1e5 <= Re_x <= 1e9``.
    wall_temperature_k : float, optional
        Wall temperature. Give at most one of this and
        ``wall_to_adiabatic_temperature``; giving neither means adiabatic.
    wall_to_adiabatic_temperature : float, optional
        ``Tw/Taw`` with ``Taw`` from Eq. (18) and ``recovery_factor``.
    recovery_factor : float, optional
        Hopkins recommends 0.88 (Van Driest's nomogram assumed 1.0).
    viscosity_law : {"keyes", "power"}
        ``F_theta = mu_e/mu_w`` from the Keyes formula, Eq. (10), or from
        ``mu ~ T^viscosity_exponent``.
    viscosity_exponent : float, optional
        Exponent for ``viscosity_law="power"``.

    Returns
    -------
    VanDriestIIResult
        Transformation functions and compressible skin-friction coefficients.

    Notes
    -----
    ``m = 0.2 Me^2`` is Hopkins' Eq. (17) and therefore fixes ``gamma = 1.4``.
    The incompressible law is Karman-Schoenherr, Eq. (1), with the local value
    from Eq. (2), both evaluated at ``Re_x_bar = F_x Re_x``; then
    ``cf = cf_bar/F_c``. For ``Me = 0`` Eq. (9) replaces Eq. (8).

    References
    ----------
    .. [1] HOPKINS-1972-TND6945. Eqs. (1)-(18), printed pp. 2-3.
    """

    me = _nonnegative(mach, "mach")
    if me > 10.0:
        raise DomainError("Van Driest II charts of TN D-6945 cover Me <= 10.", field="mach")
    te = _positive(edge_temperature_k, "edge_temperature_k")
    re = _positive(reynolds_x, "reynolds_x")
    r = _positive(recovery_factor, "recovery_factor")
    if r > 1.0:
        raise DomainError("recovery_factor above 1 is outside this model.", field="recovery_factor")
    m = 0.2 * me * me
    taw = te * (1.0 + r * m)
    if wall_temperature_k is not None and wall_to_adiabatic_temperature is not None:
        raise InputError("Give wall_temperature_k or wall_to_adiabatic_temperature, not both.",
                         field="wall_temperature_k")
    if wall_temperature_k is not None:
        tw = _positive(wall_temperature_k, "wall_temperature_k")
    elif wall_to_adiabatic_temperature is not None:
        tw = _positive(wall_to_adiabatic_temperature, "wall_to_adiabatic_temperature") * taw
    else:
        tw = taw
    f = tw / te
    if m == 0.0:
        f_c = ((1.0 + sqrt(f)) / 2.0) ** 2
    else:
        a_coef = sqrt(r * m / f)
        b_coef = (1.0 + r * m - f) / f
        denominator = sqrt(4.0 * a_coef * a_coef + b_coef * b_coef)
        alpha = (2.0 * a_coef * a_coef - b_coef) / denominator
        beta = b_coef / denominator
        if not (-1.0 <= alpha <= 1.0 and -1.0 <= beta <= 1.0):
            raise DomainError("Van Driest II arcsine arguments out of range.", field="wall_temperature_k")
        f_c = r * m / (asin(alpha) + asin(beta)) ** 2
    if viscosity_law == "keyes":
        f_theta = _keyes_viscosity_ratio(te, tw)
    elif viscosity_law == "power":
        exponent = _positive(viscosity_exponent, "viscosity_exponent")
        f_theta = (te / tw) ** exponent
    else:
        raise InputError("viscosity_law must be 'keyes' or 'power'.", field="viscosity_law")
    f_x = f_theta / f_c
    cf_bar_average = _karman_schoenherr_average(f_x * re)
    cf_bar_local = 0.242 * cf_bar_average / (0.242 + 0.8686 * sqrt(cf_bar_average))
    average = cf_bar_average / f_c
    local = cf_bar_local / f_c
    if not (isfinite(average) and isfinite(local)):
        raise ConvergenceError("Van Driest II produced a non-finite coefficient.")
    return VanDriestIIResult(
        mach=me,
        edge_temperature_k=te,
        wall_temperature_k=tw,
        adiabatic_wall_temperature_k=taw,
        reynolds_x=re,
        f_c=f_c,
        f_theta=f_theta,
        f_x=f_x,
        local_skin_friction=local,
        average_skin_friction=average,
        momentum_reynolds=0.5 * average * re,
        incompressible_average_skin_friction=cf_bar_average,
    )

