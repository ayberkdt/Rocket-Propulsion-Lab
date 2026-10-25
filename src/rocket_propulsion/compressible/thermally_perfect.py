"""Isentropic flow, normal and oblique shocks of a thermally perfect gas.

A thermally perfect gas obeys ``p = rho R T`` but its specific heats depend
on temperature (a calorically imperfect gas in the words of NACA Report 1135
[1]_). Two caloric models are provided:

* :class:`HarmonicVibratorGas` - the NACA 1135 model of a diatomic gas whose
  vibrational heat capacity is that of a simple harmonic vibrator with
  characteristic temperature ``Theta`` (5500 degR for air) [1]_.
* :class:`PolynomialGas` - any repository species or mixture evaluated from
  NASA Glenn polynomials [2]_ (``rocket_propulsion.thermochemistry``).

The flow relations are the conservation laws written with ``h(T)`` and the
entropy function ``phi(T) = integral cp dT / T``; NACA 1135 gives them in
closed form for the harmonic vibrator, equations (180)-(198) [1]_.
Dissociation and other gaseous imperfections are not modelled; NACA 1135
places their onset near 5000 degR [1]_.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953, "Imperfect-gas effects",
       pp. 626-630, equations (175)-(200), charts 9-24.
.. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, degrees, exp, expm1, isfinite, log, pi, radians, sin, sqrt, tan
from typing import Any, Protocol

from rocket_propulsion.core.constants import AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
from rocket_propulsion.core.errors import DomainError

from ._advanced_common import bracketed_root, maximize_unimodal

#: Vibrational characteristic temperature used by NACA 1135 for air,
#: ``Theta = 5500 degR`` (p. 627), converted with ``1 degR = 5/9 K``.
NACA1135_AIR_VIBRATIONAL_TEMPERATURE_K = 5500.0 * 5.0 / 9.0


class ThermallyPerfectGas(Protocol):
    """Caloric model accepted by the thermally perfect flow functions.

    Attributes
    ----------
    gas_constant_j_kg_k : float
        Specific gas constant ``R`` in J/(kg K).
    minimum_temperature_k, maximum_temperature_k : float
        Temperature range of the caloric model.

    Notes
    -----
    ``enthalpy_j_kg`` and ``entropy_function_j_kg_k`` may carry arbitrary
    additive constants; only differences are used.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953.
    """

    gas_constant_j_kg_k: float

    @property
    def minimum_temperature_k(self) -> float:  # pragma: no cover
        """Lowest temperature of the caloric model, K."""
        ...

    @property
    def maximum_temperature_k(self) -> float:  # pragma: no cover
        """Highest temperature of the caloric model, K."""
        ...

    def cp_j_kg_k(self, temperature_k: float) -> float:  # pragma: no cover
        """Specific heat at constant pressure, J/(kg K)."""
        ...

    def enthalpy_j_kg(self, temperature_k: float) -> float:  # pragma: no cover
        """Specific enthalpy (arbitrary datum), J/kg."""
        ...

    def entropy_function_j_kg_k(self, temperature_k: float) -> float:  # pragma: no cover
        """Entropy function ``integral cp dT/T`` (arbitrary datum), J/(kg K)."""
        ...


@dataclass(frozen=True, slots=True)
class HarmonicVibratorGas:
    """NACA 1135 calorically imperfect, thermally perfect diatomic gas [1]_.

    Attributes
    ----------
    gamma_perfect : float, default 1.4
        Ratio of specific heats with the vibrational mode frozen (``gamma_perf``).
    vibrational_temperature_k : float
        Characteristic vibrational temperature ``Theta`` in K; the default is
        the 5500 degR NACA 1135 value for air.
    gas_constant_j_kg_k : float
        Specific gas constant in J/(kg K). Flow ratios do not depend on it.

    Raises
    ------
    DomainError
        If a parameter is not finite or out of range.

    Notes
    -----
    With ``x = Theta/T`` and ``cp_perf = gamma_perf R/(gamma_perf - 1)``
    (NACA 1135 eqs. (176), (178), (180), (184), (189)) [1]_:

    * ``cp = cp_perf + R x^2 e^x/(e^x - 1)^2``
    * ``gamma = 1 + (gamma_perf - 1)/(1 + (gamma_perf - 1) x^2 e^x/(e^x - 1)^2)``
    * ``h = cp_perf T + R Theta/(e^x - 1)``
    * ``phi = cp_perf ln T + R [x e^x/(e^x - 1) - ln(e^x - 1)]``

    The forms are evaluated with ``exp(-x)`` and ``expm1`` so that they stay
    finite for any positive temperature.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, pp. 627-628.
    """

    gamma_perfect: float = 1.4
    vibrational_temperature_k: float = NACA1135_AIR_VIBRATIONAL_TEMPERATURE_K
    gas_constant_j_kg_k: float = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K

    def __post_init__(self) -> None:
        if not isfinite(self.gamma_perfect) or self.gamma_perfect <= 1.0:
            raise DomainError("gamma_perfect must be finite and greater than 1.",
                              field="gamma_perfect")
        if not isfinite(self.vibrational_temperature_k) or self.vibrational_temperature_k <= 0.0:
            raise DomainError("Vibrational temperature must be finite and positive.",
                              field="vibrational_temperature_k")
        if not isfinite(self.gas_constant_j_kg_k) or self.gas_constant_j_kg_k <= 0.0:
            raise DomainError("Gas constant must be finite and positive.",
                              field="gas_constant_j_kg_k")

    @property
    def minimum_temperature_k(self) -> float:
        """Lowest temperature accepted, K (a tiny positive floor) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953.
        """
        return 1.0e-6 * self.vibrational_temperature_k

    @property
    def maximum_temperature_k(self) -> float:
        """Highest temperature accepted, K (no dissociation is modelled) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953.
        """
        return 1.0e3 * self.vibrational_temperature_k

    def _cp_perfect(self) -> float:
        return self.gamma_perfect * self.gas_constant_j_kg_k / (self.gamma_perfect - 1.0)

    def _vibrational_factor(self, temperature_k: float) -> float:
        x = self.vibrational_temperature_k / temperature_k
        decay = exp(-x)
        return x * x * decay / expm1(-x) ** 2

    def cp_j_kg_k(self, temperature_k: float) -> float:
        """Return ``cp(T)`` in J/(kg K), NACA 1135 eq. (178) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953.
        """
        _check_temperature(self, temperature_k)
        return self._cp_perfect() + self.gas_constant_j_kg_k * self._vibrational_factor(
            temperature_k
        )

    def gamma(self, temperature_k: float) -> float:
        """Return the ratio of specific heats, NACA 1135 eq. (180) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953.
        """
        _check_temperature(self, temperature_k)
        g = self.gamma_perfect
        return 1.0 + (g - 1.0) / (1.0 + (g - 1.0) * self._vibrational_factor(temperature_k))

    def enthalpy_j_kg(self, temperature_k: float) -> float:
        """Return ``h(T)`` in J/kg (zero-point of the vibrator at 0 K) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953, eq. (184).
        """
        _check_temperature(self, temperature_k)
        x = self.vibrational_temperature_k / temperature_k
        vibrational = exp(-x) / -expm1(-x)
        return self._cp_perfect() * temperature_k + (
            self.gas_constant_j_kg_k * self.vibrational_temperature_k * vibrational
        )

    def entropy_function_j_kg_k(self, temperature_k: float) -> float:
        """Return ``phi(T) = integral cp dT/T`` in J/(kg K) [1]_.

        References
        ----------
        .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and
               Charts for Compressible Flow", NACA Report 1135, 1953, eq. (189).
        """
        _check_temperature(self, temperature_k)
        x = self.vibrational_temperature_k / temperature_k
        decay = exp(-x)
        one_minus = -expm1(-x)
        vibrational = x * decay / one_minus - log(one_minus)
        return self._cp_perfect() * log(temperature_k) + self.gas_constant_j_kg_k * vibrational


@dataclass(frozen=True, slots=True)
class PolynomialGas:
    """Thermally perfect gas from a NASA Glenn species or mixture [1]_.

    Attributes
    ----------
    source : object
        A :class:`rocket_propulsion.thermochemistry.Species` or
        :class:`rocket_propulsion.thermochemistry.Mixture` (any object whose
        ``properties(T)`` returns ``cp_j_kg_k``, ``enthalpy_j_kg``,
        ``entropy_j_kg_k`` and ``gas_constant_j_kg_k`` and that exposes
        ``minimum_temperature_k`` and ``maximum_temperature_k``).

    Notes
    -----
    The composition is frozen. The mixture entropy at the 1 bar reference
    pressure is used as ``phi(T)``; its constant mixing term cancels in every
    difference. Temperatures outside the polynomial range raise
    :class:`DomainError` (fail closed).

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    source: Any

    @property
    def gas_constant_j_kg_k(self) -> float:
        """Specific gas constant ``R/M`` in J/(kg K) [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        reference = 0.5 * (self.minimum_temperature_k + self.maximum_temperature_k)
        return self.source.properties(reference).gas_constant_j_kg_k

    @property
    def minimum_temperature_k(self) -> float:
        """Lowest temperature of the polynomial fit, K [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        return float(self.source.minimum_temperature_k)

    @property
    def maximum_temperature_k(self) -> float:
        """Highest temperature of the polynomial fit, K [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        return float(self.source.maximum_temperature_k)

    def cp_j_kg_k(self, temperature_k: float) -> float:
        """Return frozen ``cp(T)`` in J/(kg K) [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        _check_temperature(self, temperature_k)
        return self.source.properties(temperature_k).cp_j_kg_k

    def gamma(self, temperature_k: float) -> float:
        """Return frozen ``cp/cv`` at ``T`` [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        cp = self.cp_j_kg_k(temperature_k)
        return cp / (cp - self.gas_constant_j_kg_k)

    def enthalpy_j_kg(self, temperature_k: float) -> float:
        """Return ``h(T)`` in J/kg [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        _check_temperature(self, temperature_k)
        return self.source.properties(temperature_k).enthalpy_j_kg

    def entropy_function_j_kg_k(self, temperature_k: float) -> float:
        """Return ``s(T, 1 bar)`` in J/(kg K), used as ``phi(T)`` [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        _check_temperature(self, temperature_k)
        return self.source.properties(temperature_k).entropy_j_kg_k


@dataclass(frozen=True, slots=True)
class ThermallyPerfectIsentropicResult:
    """Isentropic state of a thermally perfect gas at one Mach number [1]_.

    Attributes
    ----------
    mach : float
        Local Mach number.
    total_temperature_k, temperature_k : float
        Stagnation and static temperatures, K.
    gamma : float
        Local ratio of specific heats ``gamma(T)``.
    temperature_ratio, pressure_ratio, density_ratio : float
        ``T/Tt``, ``p/pt`` and ``rho/rho_t``.
    dynamic_pressure_ratio : float
        ``q/pt = (gamma/2) M^2 p/pt``.
    area_ratio : float
        ``A/A*`` (stream-tube area relative to the sonic section).
    velocity_critical_ratio : float
        ``V/a*``.
    critical_temperature_k : float
        Static temperature at the sonic point, K.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eqs. (182)-(191).
    """

    mach: float
    total_temperature_k: float
    temperature_k: float
    gamma: float
    temperature_ratio: float
    pressure_ratio: float
    density_ratio: float
    dynamic_pressure_ratio: float
    area_ratio: float
    velocity_critical_ratio: float
    critical_temperature_k: float


@dataclass(frozen=True, slots=True)
class ThermallyPerfectNormalShockResult:
    """Jump across a normal shock in a thermally perfect gas [1]_.

    Attributes
    ----------
    upstream_mach, downstream_mach : float
        Mach numbers before and after the shock.
    upstream_temperature_k, downstream_temperature_k : float
        Static temperatures, K.
    total_temperature_k : float
        Stagnation temperature (unchanged across the shock), K.
    pressure_ratio, density_ratio, temperature_ratio : float
        ``p2/p1``, ``rho2/rho1`` and ``T2/T1``.
    stagnation_pressure_ratio : float
        ``pt2/pt1``.
    upstream_static_to_downstream_total : float
        ``p1/pt2`` (pitot relation).

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eqs. (192)-(194).
    """

    upstream_mach: float
    downstream_mach: float
    upstream_temperature_k: float
    downstream_temperature_k: float
    total_temperature_k: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    stagnation_pressure_ratio: float
    upstream_static_to_downstream_total: float


@dataclass(frozen=True, slots=True)
class ThermallyPerfectObliqueShockResult:
    """Attached weak oblique shock in a thermally perfect gas [1]_.

    Attributes
    ----------
    upstream_mach, downstream_mach : float
        Mach numbers before and after the wave.
    turning_angle_deg, wave_angle_deg : float
        Flow deflection and shock-wave angle, degrees.
    maximum_turning_angle_deg : float
        Largest attached deflection at this upstream state, degrees.
    normal_shock : ThermallyPerfectNormalShockResult
        Jump of the normal velocity component (``M1 sin(beta)``).
    pressure_ratio, density_ratio, temperature_ratio, stagnation_pressure_ratio : float
        Ratios across the wave.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eqs. (195)-(198),
           pp. 629-630.
    """

    upstream_mach: float
    downstream_mach: float
    turning_angle_deg: float
    wave_angle_deg: float
    maximum_turning_angle_deg: float
    normal_shock: ThermallyPerfectNormalShockResult
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    stagnation_pressure_ratio: float


def _check_temperature(gas: Any, temperature_k: float) -> None:
    """Raise :class:`DomainError` if ``T`` is outside the caloric model range."""
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise DomainError("Temperature must be finite and positive.", field="temperature_k")
    if not gas.minimum_temperature_k <= temperature_k <= gas.maximum_temperature_k:
        raise DomainError(
            f"Temperature {temperature_k:.6g} K is outside the gas model range "
            f"[{gas.minimum_temperature_k:.6g}, {gas.maximum_temperature_k:.6g}] K.",
            code="temperature_out_of_model_range",
            field="temperature_k",
        )


def _gamma(gas: ThermallyPerfectGas, temperature_k: float) -> float:
    """Return ``cp/(cp - R)`` at ``T``."""
    cp = gas.cp_j_kg_k(temperature_k)
    return cp / (cp - gas.gas_constant_j_kg_k)


def _sound_speed(gas: ThermallyPerfectGas, temperature_k: float) -> float:
    """Return ``a = sqrt(gamma(T) R T)`` (NACA 1135 eq. (182))."""
    return sqrt(_gamma(gas, temperature_k) * gas.gas_constant_j_kg_k * temperature_k)


def _static_temperature(gas: ThermallyPerfectGas, mach: float, total_temperature_k: float) -> float:
    """Solve ``h(Tt) - h(T) = M^2 a(T)^2 / 2`` for the static temperature."""
    total_enthalpy = gas.enthalpy_j_kg(total_temperature_k)
    scale = gas.cp_j_kg_k(total_temperature_k) * total_temperature_k

    def residual(temperature: float) -> float:
        kinetic = 0.5 * mach * mach * _gamma(gas, temperature) * gas.gas_constant_j_kg_k * temperature
        return (total_enthalpy - gas.enthalpy_j_kg(temperature) - kinetic) / scale

    lower = gas.minimum_temperature_k
    if residual(lower) <= 0.0:
        raise DomainError(
            "The static temperature at this Mach number lies below the gas model range.",
            code="temperature_out_of_model_range",
            details={"mach": mach, "total_temperature_k": total_temperature_k},
        )
    return bracketed_root(residual, lower, total_temperature_k, label="static temperature")


def _total_temperature(gas: ThermallyPerfectGas, temperature_k: float, velocity_m_s: float) -> float:
    """Solve ``h(Tt) = h(T) + V^2/2`` for the stagnation temperature."""
    target = gas.enthalpy_j_kg(temperature_k) + 0.5 * velocity_m_s**2
    if velocity_m_s == 0.0:
        return temperature_k
    upper = temperature_k
    while gas.enthalpy_j_kg(upper) < target:
        upper = min(2.0 * upper, gas.maximum_temperature_k)
        if upper >= gas.maximum_temperature_k and gas.enthalpy_j_kg(upper) < target:
            raise DomainError(
                "The stagnation temperature exceeds the gas model range.",
                code="temperature_out_of_model_range",
            )
    scale = gas.cp_j_kg_k(temperature_k) * temperature_k
    return bracketed_root(
        lambda candidate: (gas.enthalpy_j_kg(candidate) - target) / scale,
        temperature_k,
        upper,
        label="stagnation temperature",
    )


def _validate_mach(mach: float, *, supersonic: bool = False) -> None:
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and positive.", field="mach")
    if supersonic and mach <= 1.0:
        raise DomainError("Upstream Mach number must be greater than 1.", field="upstream_mach")


def thermally_perfect_isentropic(
    mach: float, total_temperature_k: float, gas: ThermallyPerfectGas | None = None
) -> ThermallyPerfectIsentropicResult:
    """Isentropic one-dimensional flow of a thermally perfect gas [1]_.

    Parameters
    ----------
    mach : float
        Local Mach number, finite and positive.
    total_temperature_k : float
        Stagnation temperature in K, inside the gas model range.
    gas : ThermallyPerfectGas, optional
        Caloric model; defaults to the NACA 1135 harmonic-vibrator air.

    Returns
    -------
    ThermallyPerfectIsentropicResult
        Static state and ratios to stagnation and sonic conditions.

    Raises
    ------
    DomainError
        If an input is invalid or the static temperature falls outside the
        caloric model range.

    Notes
    -----
    The static temperature solves the adiabatic energy equation
    ``h(Tt) - h(T) = M^2 gamma(T) R T / 2`` (NACA 1135 eq. (186) for the
    harmonic vibrator) [1]_. Along the isentrope
    ``p/pt = exp((phi(T) - phi(Tt))/R)`` (eq. (189)), ``rho/rho_t =
    (p/pt)(Tt/T)`` (thermal equation of state), and the stream-tube area
    follows from continuity, ``A/A* = rho* a*/(rho V)`` (eq. (191)).

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, pp. 627-628.
    """

    gas = HarmonicVibratorGas() if gas is None else gas
    _validate_mach(mach)
    _check_temperature(gas, total_temperature_k)
    r = gas.gas_constant_j_kg_k
    temperature = _static_temperature(gas, mach, total_temperature_k)
    critical_temperature = _static_temperature(gas, 1.0, total_temperature_k)
    phi_total = gas.entropy_function_j_kg_k(total_temperature_k)
    pressure_ratio = exp((gas.entropy_function_j_kg_k(temperature) - phi_total) / r)
    critical_pressure_ratio = exp(
        (gas.entropy_function_j_kg_k(critical_temperature) - phi_total) / r
    )
    density_ratio = pressure_ratio * total_temperature_k / temperature
    critical_density_ratio = critical_pressure_ratio * total_temperature_k / critical_temperature
    sound = _sound_speed(gas, temperature)
    critical_sound = _sound_speed(gas, critical_temperature)
    velocity = mach * sound
    gamma_local = _gamma(gas, temperature)
    return ThermallyPerfectIsentropicResult(
        mach=mach,
        total_temperature_k=total_temperature_k,
        temperature_k=temperature,
        gamma=gamma_local,
        temperature_ratio=temperature / total_temperature_k,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        dynamic_pressure_ratio=0.5 * gamma_local * mach * mach * pressure_ratio,
        area_ratio=(critical_density_ratio * critical_sound) / (density_ratio * velocity),
        velocity_critical_ratio=velocity / critical_sound,
        critical_temperature_k=critical_temperature,
    )


def thermally_perfect_normal_shock(
    upstream_mach: float,
    upstream_temperature_k: float,
    gas: ThermallyPerfectGas | None = None,
) -> ThermallyPerfectNormalShockResult:
    """Normal shock in a thermally perfect gas [1]_.

    Parameters
    ----------
    upstream_mach : float
        Mach number ahead of the shock, finite and greater than 1.
    upstream_temperature_k : float
        Static temperature ahead of the shock, K.
    gas : ThermallyPerfectGas, optional
        Caloric model; defaults to the NACA 1135 harmonic-vibrator air.

    Returns
    -------
    ThermallyPerfectNormalShockResult
        Downstream state and jump ratios.

    Raises
    ------
    DomainError
        If an input is invalid or a temperature leaves the model range.

    Notes
    -----
    Continuity and momentum with ``p = rho R T`` give, for a trial
    downstream temperature, the downstream velocity as the smaller root of
    ``u2^2 - (u1 + R T1/u1) u2 + R T2 = 0`` (the structure of NACA 1135
    eq. (194)) [1]_. The energy equation ``h(T2) + u2^2/2 = h(T1) + u1^2/2``
    (eq. (193)) is then solved for ``T2`` by bisection. The stagnation
    pressure ratio uses ``pt/p = exp((phi(Tt) - phi(T))/R)`` on each side
    (eq. (189)).

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, p. 629.
    """

    gas = HarmonicVibratorGas() if gas is None else gas
    _validate_mach(upstream_mach, supersonic=True)
    _check_temperature(gas, upstream_temperature_k)
    r = gas.gas_constant_j_kg_k
    t1 = upstream_temperature_k
    u1 = upstream_mach * _sound_speed(gas, t1)
    total_temperature = _total_temperature(gas, t1, u1)
    h1 = gas.enthalpy_j_kg(t1)
    momentum_sum = r * t1 / u1 + u1
    upper = min(total_temperature, momentum_sum * momentum_sum / (4.0 * r))

    def downstream_velocity(t2: float) -> float:
        discriminant = max(0.0, momentum_sum * momentum_sum - 4.0 * r * t2)
        return 0.5 * (momentum_sum - sqrt(discriminant))

    def energy(t2: float) -> float:
        u2 = downstream_velocity(t2)
        return (gas.enthalpy_j_kg(t2) + 0.5 * u2 * u2 - h1 - 0.5 * u1 * u1) / (0.5 * u1 * u1)

    t2 = bracketed_root(energy, t1, upper, label="normal-shock temperature")
    u2 = downstream_velocity(t2)
    density_ratio = u1 / u2
    pressure_ratio = density_ratio * t2 / t1
    phi1 = gas.entropy_function_j_kg_k(t1)
    phi2 = gas.entropy_function_j_kg_k(t2)
    phi_t = gas.entropy_function_j_kg_k(total_temperature)
    stagnation_ratio = pressure_ratio * exp((phi1 - phi2) / r)
    total_to_static_2 = exp((phi_t - phi2) / r)
    return ThermallyPerfectNormalShockResult(
        upstream_mach=upstream_mach,
        downstream_mach=u2 / _sound_speed(gas, t2),
        upstream_temperature_k=t1,
        downstream_temperature_k=t2,
        total_temperature_k=total_temperature,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        temperature_ratio=t2 / t1,
        stagnation_pressure_ratio=stagnation_ratio,
        upstream_static_to_downstream_total=1.0 / (pressure_ratio * total_to_static_2),
    )


def _deflection_rad(
    gas: ThermallyPerfectGas, mach: float, temperature_k: float, wave_angle: float
) -> float:
    """Deflection produced by a wave at ``wave_angle``: ``beta - atan(tan(beta)/(rho2/rho1))``."""
    normal_mach = mach * sin(wave_angle)
    if normal_mach <= 1.0 + 1.0e-7:
        return 0.0
    jump = thermally_perfect_normal_shock(normal_mach, temperature_k, gas)
    return wave_angle - atan(tan(wave_angle) / jump.density_ratio)


def thermally_perfect_oblique_shock(
    upstream_mach: float,
    turning_angle_deg: float,
    upstream_temperature_k: float,
    gas: ThermallyPerfectGas | None = None,
) -> ThermallyPerfectObliqueShockResult:
    """Attached weak oblique shock in a thermally perfect gas [1]_.

    Parameters
    ----------
    upstream_mach : float
        Mach number ahead of the wave, greater than 1.
    turning_angle_deg : float
        Flow deflection in degrees, positive.
    upstream_temperature_k : float
        Static temperature ahead of the wave, K.
    gas : ThermallyPerfectGas, optional
        Caloric model; defaults to the NACA 1135 harmonic-vibrator air.

    Returns
    -------
    ThermallyPerfectObliqueShockResult
        Weak-branch wave angle, downstream Mach number and ratios.

    Raises
    ------
    DomainError
        If an input is invalid or the deflection exceeds the attached limit.

    Notes
    -----
    NACA 1135 (pp. 629-630) states that the normal-shock relations apply to
    the oblique wave when ``M1 sin(beta)`` replaces ``M1`` and the same
    upstream static temperature is used [1]_. The tangential velocity is
    conserved, so ``tan(beta - delta) = tan(beta)/(rho2/rho1)`` (continuity
    normal to the wave) and ``M2 = M2n/sin(beta - delta)``. The maximum
    deflection is located by golden-section search and the weak-branch wave
    angle by bisection.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eqs. (195)-(198).
    """

    gas = HarmonicVibratorGas() if gas is None else gas
    _validate_mach(upstream_mach, supersonic=True)
    if not isfinite(turning_angle_deg) or turning_angle_deg <= 0.0:
        raise DomainError("Turning angle must be finite and positive.", field="turning_angle_deg")
    _check_temperature(gas, upstream_temperature_k)
    mach_angle = asin(1.0 / upstream_mach)
    peak_beta, max_deflection = maximize_unimodal(
        lambda beta: _deflection_rad(gas, upstream_mach, upstream_temperature_k, beta),
        mach_angle + 1.0e-9,
        0.5 * pi - 1.0e-9,
        iterations=60,
    )
    target = radians(turning_angle_deg)
    if target > max_deflection:
        raise DomainError(
            f"Detached shock: maximum attached turn is {degrees(max_deflection):.6g} degrees.",
            code="detached_shock",
        )
    wave_angle = bracketed_root(
        lambda beta: _deflection_rad(gas, upstream_mach, upstream_temperature_k, beta) - target,
        mach_angle + 1.0e-9,
        peak_beta,
        x_tolerance=1.0e-13,
        label="oblique-shock wave angle",
    )
    normal = thermally_perfect_normal_shock(
        upstream_mach * sin(wave_angle), upstream_temperature_k, gas
    )
    return ThermallyPerfectObliqueShockResult(
        upstream_mach=upstream_mach,
        downstream_mach=normal.downstream_mach / sin(wave_angle - target),
        turning_angle_deg=turning_angle_deg,
        wave_angle_deg=degrees(wave_angle),
        maximum_turning_angle_deg=degrees(max_deflection),
        normal_shock=normal,
        pressure_ratio=normal.pressure_ratio,
        density_ratio=normal.density_ratio,
        temperature_ratio=normal.temperature_ratio,
        stagnation_pressure_ratio=normal.stagnation_pressure_ratio,
    )
