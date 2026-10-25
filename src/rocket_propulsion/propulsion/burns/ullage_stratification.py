"""Thermally stratified (n-layer) ullage model with wall and interface convection.

The ullage of a vertical cylindrical tank is represented by gas layers that
keep their identity: pressurant enters at the top as a new layer at the inlet
temperature, older layers are pushed down as liquid drains, and each layer
exchanges heat with the tank-wall strip it covers and, for the bottom layer,
with the liquid surface.  This is the one-dimensional, no-axial-conduction
picture of tank pressurization during outflow of NASA TN D-2797 [1]_, with the
ullage discretized into layers as in NASA's multi-node ullage models [2]_.
Heat-transfer coefficients come from the natural-convection correlation
``Nu = C (Gr Pr)**n`` with the ullage-to-wall and ullage-to-liquid constants
used by NASA MSFC after Ring [3]_, [2]_, and the Grashof number of [4]_.
Pressure is uniform; with ``m_i`` and ``T_i`` the layer masses and
temperatures it follows from ``sum(m_i R T_i) = p V_ullage``.

The same solver with every layer merged after each step is the well-mixed,
single-node ideal-gas ullage assumed by the existing regulated feed model; the
comparison exposes the pressurant-mass penalty of that assumption.  SP-8112
notes that stratification lowers heat transfer to the propellant and so the
required pressurant mass [5]_.

References
----------
.. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H. Roudebush,
       "Parametric Investigation of Liquid Hydrogen Tank Pressurization
       During Outflow", NASA TN D-2797, 1965.
       https://ntrs.nasa.gov/citations/19650012768
.. [2] NASA-20240006787-MULTINODE-ULLAGE. A. LeClair and A. Majumdar, "Nodal
       Modeling of Tank Pressurization and Draining using a Multi-Node-Ullage
       Approach", AIAA SciTech 2025 extended abstract.
       https://ntrs.nasa.gov/citations/20240006787
.. [3] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
       "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
       99-0879 draft. https://ntrs.nasa.gov/citations/19980236009
.. [4] NASA-20190030427-GRASHOF. A. Bandyopadhyay, A. K. Majumdar, A. C.
       LeClair and J. G. Valenzuela, "Multi-Node Modeling of Cryogenic Tank
       Pressurization System using Generalized Fluid System Simulation
       Program". https://ntrs.nasa.gov/citations/20190030427
.. [5] NASA-SP8112-STRATIFICATION. NASA SP-8112, "Pressurization Systems for
       Liquid Rockets", 1975. https://ntrs.nasa.gov/citations/19760015212
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, isfinite, pi

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, real

from .feed_system import IdealPressurantGas

ULLAGE_STRATIFICATION_SCHEMA = "rocket_propulsion_ullage_stratification_v1"
ULLAGE_STRATIFICATION_REFERENCE_IDS = (
    "NASA-TND2797-STRATIFIED-OUTFLOW",
    "NASA-20240006787-MULTINODE-ULLAGE",
    "NASA-19980236009-PRESSURIZATION",
    "NASA-20190030427-GRASHOF",
    "NASA-SP8112-STRATIFICATION",
)
#: Upper bounds that keep one run bounded in time and memory.
MAXIMUM_STEP_COUNT = 200_000
MAXIMUM_SEGMENT_COUNT = 10_000
_ROOT_ITERATIONS = 200


def _positive(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and positive."""
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.", field=label)
    return number


def _nonnegative(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and non-negative."""
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class NaturalConvectionCorrelation:
    """Constants of ``h = K_H C (k / L) (Gr Pr)**n``.

    Parameters
    ----------
    coefficient : float
        Positive constant ``C``.
    exponent : float
        Exponent ``n`` in ``(0, 1)``.
    adjustment_factor : float, default 1.0
        Positive heat-transfer adjustment factor ``K_H``.
    label : str, default "natural-convection"
        Non-blank description.

    Raises
    ------
    DomainError
        If a value is outside its domain.

    Notes
    -----
    Majumdar and Steadman write the ullage heat-transfer coefficient as
    ``h_c = K_H C (k_f / L_s) X**n`` with ``X = Gr Pr`` and take ``K_H = 1``
    [1]_.  ``K_H`` is a calibration factor: a later NASA model set it to 0.5
    to match a measured boil-off rate [2]_.

    References
    ----------
    .. [1] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, p. 4, eqs. 10-13, and p. 5.
           https://ntrs.nasa.gov/citations/19980236009
    .. [2] NASA-20190030427-GRASHOF. A. Bandyopadhyay et al., "Multi-Node
           Modeling of Cryogenic Tank Pressurization System using GFSSP",
           PDF p. 7, below eq. 10. https://ntrs.nasa.gov/citations/20190030427
    """

    coefficient: float
    exponent: float
    adjustment_factor: float = 1.0
    label: str = "natural-convection"

    def __post_init__(self) -> None:
        coerce_fields(self, "Natural-convection correlation")
        _positive(self.coefficient, "Correlation coefficient")
        _positive(self.adjustment_factor, "Heat-transfer adjustment factor")
        if not 0.0 < self.exponent < 1.0:
            raise DomainError("Correlation exponent must be in (0, 1).", field="exponent")


#: Ullage-to-wall natural convection, ``C = 0.54``, ``n = 0.25`` (Ring, as
#: quoted by Majumdar and Steadman, AIAA 99-0879 draft p. 5; LeClair and
#: Majumdar, SciTech 2025 extended abstract eq. 1).
ULLAGE_WALL_CORRELATION = NaturalConvectionCorrelation(
    0.54, 0.25, 1.0, "ullage-to-wall natural convection (Ring 1964 via NASA MSFC)")
#: Ullage-to-liquid natural convection, ``C = 0.27``, ``n = 0.25``, ``K_H = 1``
#: (same sources, p. 5 and eq. 2 respectively).
ULLAGE_LIQUID_CORRELATION = NaturalConvectionCorrelation(
    0.27, 0.25, 1.0, "ullage-to-liquid natural convection (Ring 1964 via NASA MSFC)")


@dataclass(frozen=True, slots=True)
class UllageGasTransport:
    """Constant transport properties of the ullage gas.

    Parameters
    ----------
    thermal_conductivity_w_m_k : float
        Positive thermal conductivity ``k``, W/(m K).
    dynamic_viscosity_pa_s : float
        Positive dynamic viscosity ``mu``, Pa s.

    Raises
    ------
    DomainError
        If a value is not finite and positive.

    Notes
    -----
    Both enter the Grashof and Prandtl numbers of the correlation [1]_.
    They are inputs: this module carries no property data, and their
    temperature dependence is not modelled.

    References
    ----------
    .. [1] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, p. 4, eqs. 12-13.
           https://ntrs.nasa.gov/citations/19980236009
    """

    thermal_conductivity_w_m_k: float
    dynamic_viscosity_pa_s: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Ullage gas transport")
        _positive(self.thermal_conductivity_w_m_k, "Gas thermal conductivity")
        _positive(self.dynamic_viscosity_pa_s, "Gas dynamic viscosity")


def natural_convection_coefficient_w_m2_k(
    correlation: NaturalConvectionCorrelation,
    *,
    gas: IdealPressurantGas,
    transport: UllageGasTransport,
    pressure_pa: float,
    gas_temperature_k: float,
    surface_temperature_k: float,
    length_scale_m: float,
    acceleration_m_s2: float,
) -> float:
    """Return a natural-convection heat-transfer coefficient.

    Parameters
    ----------
    correlation : NaturalConvectionCorrelation
        Constants ``C``, ``n``, ``K_H``.
    gas : IdealPressurantGas
        Ullage gas (``R``, ``cp``).
    transport : UllageGasTransport
        Conductivity and viscosity.
    pressure_pa : float
        Positive ullage pressure, Pa.
    gas_temperature_k : float
        Positive gas temperature at which properties are evaluated, K.
    surface_temperature_k : float
        Positive wall or liquid-surface temperature, K.
    length_scale_m : float
        Positive length scale ``L``; NASA MSFC used the tank diameter, m.
    acceleration_m_s2 : float
        Non-negative acceleration driving buoyancy, m/s^2.

    Returns
    -------
    float
        ``h = K_H C (k / L) (Gr Pr)**n``, W/(m^2 K); zero without a
        temperature difference or acceleration.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    ``Gr = g beta dT L**3 rho**2 / mu**2`` and ``Pr = mu cp / k`` [2]_,
    [1]_, with ideal-gas density ``rho = p / (R T)`` and expansion
    coefficient ``beta = 1 / T`` (both from ``p v = R T`` [3]_).  The
    exponent 3 on ``L`` is printed in [2]_; it is not legible in the scanned
    draft [1]_, and is the one that makes ``Gr`` dimensionless.  Evaluating
    properties at the gas temperature is a choice of this implementation;
    the sources do not state the property temperature.

    References
    ----------
    .. [1] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, p. 4, eqs. 10-13; p. 5 (length scale).
           https://ntrs.nasa.gov/citations/19980236009
    .. [2] NASA-20190030427-GRASHOF. A. Bandyopadhyay et al., "Multi-Node
           Modeling of Cryogenic Tank Pressurization System using GFSSP",
           PDF p. 1 (nomenclature) and p. 6.
           https://ntrs.nasa.gov/citations/20190030427
    .. [3] NASA Glenn Research Center, "Specific Heats".
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    if not isinstance(correlation, NaturalConvectionCorrelation):
        raise DomainError("Correlation must be a NaturalConvectionCorrelation.")
    if not isinstance(gas, IdealPressurantGas) or not isinstance(transport, UllageGasTransport):
        raise DomainError("Gas and transport properties have the wrong type.")
    pressure = _positive(pressure_pa, "Ullage pressure")
    gas_temperature = _positive(gas_temperature_k, "Gas temperature")
    surface = _positive(surface_temperature_k, "Surface temperature")
    length = _positive(length_scale_m, "Convection length scale")
    acceleration = _nonnegative(acceleration_m_s2, "Acceleration")
    return _coefficient(correlation, gas.specific_gas_constant_j_kg_k, gas.specific_heat_cp_j_kg_k,
                        transport, pressure, gas_temperature, surface, length, acceleration)


def _coefficient(correlation: NaturalConvectionCorrelation, gas_constant: float, cp: float,
                 transport: UllageGasTransport, pressure: float, gas_temperature: float,
                 surface: float, length: float, acceleration: float) -> float:
    """Unchecked kernel of `natural_convection_coefficient_w_m2_k`."""
    difference = abs(gas_temperature - surface)
    if difference == 0.0 or acceleration == 0.0:
        return 0.0
    density = pressure / (gas_constant * gas_temperature)
    mu = transport.dynamic_viscosity_pa_s
    k = transport.thermal_conductivity_w_m_k
    grashof = acceleration * difference / gas_temperature * length ** 3 * density * density / (mu * mu)
    prandtl = mu * cp / k
    return (correlation.adjustment_factor * correlation.coefficient * k / length
            * (grashof * prandtl) ** correlation.exponent)


@dataclass(frozen=True, slots=True)
class CylindricalTankWall:
    """Vertical cylindrical tank with a thin, axially non-conducting wall.

    Parameters
    ----------
    inner_radius_m : float
        Positive inner radius ``r``, m.
    cylinder_height_m : float
        Positive internal height, m.
    wall_thickness_m : float
        Positive wall thickness ``delta``, m.
    wall_density_kg_m3 : float
        Positive wall density, kg/m^3.
    wall_specific_heat_j_kg_k : float
        Positive wall specific heat, J/(kg K).

    Raises
    ------
    DomainError
        If a value is not finite and positive.

    Notes
    -----
    Wall mass per strip is ``rho_w (2 pi r dz) delta`` as in [1]_.  Axial
    wall conduction is neglected and the wall wetted by liquid is held at
    the liquid temperature, assumptions (6) and (13) of [2]_.  End domes
    and their heat transfer are not represented [2]_.

    References
    ----------
    .. [1] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, p. 5, eqs. 15-16.
           https://ntrs.nasa.gov/citations/19980236009
    .. [2] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed pp. 2-3.
           https://ntrs.nasa.gov/citations/19650012768
    """

    inner_radius_m: float
    cylinder_height_m: float
    wall_thickness_m: float
    wall_density_kg_m3: float
    wall_specific_heat_j_kg_k: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Tank wall")
        for value, label in (
            (self.inner_radius_m, "Tank radius"),
            (self.cylinder_height_m, "Tank height"),
            (self.wall_thickness_m, "Wall thickness"),
            (self.wall_density_kg_m3, "Wall density"),
            (self.wall_specific_heat_j_kg_k, "Wall specific heat"),
        ):
            _positive(value, label)

    @property
    def cross_section_m2(self) -> float:
        """Internal cross-section ``pi r**2``, m^2 [1]_.

        References
        ----------
        .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
               Roudebush, NASA TN D-2797, 1965, printed p. 8, eq. 23.
               https://ntrs.nasa.gov/citations/19650012768
        """

        return pi * self.inner_radius_m * self.inner_radius_m


@dataclass(frozen=True, slots=True)
class UllageConfiguration:
    """Tank, gas, initial state, outflow and pressurization of one run.

    Parameters
    ----------
    gas : IdealPressurantGas
        Pressurant (calorically perfect ideal gas).
    transport : UllageGasTransport
        Gas conductivity and viscosity.
    tank : CylindricalTankWall
        Tank geometry and wall.
    initial_ullage_height_m : float
        Initial ullage height, in ``(0, cylinder height)``, m.
    initial_pressure_pa : float
        Positive initial ullage pressure, Pa.
    initial_gas_temperature_k : float
        Positive uniform initial gas temperature, K.
    initial_wall_temperature_k : float
        Positive initial temperature of the dry wall, K.
    liquid_temperature_k : float
        Positive liquid (and wetted-wall) temperature, K.
    inlet_temperature_k : float
        Positive pressurant inlet temperature, K.
    liquid_outflow_m3_s : float
        Non-negative constant liquid volume outflow, m^3/s.
    acceleration_m_s2 : float
        Non-negative axial acceleration driving natural convection, m/s^2.
    duration_s : float
        Positive simulated time, s.
    step_count : int
        Number of equal time steps (1 to `MAXIMUM_STEP_COUNT`).
    set_pressure_pa : float or None, default None
        Regulated ullage pressure, Pa.  When given, pressurant is added to
        hold it (never removed).
    prescribed_inflow_kg_s : float, default 0.0
        Non-negative constant pressurant inflow when ``set_pressure_pa`` is
        ``None``, kg/s.
    wall_correlation : NaturalConvectionCorrelation or None
        Ullage-to-wall correlation; default `ULLAGE_WALL_CORRELATION`.
    interface_correlation : NaturalConvectionCorrelation or None
        Ullage-to-liquid correlation; default `ULLAGE_LIQUID_CORRELATION`;
        ``None`` gives an adiabatic interface.
    wall_heat_transfer_coefficient_w_m2_k : float or None, default None
        Constant gas-to-wall coefficient replacing ``wall_correlation``
        (requires ``wall_correlation=None``), W/(m^2 K).
    wall_segment_count : int, default 40
        Number of equal wall strips over the cylinder height.
    initial_layer_count : int, default 4
        Number of equal-mass layers of the initial ullage gas.
    maximum_layer_count : int, default 400
        Cap on gas layers; beyond it the adjacent pair with the smallest
        temperature difference is merged at constant pressure.
    sample_every : int, default 1
        Record a sample every this many steps (the final step is always
        recorded).
    schema : str, default ULLAGE_STRATIFICATION_SCHEMA
        Only the current schema is accepted.

    Raises
    ------
    DomainError
        If a value is outside its domain, both a correlation and a constant
        wall coefficient are given, both pressurization modes are given, or
        the outflow would empty the tank within ``duration_s``.

    Notes
    -----
    Constant inlet temperature, outflow and (when regulated) tank pressure
    follow assumption (10) of [1]_; the default convection constants are
    those of [2]_, [3]_.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed pp. 2-3.
           https://ntrs.nasa.gov/citations/19650012768
    .. [2] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, p. 5. https://ntrs.nasa.gov/citations/19980236009
    .. [3] NASA-20240006787-MULTINODE-ULLAGE. A. LeClair and A. Majumdar,
           AIAA SciTech 2025 extended abstract, p. 4, eqs. 1-2.
           https://ntrs.nasa.gov/citations/20240006787
    """

    gas: IdealPressurantGas
    transport: UllageGasTransport
    tank: CylindricalTankWall
    initial_ullage_height_m: float
    initial_pressure_pa: float
    initial_gas_temperature_k: float
    initial_wall_temperature_k: float
    liquid_temperature_k: float
    inlet_temperature_k: float
    liquid_outflow_m3_s: float
    acceleration_m_s2: float
    duration_s: float
    step_count: int
    set_pressure_pa: float | None = None
    prescribed_inflow_kg_s: float = 0.0
    wall_correlation: NaturalConvectionCorrelation | None = ULLAGE_WALL_CORRELATION
    interface_correlation: NaturalConvectionCorrelation | None = ULLAGE_LIQUID_CORRELATION
    wall_heat_transfer_coefficient_w_m2_k: float | None = None
    wall_segment_count: int = 40
    initial_layer_count: int = 4
    maximum_layer_count: int = 400
    sample_every: int = 1
    schema: str = ULLAGE_STRATIFICATION_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Ullage configuration")
        if self.schema != ULLAGE_STRATIFICATION_SCHEMA:
            raise DomainError("Unsupported ullage schema.", field="schema")
        if not isinstance(self.gas, IdealPressurantGas):
            raise DomainError("Ullage gas must be an IdealPressurantGas.", field="gas")
        if not isinstance(self.transport, UllageGasTransport):
            raise DomainError("Transport must be UllageGasTransport.", field="transport")
        if not isinstance(self.tank, CylindricalTankWall):
            raise DomainError("Tank must be a CylindricalTankWall.", field="tank")
        for value, label in (
            (self.initial_ullage_height_m, "Initial ullage height"),
            (self.initial_pressure_pa, "Initial pressure"),
            (self.initial_gas_temperature_k, "Initial gas temperature"),
            (self.initial_wall_temperature_k, "Initial wall temperature"),
            (self.liquid_temperature_k, "Liquid temperature"),
            (self.inlet_temperature_k, "Inlet temperature"),
            (self.duration_s, "Duration"),
        ):
            _positive(value, label)
        _nonnegative(self.liquid_outflow_m3_s, "Liquid outflow")
        _nonnegative(self.acceleration_m_s2, "Acceleration")
        _nonnegative(self.prescribed_inflow_kg_s, "Prescribed inflow")
        if self.set_pressure_pa is not None:
            _positive(self.set_pressure_pa, "Set pressure")
            if self.prescribed_inflow_kg_s != 0.0:
                raise DomainError("Give either a set pressure or a prescribed inflow.")
        for value, label in (
            (self.wall_correlation, "wall_correlation"),
            (self.interface_correlation, "interface_correlation"),
        ):
            if value is not None and not isinstance(value, NaturalConvectionCorrelation):
                raise DomainError(f"{label} must be a NaturalConvectionCorrelation or None.",
                                  field=label)
        if self.wall_heat_transfer_coefficient_w_m2_k is not None:
            _nonnegative(self.wall_heat_transfer_coefficient_w_m2_k, "Wall coefficient")
            if self.wall_correlation is not None:
                raise DomainError("Use either a wall correlation or a constant wall coefficient.")
        if self.initial_ullage_height_m >= self.tank.cylinder_height_m:
            raise DomainError("Initial ullage height must be below the tank height.")
        for value, label, upper in (
            (self.step_count, "Step count", MAXIMUM_STEP_COUNT),
            (self.wall_segment_count, "Wall segment count", MAXIMUM_SEGMENT_COUNT),
            (self.initial_layer_count, "Initial layer count", MAXIMUM_SEGMENT_COUNT),
            (self.maximum_layer_count, "Maximum layer count", MAXIMUM_SEGMENT_COUNT),
            (self.sample_every, "Sample interval", MAXIMUM_STEP_COUNT),
        ):
            if not 1 <= value <= upper:
                raise DomainError(f"{label} must be between 1 and {upper}.", field=label)
        if self.maximum_layer_count < self.initial_layer_count + 1:
            raise DomainError("Maximum layer count must exceed the initial layer count.")
        drained = self.liquid_outflow_m3_s * self.duration_s
        liquid = self.tank.cross_section_m2 * (
            self.tank.cylinder_height_m - self.initial_ullage_height_m)
        if drained >= liquid:
            raise DomainError("Liquid outflow empties the tank within the duration.")


@dataclass(frozen=True, slots=True)
class UllageSample:
    """Ullage state at one sample time.

    Attributes
    ----------
    time_s : float
        Time, s.
    pressure_pa : float
        Uniform ullage pressure, Pa.
    ullage_volume_m3 : float
        Ullage volume, m^3.
    gas_mass_kg : float
        Pressurant mass in the ullage, kg.
    added_mass_kg : float
        Cumulative pressurant added since ``t = 0``, kg.
    inflow_kg_s : float
        Mean pressurant inflow over the last step, kg/s.
    mean_temperature_k : float
        Mass-mean temperature ``p V / (m R)``, K.
    top_temperature_k : float
        Temperature of the top layer, K.
    bottom_temperature_k : float
        Temperature of the layer at the liquid surface, K.
    heat_to_wall_j : float
        Cumulative heat from gas to wall, J.
    heat_to_liquid_j : float
        Cumulative heat from gas to liquid, J.
    layer_count : int
        Number of gas layers.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed p. 8, eq. 23.
           https://ntrs.nasa.gov/citations/19650012768
    """

    time_s: float
    pressure_pa: float
    ullage_volume_m3: float
    gas_mass_kg: float
    added_mass_kg: float
    inflow_kg_s: float
    mean_temperature_k: float
    top_temperature_k: float
    bottom_temperature_k: float
    heat_to_wall_j: float
    heat_to_liquid_j: float
    layer_count: int


@dataclass(frozen=True, slots=True)
class UllageHistory:
    """Result of one stratified or well-mixed ullage run.

    Attributes
    ----------
    well_mixed : bool
        Whether layers were merged after every step.
    samples : tuple of UllageSample
        Time history.
    final_layer_depths_m : tuple of float
        Depth of each final layer centre below the tank top, m.
    final_layer_temperatures_k : tuple of float
        Final layer temperatures, top to bottom, K.
    wall_segment_depths_m : tuple of float
        Depth of each wall-strip centre, m.
    final_wall_temperatures_k : tuple of float
        Final wall-strip temperatures, K.
    added_mass_kg : float
        Total pressurant added, kg.
    layer_merges : int
        Number of merges forced by ``maximum_layer_count`` (excluding the
        well-mixed merges).
    energy_balance_residual_j : float
        ``dE_gas + dE_wall - H_in + W_liquid + Q_liquid``, J.
    energy_scale_j : float
        Inflow enthalpy plus initial gas internal energy, J, for scaling the
        residual.

    Notes
    -----
    Gas internal energy is ``cv sum(m T)``, inflow enthalpy ``cp T_in dm``
    and boundary work ``p dV`` with the trapezoidal pressure of each step
    [2]_.  The residual is a time-discretization diagnostic: it vanishes to
    round-off for steps held at constant pressure (regulated operation) and
    is first order in the time step when the pressure changes, because the
    new parcel enters at the end-of-step pressure.  The total added mass
    corresponds to eq. (23) of [1]_.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed p. 8, eq. 23.
           https://ntrs.nasa.gov/citations/19650012768
    .. [2] NASA Glenn Research Center, "Specific Heats".
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    well_mixed: bool
    samples: tuple[UllageSample, ...]
    final_layer_depths_m: tuple[float, ...]
    final_layer_temperatures_k: tuple[float, ...]
    wall_segment_depths_m: tuple[float, ...]
    final_wall_temperatures_k: tuple[float, ...]
    added_mass_kg: float
    layer_merges: int
    energy_balance_residual_j: float
    energy_scale_j: float
    reference_ids: tuple[str, ...] = field(default=ULLAGE_STRATIFICATION_REFERENCE_IDS)


def _solve_pressure_ratio(volume: float, pressure: float, gas_constant: float,
                          old_sum: float, kappa: float, new_term: float) -> float:
    """Positive root ``s`` of ``V p s / R - old_sum s**kappa - new_term = 0``.

    Old layers are compressed isentropically (``T ~ s**kappa``) and the new
    layer enters at fixed temperature.  The left side is convex in ``s`` and
    negative at ``s = 0`` (or zero only when ``new_term = 0``, handled in
    closed form), so the positive root is unique and is bracketed by
    doubling.
    """
    scale = volume * pressure / gas_constant
    if new_term == 0.0:
        if old_sum == 0.0:
            raise DomainError("Ullage contains no gas.")
        return (old_sum / scale) ** (1.0 / (1.0 - kappa))

    def residual(s: float) -> float:
        return scale * s - old_sum * s ** kappa - new_term

    high = 1.0
    while residual(high) <= 0.0:
        high *= 2.0
        if high > 1e12:
            raise DomainError("Ullage pressure root not bracketed.")
    low = 0.0
    for _ in range(_ROOT_ITERATIONS):
        middle = 0.5 * (low + high)
        if middle in (low, high):
            break
        if residual(middle) > 0.0:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def _merge_closest(layers: list[list[float]]) -> None:
    """Merge the adjacent layer pair with the smallest temperature difference."""
    index = min(range(len(layers) - 1), key=lambda i: abs(layers[i][1] - layers[i + 1][1]))
    (m1, t1), (m2, t2) = layers[index], layers[index + 1]
    layers[index:index + 2] = [[m1 + m2, (m1 * t1 + m2 * t2) / (m1 + m2)]]


def simulate_ullage(configuration: UllageConfiguration, *, well_mixed: bool = False,
                    ) -> UllageHistory:
    """Integrate the layered (or well-mixed) ullage during outflow.

    Parameters
    ----------
    configuration : UllageConfiguration
        Tank, gas, initial state and pressurization.
    well_mixed : bool, default False
        Merge all layers after every step: the single-node ideal-gas ullage.

    Returns
    -------
    UllageHistory
        Samples, final profiles, added mass and energy diagnostics.

    Raises
    ------
    DomainError
        If the configuration is invalid or a state becomes non-finite or
        non-positive.

    Notes
    -----
    Each step of length ``dt``:

    1. Layer volumes ``m_i R T_i / p`` stack downward from the tank top; each
       layer couples to the dry wall strips it overlaps with conductance
       ``h (2 pi r) dz`` and the bottom layer to the liquid with
       ``h pi r**2``; ``h`` is from the correlation with the tank diameter
       as length scale [2]_, [3]_.
    2. At constant pressure a closed layer obeys ``m cp dT = dQ`` [4]_; with
       sink temperatures frozen over the step it relaxes exponentially
       toward their conductance-weighted mean, which is the heat-transfer
       part of eq. (1) of [1]_ and is unconditionally stable.  The wall
       strips receive exactly the heat the gas loses (eq. 2 of [1]_ in its
       lumped form, eq. 15 of [2]_).
    3. The liquid surface descends by ``Q_L dt / (pi r**2)``.  Old layers
       are compressed isentropically, ``T ~ p**((gamma-1)/gamma)`` [5]_, a
       new layer of mass ``dm`` enters at the inlet temperature, and the new
       pressure satisfies ``sum(m_i R T_i) = p V``.  In regulated mode
       ``dm`` is the mass that restores the set pressure (zero if heating
       alone exceeds it: regulator lockup); otherwise ``dm = mdot dt``.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed p. 3, eqs. 1-3.
           https://ntrs.nasa.gov/citations/19650012768
    .. [2] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank", AIAA
           99-0879 draft, pp. 4-5, eqs. 9-16.
           https://ntrs.nasa.gov/citations/19980236009
    .. [3] NASA-20240006787-MULTINODE-ULLAGE. A. LeClair and A. Majumdar,
           AIAA SciTech 2025 extended abstract, pp. 2-4.
           https://ntrs.nasa.gov/citations/20240006787
    .. [4] NASA Glenn Research Center, "Specific Heats".
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    .. [5] NASA-GRC-ISENTROPIC. NASA Glenn Research Center, "Isentropic
           Flow", eqs. 6-7. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    """

    if not isinstance(configuration, UllageConfiguration):
        raise DomainError("Configuration must be an UllageConfiguration.")
    c = configuration
    gas, tank = c.gas, c.tank
    r_gas = gas.specific_gas_constant_j_kg_k
    cp, cv = gas.specific_heat_cp_j_kg_k, gas.specific_heat_cv_j_kg_k
    kappa = r_gas / cp
    area = tank.cross_section_m2
    perimeter = 2.0 * pi * tank.inner_radius_m
    diameter = 2.0 * tank.inner_radius_m
    dt = c.duration_s / c.step_count
    dz = tank.cylinder_height_m / c.wall_segment_count
    wall_mass = tank.wall_density_kg_m3 * perimeter * dz * tank.wall_thickness_m
    wall_capacity = wall_mass * tank.wall_specific_heat_j_kg_k
    segment_tops = [j * dz for j in range(c.wall_segment_count)]
    walls = [c.initial_wall_temperature_k if top < c.initial_ullage_height_m
             else c.liquid_temperature_k for top in segment_tops]

    pressure = c.initial_pressure_pa
    volume = area * c.initial_ullage_height_m
    initial_mass = pressure * volume / (r_gas * c.initial_gas_temperature_k)
    count = 1 if well_mixed else c.initial_layer_count
    layers = [[initial_mass / count, c.initial_gas_temperature_k] for _ in range(count)]

    def gas_energy() -> float:
        return cv * sum(m * t for m, t in layers)

    energy0 = gas_energy() + wall_capacity * sum(walls)
    added = enthalpy_in = work_out = heat_to_wall = heat_to_liquid = 0.0
    merges = 0
    samples = []

    def sample(time: float, inflow: float) -> UllageSample:
        mass = sum(m for m, _ in layers)
        return UllageSample(
            time_s=time, pressure_pa=pressure, ullage_volume_m3=volume, gas_mass_kg=mass,
            added_mass_kg=added, inflow_kg_s=inflow,
            mean_temperature_k=pressure * volume / (mass * r_gas),
            top_temperature_k=layers[0][1], bottom_temperature_k=layers[-1][1],
            heat_to_wall_j=heat_to_wall, heat_to_liquid_j=heat_to_liquid,
            layer_count=len(layers))

    samples.append(sample(0.0, 0.0))

    def h_wall(t_gas: float, t_wall: float) -> float:
        if c.wall_heat_transfer_coefficient_w_m2_k is not None:
            return c.wall_heat_transfer_coefficient_w_m2_k
        if c.wall_correlation is None:
            return 0.0
        return _coefficient(c.wall_correlation, r_gas, cp, c.transport, pressure, t_gas,
                            t_wall, diameter, c.acceleration_m_s2)

    for step in range(1, c.step_count + 1):
        # 1. Geometry and conductances (two-pointer sweep down the tank).
        depth_bottom = volume / area
        links: list[list[tuple[int, float]]] = [[] for _ in layers]
        top = 0.0
        segment = 0
        for index, (mass, temperature) in enumerate(layers):
            bottom = top + mass * r_gas * temperature / (pressure * area)
            if index == len(layers) - 1:
                bottom = depth_bottom
            while segment < c.wall_segment_count and segment_tops[segment] + dz <= top:
                segment += 1
            j = segment
            while j < c.wall_segment_count and segment_tops[j] < bottom:
                overlap = min(bottom, segment_tops[j] + dz) - max(top, segment_tops[j])
                if overlap > 0.0:
                    conductance = h_wall(temperature, walls[j]) * perimeter * overlap
                    if conductance > 0.0:
                        links[index].append((j, conductance))
                j += 1
            top = bottom
        interface = 0.0
        if c.interface_correlation is not None:
            interface = area * _coefficient(
                c.interface_correlation, r_gas, cp, c.transport, pressure, layers[-1][1],
                c.liquid_temperature_k, diameter, c.acceleration_m_s2)

        # 2. Heat exchange at constant pressure, energy-exact between gas and sinks.
        wall_gain = [0.0] * c.wall_segment_count
        step_liquid = 0.0
        for index, layer in enumerate(layers):
            mass, temperature = layer
            sinks = [(walls[j], g) for j, g in links[index]]
            if index == len(layers) - 1 and interface > 0.0:
                sinks.append((c.liquid_temperature_k, interface))
            total = sum(g for _, g in sinks)
            if total == 0.0:
                continue
            sink_temperature = sum(t * g for t, g in sinks) / total
            relaxed = sink_temperature + (temperature - sink_temperature) * exp(
                -total * dt / (mass * cp))
            heat = mass * cp * (relaxed - temperature)  # J into the gas over dt
            mean_gas = sink_temperature - heat / (total * dt)
            for j, g in links[index]:
                wall_gain[j] -= g * (walls[j] - mean_gas) * dt
            if index == len(layers) - 1 and interface > 0.0:
                step_liquid -= interface * (c.liquid_temperature_k - mean_gas) * dt
            layer[1] = relaxed
        for j, gain in enumerate(wall_gain):
            if gain:
                walls[j] += gain / wall_capacity
                heat_to_wall += gain
        heat_to_liquid += step_liquid

        # 3. Outflow, inflow and the new uniform pressure.
        new_volume = volume + c.liquid_outflow_m3_s * dt
        old_sum = sum(m * t for m, t in layers)
        if c.set_pressure_pa is not None:
            ratio = c.set_pressure_pa / pressure
            new_mass = (new_volume * c.set_pressure_pa / r_gas
                        - old_sum * ratio ** kappa) / c.inlet_temperature_k
            if new_mass < 0.0:
                new_mass = 0.0
                ratio = _solve_pressure_ratio(new_volume, pressure, r_gas, old_sum, kappa, 0.0)
        else:
            new_mass = c.prescribed_inflow_kg_s * dt
            ratio = _solve_pressure_ratio(new_volume, pressure, r_gas, old_sum, kappa,
                                          new_mass * c.inlet_temperature_k)
        new_pressure = pressure * ratio
        if not isfinite(new_pressure) or new_pressure <= 0.0:
            raise DomainError("Ullage pressure became non-finite or non-positive.")
        compression = ratio ** kappa
        for layer in layers:
            layer[1] *= compression
        if new_mass > 0.0:
            layers.insert(0, [new_mass, c.inlet_temperature_k])
        work_out += 0.5 * (pressure + new_pressure) * (new_volume - volume)
        enthalpy_in += new_mass * cp * c.inlet_temperature_k
        added += new_mass
        pressure, volume = new_pressure, new_volume
        if well_mixed and len(layers) > 1:
            mass = sum(m for m, _ in layers)
            layers[:] = [[mass, sum(m * t for m, t in layers) / mass]]
        while len(layers) > c.maximum_layer_count:
            _merge_closest(layers)
            merges += 1
        if any(not isfinite(t) or t <= 0.0 for _, t in layers):
            raise DomainError("Ullage temperature became non-finite or non-positive.")
        if step % c.sample_every == 0 or step == c.step_count:
            samples.append(sample(step * dt, new_mass / dt))

    depths = []
    top = 0.0
    for mass, temperature in layers:
        height = mass * r_gas * temperature / (pressure * area)
        depths.append(top + 0.5 * height)
        top += height
    energy1 = gas_energy() + wall_capacity * sum(walls)
    residual = energy1 - energy0 - enthalpy_in + work_out + heat_to_liquid
    scale = enthalpy_in + cv * initial_mass * c.initial_gas_temperature_k
    return UllageHistory(
        well_mixed=well_mixed,
        samples=tuple(samples),
        final_layer_depths_m=tuple(depths),
        final_layer_temperatures_k=tuple(t for _, t in layers),
        wall_segment_depths_m=tuple(top_j + 0.5 * dz for top_j in segment_tops),
        final_wall_temperatures_k=tuple(walls),
        added_mass_kg=added,
        layer_merges=merges,
        energy_balance_residual_j=residual,
        energy_scale_j=scale,
    )


@dataclass(frozen=True, slots=True)
class UllageModelComparison:
    """Stratified versus well-mixed pressurant requirement at constant pressure.

    Attributes
    ----------
    stratified : UllageHistory
        Layered run.
    well_mixed : UllageHistory
        Single-node (well-mixed) run with the same inputs.
    ideal_added_mass_kg : float
        ``p dV / (R T_in)``: pressurant needed with no heat transfer,
        eq. (24) of [1]_, kg.
    stratified_collapse_factor : float
        Stratified added mass divided by the ideal mass, eq. (25) of [1]_.
    well_mixed_collapse_factor : float
        Well-mixed added mass divided by the ideal mass.
    sp8112_ullage_mass_kg : float
        First-cut total ullage pressurant mass of SP-8112 eq. (14),
        ``V (rho_f)_u``, with the final ullage volume and the density at the
        set pressure and the lower of liquid and initial wall temperature
        [2]_, kg.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed p. 8, eqs. 24-25.
           https://ntrs.nasa.gov/citations/19650012768
    .. [2] NASA-SP8112-STRATIFICATION. NASA SP-8112, "Pressurization Systems
           for Liquid Rockets", 1975, printed p. 53, eq. 14.
           https://ntrs.nasa.gov/citations/19760015212
    """

    stratified: UllageHistory
    well_mixed: UllageHistory
    ideal_added_mass_kg: float
    stratified_collapse_factor: float
    well_mixed_collapse_factor: float
    sp8112_ullage_mass_kg: float


def compare_ullage_models(configuration: UllageConfiguration) -> UllageModelComparison:
    """Run the layered and well-mixed models and form collapse factors.

    Parameters
    ----------
    configuration : UllageConfiguration
        Regulated configuration whose initial pressure equals the set
        pressure and whose outflow is positive.

    Returns
    -------
    UllageModelComparison
        Both histories, the ideal (adiabatic) added mass and the ratios.

    Raises
    ------
    DomainError
        If the configuration is not regulated at its initial pressure or has
        no outflow.

    Notes
    -----
    The ideal mass ``m_i = p (V_f - V_0) / (R T_g)`` and the mass ratio
    ``m / m_i`` (the collapse factor) are eqs. (24)-(25) of [1]_.  The
    well-mixed run reproduces the single-node ullage of the existing
    regulated feed model; SP-8112 states that stratification reduces heat
    transfer to the propellant and the pressurant mass required [2]_.

    References
    ----------
    .. [1] NASA-TND2797-STRATIFIED-OUTFLOW. D. A. Mandell and W. H.
           Roudebush, NASA TN D-2797, 1965, printed p. 8.
           https://ntrs.nasa.gov/citations/19650012768
    .. [2] NASA-SP8112-STRATIFICATION. NASA SP-8112, "Pressurization Systems
           for Liquid Rockets", 1975, printed p. 90, section 2.3.1.2, and
           p. 53, eq. 14. https://ntrs.nasa.gov/citations/19760015212
    """

    if not isinstance(configuration, UllageConfiguration):
        raise DomainError("Configuration must be an UllageConfiguration.")
    c = configuration
    if c.set_pressure_pa is None or c.set_pressure_pa != c.initial_pressure_pa:
        raise DomainError("Comparison requires regulation at the initial pressure.")
    if c.liquid_outflow_m3_s == 0.0:
        raise DomainError("Comparison requires a positive liquid outflow.")
    stratified = simulate_ullage(c)
    mixed = simulate_ullage(c, well_mixed=True)
    r_gas = c.gas.specific_gas_constant_j_kg_k
    displaced = c.liquid_outflow_m3_s * c.duration_s
    ideal = c.set_pressure_pa * displaced / (r_gas * c.inlet_temperature_k)
    final_volume = c.tank.cross_section_m2 * c.initial_ullage_height_m + displaced
    coldest = min(c.liquid_temperature_k, c.initial_wall_temperature_k)
    return UllageModelComparison(
        stratified=stratified,
        well_mixed=mixed,
        ideal_added_mass_kg=ideal,
        stratified_collapse_factor=stratified.added_mass_kg / ideal,
        well_mixed_collapse_factor=mixed.added_mass_kg / ideal,
        sp8112_ullage_mass_kg=final_volume * c.set_pressure_pa / (r_gas * coldest),
    )
