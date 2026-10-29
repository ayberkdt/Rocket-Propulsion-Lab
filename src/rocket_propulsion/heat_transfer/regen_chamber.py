"""Coupled gas-wall-coolant thermal analysis of a regeneratively cooled thrust chamber.

The solver follows the station-marching scheme of NASA's Rocket Thermal
Evaluation code as documented by Wadel [1]_: the coolant is marched along its
flow direction and, at every wall station, the gas-side film coefficient, the
conduction through the hot wall and the coolant-side film coefficient (with
the channel lands treated as fins) are solved simultaneously for the hot-wall
temperature. Each piece reuses a verified module of this package:

* gas side: Bartz's correlation with the local wall temperature [2]_ — the same
  expression as :func:`rocket_propulsion.thermodynamics.calculate_bartz_heat_transfer`,
  evaluated per station — and an adiabatic-wall (recovery) temperature with
  recovery factor ``Pr^(1/3)`` [3]_;
* optional gas radiation from a weighted-sum-of-gray-gases model, only where
  that model's fit range holds (:mod:`.radiation`);
* hot-wall conduction ``q = k (T_wg - T_wc) / t`` and straight-fin lands
  (:mod:`.conduction`) [4]_;
* coolant side: Gnielinski with the Chen friction factor, or Taylor's hydrogen
  correlation (:mod:`.convection`), with density, heat capacity, viscosity and
  conductivity from the NIST-verified real-fluid equations of state and
  transport correlations (:mod:`rocket_propulsion.thermodynamics.real_fluids`);
* coolant energy ``dh = q dA / mdot`` and Darcy pressure loss along the channel.

The march is explicit in the coolant direction (first order in station
spacing); refinement studies are in ``tests/test_regen_chamber.py``.

References
----------
.. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, "Comparison of High Aspect Ratio Cooling
       Channel Designs for a Rocket Combustion Chamber With Development of an
       Optimized Design", NASA/TM-1998-206313.
.. [2] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101, equations
       4-10, 4-13 and 4-14,
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [3] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
       equations 1, 2 and 4,
       https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
.. [4] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
       Transfer Textbook", v6.00, 2024, Ch. 2 (resistances) and Ch. 4 (fins).
"""

from __future__ import annotations

import itertools
from bisect import bisect_right
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from decimal import Context, Decimal, localcontext
from fractions import Fraction
from math import cos, isfinite, pi, radians, sqrt
from sys import float_info
from typing import Literal

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import (
    coerce_fields,
    integer,
    numeric_range_guard,
    real,
    real_tuple,
)
from rocket_propulsion.geometry import NozzleContour
from rocket_propulsion.heat_transfer import conduction, convection, radiation
from rocket_propulsion.heat_transfer.channel_wall import solve_channel_wall_section
from rocket_propulsion.thermodynamics.materials import get_wall_material
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid
from rocket_propulsion.thermodynamics.real_fluids.transport import transport_tp

CoolantCorrelation = Literal["auto", "gnielinski", "taylor"]
GasRadiation = Literal["none", "wsgg"]
_HYDROGEN = {"parahydrogen", "normal_hydrogen"}


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise DomainError(f"{label} must be finite and positive.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class HotGas:
    """Chamber stagnation state and frozen gas properties for the gas side [1]_.

    Attributes
    ----------
    chamber_pressure_pa, chamber_temperature_k : float
        Stagnation pressure and temperature.
    characteristic_velocity_m_s : float
        ``c*`` of the propellant.
    specific_heat_j_kg_k, viscosity_pa_s, prandtl_number, gamma : float
        Frozen chamber properties used by Bartz's correlation.
    throat_radius_of_curvature_m : float
        Throat wall radius of curvature.
    h2o_mole_fraction, co2_mole_fraction : float
        Frozen chamber mole fractions used by the optional gas radiation.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """

    chamber_pressure_pa: float
    chamber_temperature_k: float
    characteristic_velocity_m_s: float
    specific_heat_j_kg_k: float
    viscosity_pa_s: float
    prandtl_number: float
    gamma: float
    throat_radius_of_curvature_m: float
    h2o_mole_fraction: float = 0.0
    co2_mole_fraction: float = 0.0

    def __post_init__(self) -> None:
        """Own the frozen scalar inputs and reject impossible gas-state domains."""
        coerce_fields(self, 'Hot gas')
        for name in ('chamber_pressure_pa', 'chamber_temperature_k',
                     'characteristic_velocity_m_s', 'specific_heat_j_kg_k',
                     'viscosity_pa_s', 'prandtl_number', 'throat_radius_of_curvature_m'):
            _positive(getattr(self, name), name)
        if self.gamma <= 1.0:
            raise DomainError('Frozen gas gamma must exceed one.', field='gamma')
        for name in ('h2o_mole_fraction', 'co2_mole_fraction'):
            if not 0.0 <= getattr(self, name) <= 1.0:
                raise DomainError('Mole fraction must lie in [0, 1].', field=name)
        if self.h2o_mole_fraction + self.co2_mole_fraction > 1.0:
            raise DomainError('H2O and CO2 mole fractions cannot sum above one.',
                              field='gas_mole_fractions')


@dataclass(frozen=True, slots=True)
class CoolingChannels:
    """Rectangular milled channels around the circumference [1]_.

    Attributes
    ----------
    count : int
        Number of channels.
    width_m, height_m : float
        Channel width (circumferential) and depth (radial).
    hot_wall_thickness_m : float
        Wall between the hot gas and the channel bottom.
    roughness_m : float
        Absolute surface roughness of the channel.
    helix_angle_deg : float, default 0
        Angle between the channels and the chamber axis; 0 is axial.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313 (channel geometry).
    """

    count: int
    width_m: float
    height_m: float
    hot_wall_thickness_m: float
    roughness_m: float = 0.0
    helix_angle_deg: float = 0.0

    def __post_init__(self) -> None:
        coerce_fields(self, 'Cooling channels')
        _channel_scalars(self.count, self.roughness_m, self.helix_angle_deg)
        for name in ('width_m', 'height_m', 'hot_wall_thickness_m'):
            _positive(getattr(self, name), name)


#: Largest accepted helix angle (library choice): beyond it the channel-normal
#: pitch, ``2 pi r cos(beta)/N``, is a small fraction of the circumference.
MAXIMUM_HELIX_ANGLE_DEG = 75.0


@dataclass(frozen=True, slots=True)
class ChannelProfile:
    """Channels whose cross-section varies along the chamber [1]_.

    Real regenerative chambers narrow and deepen the channels toward the
    throat (high-aspect-ratio cooling) to raise the coolant velocity where
    the heat flux peaks. Geometry is tabulated at axial stations and
    interpolated linearly in ``x``; every contour station must lie inside
    the table.

    Attributes
    ----------
    x_m : tuple of float
        Strictly increasing axial stations (same origin as the contour).
    width_m, height_m, hot_wall_thickness_m : tuple of float
        Channel width, depth and hot-wall thickness at each station.
    count : int or tuple of int
        Number of channels, or one count per table row: ``count[i]`` holds
        from ``x_m[i]`` up to the next row, so a change between rows is a
        bifurcation (or merge) at the later row's ``x``.
    roughness_m : float, default 0
        Absolute surface roughness.
    helix_angle_deg : float, default 0
        Angle between the channels and the chamber axis.

    Raises
    ------
    DomainError
        If the table is shorter than two stations, not strictly increasing,
        not finite and positive, or the helix angle is outside [0, 75) deg.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313
           (channel designs varying along the chamber, Tables A-I to A-XV).
    """

    x_m: tuple[float, ...]
    width_m: tuple[float, ...]
    height_m: tuple[float, ...]
    hot_wall_thickness_m: tuple[float, ...]
    count: int | tuple[int, ...]
    roughness_m: float = 0.0
    helix_angle_deg: float = 0.0

    def __post_init__(self) -> None:
        coerce_fields(self, 'Channel profile')
        columns = {name: real_tuple(getattr(self, name), name) for name in
                   ('x_m', 'width_m', 'height_m', 'hot_wall_thickness_m')}
        sizes = {len(column) for column in columns.values()}
        if len(sizes) != 1 or sizes.pop() < 2:
            raise DomainError("Channel profile columns need equal lengths of at least two.", field="x_m")
        for name, values in columns.items():
            if name != "x_m" and min(values) <= 0.0:
                raise DomainError(f"{name} must be positive.", field=name)
            object.__setattr__(self, name, values)
        if any(b <= a for a, b in zip(self.x_m, self.x_m[1:], strict=False)):
            raise DomainError("x_m must be strictly increasing.", field="x_m")
        if isinstance(self.count, tuple | list):
            counts = tuple(integer(value, 'count') for value in self.count)
            if len(counts) != len(self.x_m):
                raise DomainError("count needs one entry per table row.", field="count")
            for value in counts:
                _channel_scalars(value, self.roughness_m, self.helix_angle_deg)
            object.__setattr__(self, "count", counts)
        else:
            object.__setattr__(self, 'count', integer(self.count, 'count'))
            _channel_scalars(self.count, self.roughness_m, self.helix_angle_deg)

    def count_at(self, x_m: float) -> int:
        """Return the channel count at ``x_m`` (piecewise constant, right-continuous)."""
        x_m = real(x_m, 'x_m')
        if isinstance(self.count, int):
            return self.count
        index = min(max(bisect_right(self.x_m, x_m) - 1, 0), len(self.x_m) - 1)
        return self.count[index]

    def geometry_at(self, x_m: float) -> tuple[float, float, float]:
        """Return ``(width, height, hot-wall thickness)`` at ``x_m`` by linear interpolation."""
        x_m = real(x_m, 'x_m')
        span = self.x_m[-1] - self.x_m[0]
        if isfinite(span):
            within = self.x_m[0] - 1.0e-9 * span <= x_m <= self.x_m[-1] + 1.0e-9 * span
        else:
            low, high, coordinate = map(Fraction, (self.x_m[0], self.x_m[-1], x_m))
            tolerance = Fraction(1e-9) * (high-low)
            within = low-tolerance <= coordinate <= high+tolerance
        if not within:
            raise DomainError(f"x = {x_m:.6g} m lies outside the channel profile "
                              f"[{self.x_m[0]:.6g}, {self.x_m[-1]:.6g}] m.", field="x_m")
        index = min(max(bisect_right(self.x_m, x_m), 1), len(self.x_m) - 1)
        x0, x1 = map(Fraction, (self.x_m[index - 1], self.x_m[index]))
        weight = min(Fraction(1), max(Fraction(0), (Fraction(x_m)-x0)/(x1-x0)))
        # Evaluate each complete convex combination before exporting it. Both
        # axial span and an individual weight can be unrepresentable while
        # the positive interpolated geometry is finite. The result lies
        # between the owned, positive finite endpoint dimensions.
        return tuple(float((1-weight)*Fraction(column[index-1])
                           + weight*Fraction(column[index]))  # type: ignore[return-value]
                     for column in (self.width_m, self.height_m, self.hot_wall_thickness_m))


def sudden_area_change_coefficient(area_ratio: float, reynolds: float, friction_factor: float) -> float:
    """Return the loss coefficient ``K1`` of a sudden area change, on the upstream velocity [1]_.

    Parameters
    ----------
    area_ratio : float
        Upstream over downstream flow area ``A1/A2`` (GFSSP's ``(D1/D2)^2``);
        below one is an expansion, above one a reduction.
    reynolds : float
        Upstream Reynolds number.
    friction_factor : float
        Upstream Darcy friction factor.

    Returns
    -------
    float
        ``K1`` such that the total-pressure loss is ``K1 rho v1^2 / 2``.

    Notes
    -----
    GFSSP v6 Eqs. (38)-(41): reduction ``(1.2 + 160/Re)((A1/A2)^2 - 1)`` for
    ``Re <= 2500`` and ``(0.6 + 0.48 f)(A1/A2)((A1/A2) - 1)`` above; expansion
    ``2 (1 - (A1/A2)^2)`` for ``Re <= 4000`` and ``(1 + 0.8 f)(1 - A1/A2)^2``
    above (the Borda-Carnot loss as ``f -> 0``). The manual writes them for
    round pipes with ``(D1/D2)^2``; for channels the flow-area ratio takes
    its place.

    References
    ----------
    .. [1] NASA-GFSSP6-AREA-CHANGE. GFSSP v6 user manual, Sects. 3.1.7.7 and
           3.1.7.8, Eqs. (33), (38)-(41), printed pp. 42-45.
    """
    if not (isfinite(area_ratio) and area_ratio > 0.0 and isfinite(reynolds) and reynolds > 0.0
            and isfinite(friction_factor) and friction_factor >= 0.0):
        raise DomainError("Area ratio, Reynolds number and friction factor must be finite and positive.",
                          field="area_ratio")
    if area_ratio == 1.0:
        return 0.0
    if area_ratio > 1.0:  # reduction
        if reynolds <= 2500.0:
            return (1.2 + 160.0 / reynolds) * (area_ratio**2 - 1.0)
        return (0.6 + 0.48 * friction_factor) * area_ratio * (area_ratio - 1.0)
    if reynolds <= 4000.0:  # expansion
        return 2.0 * (1.0 - area_ratio**2)
    return (1.0 + 0.8 * friction_factor) * (1.0 - area_ratio) ** 2


def _channel_scalars(count: object, roughness: object, helix: object) -> tuple[int, float, float]:
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise DomainError("At least one channel is required.", field="count")
    roughness_value, helix_value = float(roughness), float(helix)  # type: ignore[arg-type]
    if not isfinite(roughness_value) or roughness_value < 0.0:
        raise DomainError("Channel roughness must be finite and non-negative.", field="roughness_m")
    if not (isfinite(helix_value) and 0.0 <= helix_value < MAXIMUM_HELIX_ANGLE_DEG):
        raise DomainError(f"Helix angle must be in [0, {MAXIMUM_HELIX_ANGLE_DEG:g}) degrees.",
                          field="helix_angle_deg")
    return count, roughness_value, cos(radians(helix_value))


@dataclass(frozen=True, slots=True)
class Coolant:
    """Coolant inlet state and flow arrangement.

    Attributes
    ----------
    fluid : str
        Real-fluid key (``"parahydrogen"``, ``"methane"``, ``"oxygen"``, ...).
    inlet_temperature_k, inlet_pressure_pa : float
        Channel inlet state.
    mass_flow_kg_s : float
        Total coolant mass flow over all channels.
    counterflow : bool, default True
        Enter at the nozzle exit and flow toward the injector.
    """

    fluid: str
    inlet_temperature_k: float
    inlet_pressure_pa: float
    mass_flow_kg_s: float
    counterflow: bool = True

    def __post_init__(self) -> None:
        coerce_fields(self, 'Coolant')
        for name in ('inlet_temperature_k', 'inlet_pressure_pa', 'mass_flow_kg_s'):
            _positive(getattr(self, name), name)
        if not isinstance(self.counterflow, bool):
            raise InputError('counterflow must be a boolean.', field='counterflow')


@dataclass(frozen=True, slots=True)
class ThermalStation:
    """Solved temperatures, fluxes and coolant state at one wall station."""

    x_m: float
    radius_m: float
    area_ratio: float
    mach: float
    adiabatic_wall_temperature_k: float
    gas_side_coefficient_w_m2_k: float
    radiative_heat_flux_w_m2: float
    heat_flux_w_m2: float
    hot_wall_temperature_k: float
    cold_wall_temperature_k: float
    coolant_temperature_k: float
    coolant_pressure_pa: float
    coolant_side_coefficient_w_m2_k: float
    coolant_reynolds: float
    land_fin_efficiency: float
    wall_limit_margin_k: float | None
    coolant_velocity_m_s: float = 0.0
    coolant_mach: float = 0.0
    channel_width_m: float = 0.0
    channel_height_m: float = 0.0
    hot_wall_thickness_m: float = 0.0
    hot_wall_mean_temperature_k: float = 0.0
    axial_conduction_flux_w_m2: float = 0.0


#: Coolant Mach number above which the march warns that the flow nears choking (library choice).
COOLANT_MACH_WARNING = 0.5


@dataclass(frozen=True, slots=True)
class RegenChamberResult:
    """Station results and totals of a regenerative cooling solution.

    Attributes
    ----------
    stations : tuple of ThermalStation
        In contour order (injector to exit).
    total_heat_load_w : float
        Heat absorbed by the coolant.
    coolant_outlet_temperature_k, coolant_outlet_pressure_pa : float
        Channel outlet state.
    coolant_pressure_drop_pa : float
        Inlet minus outlet pressure.
    maximum_hot_wall_temperature_k : float
        Largest gas-side wall temperature.
    wall_limit_temperature_k : float or None
        Material limit when a catalog material is used.
    energy_balance_residual_w : float
        ``mdot (h_out - h_in) - Q`` with ``h`` re-evaluated from the EOS at the
        outlet state; measures the enthalpy inversion, not the march.
    warnings : tuple of WarningMessage
        Model statements and range notes.
    wall_model : {"fin", "2d"}
        Wall conduction model used.
    """

    stations: tuple[ThermalStation, ...]
    total_heat_load_w: float
    coolant_outlet_temperature_k: float
    coolant_outlet_pressure_pa: float
    coolant_pressure_drop_pa: float
    maximum_hot_wall_temperature_k: float
    wall_limit_temperature_k: float | None
    energy_balance_residual_w: float
    warnings: tuple[WarningMessage, ...]
    wall_model: str = "fin"


def bartz_coefficient(gas: HotGas, *, throat_radius_m: float, area_ratio: float, mach: float,
                      wall_temperature_k: float) -> tuple[float, float]:
    """Return ``(h_g, T_aw)`` at one station with its own wall temperature [1]_, [2]_.

    With ``C = 1 + (gamma - 1) M^2 / 2`` and ``B = (T_w/T_c C + 1) / 2``,
    ``h = 0.026 cp mu^0.2 Pr^-0.6 D^-0.2 (p_c/c*)^0.8 (D/R_c)^0.1 (A/A*)^-0.9 B^-0.68 C^-0.12``
    and ``T_aw = T_c (1 + Pr^(1/3) (C - 1)) / C`` — the expression used by
    :func:`rocket_propulsion.thermodynamics.calculate_bartz_heat_transfer`.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101, equations
           4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """

    if not isinstance(gas, HotGas):
        raise DomainError('gas must be a HotGas state.', field='gas')
    throat_radius_m = _positive(real(throat_radius_m, 'throat_radius_m'), 'throat_radius_m')
    wall_temperature_k = _positive(real(wall_temperature_k, 'wall_temperature_k'),
                                   'wall_temperature_k')
    area_ratio = real(area_ratio, 'area_ratio')
    mach = real(mach, 'mach')
    if area_ratio < 1.0:
        raise DomainError('Area ratio A/A* must be at least one.', field='area_ratio')
    if mach < 0.0:
        raise DomainError('Mach number must be non-negative.', field='mach')
    try:
        coefficient, recovery = _bartz_native(gas, throat_radius_m, area_ratio, mach,
                                               wall_temperature_k)
        if all(isfinite(v) and v > 0.0 for v in (coefficient, recovery)):
            return coefficient, recovery
    except (OverflowError, ZeroDivisionError):
        pass
    return _bartz_extended(gas, throat_radius_m, area_ratio, mach, wall_temperature_k)


def _bartz_native(gas: HotGas, radius: float, area_ratio: float, mach: float,
                  wall_temperature_k: float) -> tuple[float, float]:
    diameter = 2.0 * radius
    kinetic = 0.5 * (gas.gamma - 1.0) * mach * mach
    compression = 1.0 + kinetic
    pressure_ratio = gas.chamber_pressure_pa / gas.characteristic_velocity_m_s
    curvature_ratio = diameter / gas.throat_radius_of_curvature_m
    temperature_ratio = wall_temperature_k / gas.chamber_temperature_k
    # A nonzero subnormal intermediate can already have large relative error,
    # even if later scaling makes the complete output normal and finite.
    if any(not isfinite(v) or v < float_info.min for v in
           (diameter, pressure_ratio, curvature_ratio, temperature_ratio)):
        raise OverflowError('Bartz ratio needs extended arithmetic')
    sigma_base = 0.5 * (temperature_ratio * compression + 1.0)
    coefficient = 0.026
    for factor in (gas.specific_heat_j_kg_k, gas.viscosity_pa_s**0.2,
                   gas.prandtl_number**-0.6, diameter**-0.2, pressure_ratio**0.8,
                   curvature_ratio**0.1, area_ratio**-0.9, sigma_base**-0.68,
                   compression**-0.12):
        coefficient *= factor
        if not isfinite(coefficient) or coefficient < float_info.min:
            raise OverflowError('Bartz product needs extended arithmetic')
    recovery = gas.chamber_temperature_k * (1.0 + gas.prandtl_number ** (1.0 / 3.0) * kinetic) / compression
    return coefficient, recovery


def _bartz_extended(gas: HotGas, radius: float, area_ratio: float, mach: float,
                    wall_temperature_k: float) -> tuple[float, float]:
    """Evaluate the complete positive expressions before their binary64 export.

    This extends arithmetic range only, not the correlation's physical validity.
    Decimal intermediates retain ratios/products that can overflow or underflow
    in binary64 while the complete film coefficient and recovery remain finite.
    """
    with localcontext(Context(prec=90)):
        d = Decimal.from_float
        diameter = 2 * d(radius)
        kinetic = (d(gas.gamma) - 1) * d(mach)**2 / 2
        compression = 1 + kinetic
        sigma = (1 + d(wall_temperature_k) * compression / d(gas.chamber_temperature_k)) / 2
        coefficient = (Decimal('0.026') * d(gas.specific_heat_j_kg_k)
                       * d(gas.viscosity_pa_s)**Decimal('0.2')
                       * d(gas.prandtl_number)**Decimal('-0.6') * diameter**Decimal('-0.2')
                       * (d(gas.chamber_pressure_pa) / d(gas.characteristic_velocity_m_s))**Decimal('0.8')
                       * (diameter / d(gas.throat_radius_of_curvature_m))**Decimal('0.1')
                       * d(area_ratio)**Decimal('-0.9') * sigma**Decimal('-0.68')
                       * compression**Decimal('-0.12'))
        recovery = d(gas.chamber_temperature_k) * (
            1 + d(gas.prandtl_number)**(Decimal(1) / 3) * kinetic) / compression
        result = float(coefficient), float(recovery)
    if not all(isfinite(v) and v > 0.0 for v in result):
        raise DomainError('Bartz outputs exceed the positive finite float range.',
                          code='numeric_range', field='gas')
    return result


def _axial_network(result: RegenChamberResult, channels: CoolingChannels | ChannelProfile, conductivity: float,
                   closeout: float, counterflow: bool, closeout_conductivity: float | None
                   ) -> tuple[list[int], list[float], list[float]]:
    """Return the stations that carry march area (in axial order), their areas and the conductances between them.

    The last station in coolant order deposits no heat of its own, so it is
    left out of the axial network (an adiabatic end one segment early).
    """
    stations = result.stations
    _, _, cos_helix = _channel_scalars(
        channels.count_at(stations[0].x_m) if isinstance(channels, ChannelProfile) else channels.count,
        channels.roughness_m, channels.helix_angle_deg)
    k_close = conductivity if closeout_conductivity is None else closeout_conductivity

    def section(i: int) -> float:
        """Metal section per unit circumference times k: hot wall, land share, closeout."""
        item = stations[i]
        count = channels.count_at(item.x_m) if isinstance(channels, ChannelProfile) else channels.count
        pitch = 2.0 * pi * (item.radius_m + item.hot_wall_thickness_m) * cos_helix / count
        land_share = max(0.0, 1.0 - item.channel_width_m / pitch)
        return conductivity * (item.hot_wall_thickness_m + land_share * item.channel_height_m) + k_close * closeout

    order = list(reversed(range(len(stations)))) if counterflow else list(range(len(stations)))
    area = {}
    for step, i in enumerate(order[:-1]):
        nxt = stations[order[step + 1]]
        area[i] = 2.0 * pi * stations[i].radius_m * sqrt((nxt.x_m - stations[i].x_m) ** 2
                                                          + (nxt.radius_m - stations[i].radius_m) ** 2)
    nodes = sorted(area)
    conductances = []
    for a_index, b_index in itertools.pairwise(nodes):
        a, b = stations[a_index], stations[b_index]
        length = sqrt((b.x_m - a.x_m) ** 2 + (b.radius_m - a.radius_m) ** 2)
        conductances.append(0.5 * (section(a_index) + section(b_index)) * pi * (a.radius_m + b.radius_m) / length)
    return nodes, [area[i] for i in nodes], conductances


def _mid_wall(item: ThermalStation) -> float:
    return 0.5 * (item.hot_wall_mean_temperature_k + item.cold_wall_temperature_k)


def _implicit_axial_sources(result: RegenChamberResult, nodes: list[int], areas: list[float],
                            conductances: list[float], conductivity: float) -> dict[int, float]:
    """Solve the station network implicitly for the axial fluxes.

    Each station's mid-wall temperature responds to an imposed flux ``S`` as
    ``T = T0 + R S`` with ``R = 1/(1/r_g + 1/r_c)``: ``r_g = 1/h_g + t/(2k)``
    toward the gas and ``r_c = (T - T_b)/(q + S)`` toward the coolant, both
    read from the last march. With ``A S_i = sum G_ij (T_j - T_i)`` this is a
    tridiagonal system in ``T``, solved by the Thomas algorithm.
    """
    stations = result.stations
    n = len(nodes)
    response, base = [], []
    for i in nodes:
        item = stations[i]
        r_gas = 1.0 / item.gas_side_coefficient_w_m2_k + 0.5 * item.hot_wall_thickness_m / conductivity
        to_coolant = item.heat_flux_w_m2 + item.axial_conduction_flux_w_m2
        r_cool = max(_mid_wall(item) - item.coolant_temperature_k, 1.0e-9) / max(to_coolant, 1.0e-9)
        r = 1.0 / (1.0 / r_gas + 1.0 / r_cool)
        response.append(r)
        base.append(_mid_wall(item) - r * item.axial_conduction_flux_w_m2)
    # rows: T_i - (R_i/A_i) sum_j G_ij (T_j - T_i) = base_i
    lower, diagonal, upper, rhs = [0.0] * n, [1.0] * n, [0.0] * n, list(base)
    for k, g in enumerate(conductances):
        for row, other in ((k, k + 1), (k + 1, k)):
            weight = response[row] / areas[row] * g
            diagonal[row] += weight
            if other > row:
                upper[row] -= weight
            else:
                lower[row] -= weight
    for row in range(1, n):
        factor = lower[row] / diagonal[row - 1]
        diagonal[row] -= factor * upper[row - 1]
        rhs[row] -= factor * rhs[row - 1]
    temperature = [0.0] * n
    temperature[-1] = rhs[-1] / diagonal[-1]
    for row in range(n - 2, -1, -1):
        temperature[row] = (rhs[row] - upper[row] * temperature[row + 1]) / diagonal[row]
    net = [0.0] * n
    for k, g in enumerate(conductances):
        flow = g * (temperature[k + 1] - temperature[k])
        net[k] += flow
        net[k + 1] -= flow
    return {i: net[k] / areas[k] for k, i in enumerate(nodes)}


def _with_axial_conduction(**arguments) -> RegenChamberResult:
    """Iterate the march with implicitly solved axial wall conduction until the fluxes settle."""
    channels = arguments["channels"]
    material = arguments["wall_material"]
    conductivity = arguments["wall_thermal_conductivity_w_m_k"]
    if conductivity is None:
        catalog = get_wall_material(material).thermal_conductivity if material is not None else None
        if catalog is None:
            raise InputError("Axial conduction needs a wall conductivity.", field="wall_thermal_conductivity_w_m_k")
        conductivity = catalog.value
    sources: dict[int, float] = {}
    nodes: list[int] = []
    for _ in range(40):
        result = solve_regen_chamber(**arguments, _axial_sources=sources)
        nodes, areas, conductances = _axial_network(result, channels, conductivity, arguments["closeout_thickness_m"],
                                                    arguments["coolant"].counterflow,
                                                    arguments["closeout_conductivity_w_m_k"])
        updated = _implicit_axial_sources(result, nodes, areas, conductances, conductivity)
        scale = max(abs(item.heat_flux_w_m2) for item in result.stations)
        change = max(abs(updated[i] - sources.get(i, 0.0)) for i in updated)
        if change <= 1.0e-6 * scale:
            break
        sources = updated
    else:
        raise ConvergenceError("Axial wall conduction did not settle.", code="axial_conduction_not_converged")
    largest = max(abs(item.axial_conduction_flux_w_m2) for item in result.stations)
    note = WarningMessage(
        "axial_wall_conduction",
        f"Axial wall conduction included; the largest axial flux is {largest:.3g} W/m^2 "
        f"({largest / scale:.2%} of the peak gas flux).", severity="info")
    return replace(result, warnings=(*result.warnings, note))


def _temperature_from_enthalpy(fluid: str, enthalpy: float, pressure: float, guess: float) -> float:
    """Invert ``h(T, p) = enthalpy``: Newton on ``dh/dT = cp``, bisection as fallback."""
    model = reference_fluid(fluid)
    temperature = guess
    for _ in range(12):
        try:
            state = model.state_tp(temperature, pressure)
        except DomainError:
            break
        step = (state.enthalpy_j_kg - enthalpy) / state.cp_j_kg_k
        temperature -= step
        if abs(step) < 1.0e-9 * temperature:
            return temperature

    def residual(temperature: float) -> float:
        return model.state_tp(temperature, pressure).enthalpy_j_kg - enthalpy

    low, high = max(model.model.t_min_k * 1.0001, guess * 0.8), min(model.model.t_max_k, guess * 1.5 + 5.0)
    f_low, f_high = residual(low), residual(high)
    for _ in range(40):
        if f_low <= 0.0 <= f_high:
            break
        if f_low > 0.0:
            low = max(model.model.t_min_k * 1.0001, low * 0.8)
            f_low = residual(low)
        if f_high < 0.0:
            high = min(model.model.t_max_k, high * 1.4)
            f_high = residual(high)
    else:
        raise ConvergenceError("Coolant temperature could not be bracketed from enthalpy.",
                               code="coolant_enthalpy_inversion")
    for _ in range(200):
        middle = 0.5 * (low + high)
        value = residual(middle)
        if value > 0.0:
            high = middle
        else:
            low = middle
        if high - low < 1.0e-9 * middle:
            break
    return 0.5 * (low + high)


@numeric_range_guard("contour")
def solve_regen_chamber(
    *,
    contour: NozzleContour,
    gas: HotGas,
    channels: CoolingChannels | ChannelProfile,
    coolant: Coolant,
    wall_material: str | None = "grcop-84",
    wall_thermal_conductivity_w_m_k: float | None = None,
    coolant_correlation: CoolantCorrelation = "auto",
    gas_radiation: GasRadiation = "none",
    wall_model: str = "fin",
    closeout_thickness_m: float = 0.0,
    closeout_conductivity_w_m_k: float | None = None,
    axial_conduction: bool = False,
    _axial_sources: Mapping[int, float] | None = None,
) -> RegenChamberResult:
    """March the coupled gas-wall-coolant problem along the cooling channels [1]_, [2]_.

    Parameters
    ----------
    contour : NozzleContour
        Wall stations with radius, area ratio and Mach number.
    gas : HotGas
        Chamber state and frozen gas properties.
    channels : CoolingChannels or ChannelProfile
        Constant or axially varying channel geometry; the land width follows
        from the local channel-normal pitch.
    coolant : Coolant
        Coolant fluid, inlet state, flow and direction.
    wall_material : str or None, default "grcop-84"
        Catalog material for conductivity and the temperature limit.
    wall_thermal_conductivity_w_m_k : float, optional
        Overrides the catalog conductivity.
    coolant_correlation : {"auto", "gnielinski", "taylor"}, default "auto"
        Coolant-side correlation. ``"auto"`` uses Taylor's correlation for
        hydrogen, whose wall-to-bulk temperature-ratio term matters for
        supercritical hydrogen (NASA TN D-4332), and Gnielinski otherwise;
        Taylor applies to hydrogen only.
    gas_radiation : {"none", "wsgg"}, default "none"
        Add H2O/CO2 gas radiation where the WSGG fit range holds.
    wall_model : {"fin", "2d"}, default "fin"
        ``"fin"`` treats each land as a straight fin on a hot wall of uniform
        temperature. ``"2d"`` solves steady two-dimensional conduction in the
        half-pitch wall cross-section at every station
        (:func:`~rocket_propulsion.heat_transfer.channel_wall.solve_channel_wall_section`),
        iterating the Bartz and coolant coefficients on the field's mean wall
        temperatures; ``hot_wall_temperature_k`` is then the peak over the
        pitch (above the channel centre) and ``hot_wall_mean_temperature_k``
        the mean. About 10 ms per station.
    closeout_thickness_m : float, default 0
        Closeout over the channels (``"2d"`` only): bonded to the lands, its
        underside cooled by the coolant, its outer face adiabatic.
    closeout_conductivity_w_m_k : float, optional
        Closeout conductivity; defaults to the wall conductivity.
    axial_conduction : bool, default False
        Couple the stations through conduction along the wall: each station
        receives ``sum G_ij (T_j - T_i)`` from its neighbours, with
        ``G = k A' 2 pi r / ds``, ``A'`` the metal section per unit
        circumference (hot wall, the land share of the channel height and
        the closeout) and ``T`` the mid-wall temperature; the march is
        repeated until these fluxes change by less than 1e-6 of the peak gas
        flux. The fluxes sum to zero, so the heat load is redistributed, not
        created.

    Returns
    -------
    RegenChamberResult
        Station temperatures and fluxes, coolant state and totals.

    Raises
    ------
    InputError
        For unknown options or a hydrogen-only correlation on another fluid.
    DomainError
        For non-physical geometry (no land between channels), a coolant
        state outside its equation range, or a two-phase coolant.
    ConvergenceError
        If a station's wall temperature cannot be bracketed.

    Notes
    -----
    At each station the hot-wall temperature solves
    ``q_g(T_wg) = (T_wg - T_b) / (t/k + W / (h_c (w + 2 eta h)))``, where
    ``q_g = h_g(T_wg) (T_aw - T_wg) + q_rad``, ``W = 2 pi (r + t) cos(beta)/N``
    is the channel-normal pitch at the hot wall (``beta`` the helix angle) and
    ``eta`` the efficiency of a land fin of thickness ``W - w`` and height
    ``h`` [3]_; ``w``, ``h`` and ``t`` are the local channel values. The
    coolant then absorbs ``q_g 2 pi r ds`` over the next meridional segment
    ``ds`` and travels ``ds/cos(beta)`` along the channel. As in RTE [1]_
    (Eqs. (14)-(16)), its static pressure falls by the Darcy friction drop
    with the two stations' mean density, velocity ``G/rho`` and hydraulic
    diameter (Chen friction factor) plus the momentum change. The latter is
    Euler's ``dp = -rho v dv`` for steady flow with ``rho v A`` constant [4]_,
    integrated as ``(G_1 + G_2)/2 (v_2 - v_1)``; for a constant channel area
    it equals RTE's ``G^2/rho_2 - G^2/rho_1``, and with a varying area it
    keeps the reversible pressure change of the area change out of the loss.
    Where the channel count changes (bifurcation or merge), the sudden change
    of total flow area adds ``K1 rho_1 v_1^2/2`` with ``K1`` from
    :func:`sudden_area_change_coefficient` [5]_. The downstream state is
    found by fixed-point iteration on the pressure to 1e-12.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313.
    .. [2] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [3] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
           Transfer Textbook", v6.00, 2024, Ch. 4.
    .. [4] NASA-GLENN-CD-NOZZLE. NASA Glenn Research Center, "Nozzle Design -
           Converging/Diverging (CD) Nozzle": ``r V dV = -dp`` with
           ``mdot = r V A``. https://www.grc.nasa.gov/www/k-12/airplane/nozzled.html
    .. [5] NASA-GFSSP6-AREA-CHANGE. GFSSP v6 user manual, Eqs. (33), (38)-(41).
    """

    if coolant_correlation not in ("auto", "gnielinski", "taylor"):
        raise InputError("Coolant correlation must be auto, gnielinski or taylor.", field="coolant_correlation")
    if coolant_correlation == "auto":
        coolant_correlation = "taylor" if coolant.fluid in _HYDROGEN else "gnielinski"
    if coolant_correlation == "taylor" and coolant.fluid not in _HYDROGEN:
        raise InputError("Taylor's correlation applies to hydrogen coolant only.", field="coolant_correlation")
    if wall_model not in ("fin", "2d"):
        raise InputError("Wall model must be fin or 2d.", field="wall_model")
    closeout = float(closeout_thickness_m)
    if not isfinite(closeout) or closeout < 0.0:
        raise DomainError("Closeout thickness must be finite and non-negative.", field="closeout_thickness_m")
    if closeout > 0.0 and wall_model != "2d":
        raise InputError("A closeout needs wall_model='2d'.", field="closeout_thickness_m")
    if not isinstance(axial_conduction, bool):
        raise InputError("axial_conduction must be a boolean.", field="axial_conduction")
    if axial_conduction:
        return _with_axial_conduction(
            contour=contour, gas=gas, channels=channels, coolant=coolant, wall_material=wall_material,
            wall_thermal_conductivity_w_m_k=wall_thermal_conductivity_w_m_k,
            coolant_correlation=coolant_correlation, gas_radiation=gas_radiation, wall_model=wall_model,
            closeout_thickness_m=closeout, closeout_conductivity_w_m_k=closeout_conductivity_w_m_k)
    sources = dict(_axial_sources or {})
    if gas_radiation not in ("none", "wsgg"):
        raise InputError("Gas radiation must be none or wsgg.", field="gas_radiation")
    if isinstance(channels, ChannelProfile):
        geometry_at = channels.geometry_at
        count_at = channels.count_at
        _, roughness, cos_helix = _channel_scalars(count_at(channels.x_m[0]), channels.roughness_m,
                                                   channels.helix_angle_deg)
    else:
        fixed_count, roughness, cos_helix = _channel_scalars(channels.count, channels.roughness_m,
                                                             channels.helix_angle_deg)

        def count_at(x_m: float) -> int:
            return fixed_count

        constant = (_positive(channels.width_m, "width_m"), _positive(channels.height_m, "height_m"),
                    _positive(channels.hot_wall_thickness_m, "hot_wall_thickness_m"))

        def geometry_at(x_m: float) -> tuple[float, float, float]:
            return constant
    mass_flow = _positive(coolant.mass_flow_kg_s, "mass_flow_kg_s")
    limit = None
    if wall_thermal_conductivity_w_m_k is not None:
        conductivity = _positive(wall_thermal_conductivity_w_m_k, "wall_thermal_conductivity_w_m_k")
    elif wall_material is not None:
        catalog = get_wall_material(wall_material).thermal_conductivity
        if catalog is None:
            raise InputError(f"{wall_material} has no verified catalog conductivity; give "
                             "wall_thermal_conductivity_w_m_k.",
                             code="wall_conductivity_unavailable",
                             field="wall_thermal_conductivity_w_m_k")
        conductivity = catalog.value
    else:
        raise InputError("Give a wall material or a wall conductivity.", field="wall_material")
    if wall_material is not None:
        limit = get_wall_material(wall_material).temperature_limit()[0].value

    stations = list(contour.stations)
    order = list(reversed(range(len(stations)))) if coolant.counterflow else list(range(len(stations)))
    fluid = reference_fluid(coolant.fluid)
    geometry = {index: geometry_at(stations[index].x_m) for index in order}
    counts = {index: count_at(stations[index].x_m) for index in order}

    def passage(index: int) -> tuple[float, float]:
        """Return the hydraulic diameter and mass flux at a station."""
        local_width, local_height, _ = geometry[index]
        area = local_width * local_height
        return 2.0 * area / (local_width + local_height), mass_flow / (counts[index] * area)

    temperature, pressure = float(coolant.inlet_temperature_k), float(coolant.inlet_pressure_pa)
    junctions, junction_loss = 0, 0.0
    inlet_state = fluid.state_tp(temperature, pressure)
    enthalpy = inlet_state.enthalpy_j_kg
    total_heat = 0.0
    radiation_skipped = 0
    distance = 0.0
    solved: dict[int, ThermalStation] = {}
    warnings: list[WarningMessage] = []
    # Coolant correlations run with extrapolate=True so the march can finish. Their range
    # warnings are kept only for the converged wall temperature (not for bisection probes)
    # and are reported once per out-of-range parameter at the end.
    extrapolated: dict[str, tuple[str, set[int]]] = {}
    recording = False

    def note(result_warnings: Sequence[WarningMessage], station_index: int) -> None:
        if not recording:
            return
        for item in result_warnings:
            if item.code == "correlation_extrapolated":
                _, _, stated = item.message.partition(" lies outside the stated range of ")
                extrapolated.setdefault(str(item.field), (stated.rstrip("."), set()))[1].add(station_index)

    for step, index in enumerate(order):
        station = stations[index]
        width, height, thickness = geometry[index]
        hydraulic, mass_flux = passage(index)
        circumference_per_channel = 2.0 * pi * (station.radius_m + thickness) * cos_helix / counts[index]
        land = circumference_per_channel - width
        if land <= 0.0:
            raise DomainError("Channels do not fit the circumference at this station.",
                              field="width_m", details={"x_m": station.x_m})
        properties = transport_tp(coolant.fluid, temperature, pressure)
        if properties.phase == "liquid" and fluid.model.p_crit_pa > pressure:
            saturation = fluid.saturation_p(pressure)
            if temperature >= saturation.temperature_k * 0.999:
                raise DomainError("Coolant reaches saturation; two-phase flow is not modelled.",
                                  code="coolant_boiling", field="coolant")
        reynolds = mass_flux * hydraulic / properties.viscosity_pa_s
        friction = convection.chen_friction_factor(reynolds, roughness / hydraulic).value

        def coolant_coefficient(wall_temperature: float, *, temperature=temperature,
                                pressure=pressure, distance=distance, properties=properties,
                                index=index, reynolds=reynolds, friction=friction,
                                mass_flux=mass_flux, hydraulic=hydraulic) -> float:
            if coolant_correlation == "taylor":
                taylor = convection.taylor_hydrogen_coefficient(
                    fluid=coolant.fluid, bulk_temperature_k=temperature, wall_temperature_k=wall_temperature,
                    pressure_pa=pressure, mass_flux_kg_m2_s=mass_flux, diameter_m=hydraulic,
                    distance_from_entrance_m=max(distance, 50.0 * hydraulic),
                    bulk_viscosity_pa_s=properties.viscosity_pa_s,
                    bulk_conductivity_w_m_k=properties.thermal_conductivity_w_m_k,
                    extrapolate=True)
                note(taylor.warnings, index)
                return taylor.heat_transfer_coefficient_w_m2_k
            gnielinski = convection.gnielinski_nusselt(reynolds, properties.prandtl_number,
                                                       darcy_friction_factor=friction, extrapolate=True)
            note(gnielinski.warnings, index)
            return gnielinski.value * properties.thermal_conductivity_w_m_k / hydraulic

        def radiative(wall_temperature: float, *, station=station) -> float:
            nonlocal radiation_skipped
            if gas_radiation == "none" or (gas.h2o_mole_fraction + gas.co2_mole_fraction) <= 0.0:
                return 0.0
            compression = 1.0 + 0.5 * (gas.gamma - 1.0) * station.mach**2
            gas_temperature = gas.chamber_temperature_k / compression
            static_pressure = gas.chamber_pressure_pa * compression ** (-gas.gamma / (gas.gamma - 1.0))
            beam = radiation.mean_beam_length_for_geometry("infinite_cylinder", 2.0 * station.radius_m)
            try:
                emission = radiation.wsgg_mixture_emissivity(
                    temperature_k=gas_temperature, path_length_m=beam, total_pressure_pa=static_pressure,
                    h2o_mole_fraction=gas.h2o_mole_fraction, co2_mole_fraction=gas.co2_mole_fraction,
                    allow_pressure_extrapolation=True)
                absorption = radiation.wsgg_mixture_absorptivity(
                    gas_temperature_k=gas_temperature, wall_temperature_k=wall_temperature, path_length_m=beam,
                    total_pressure_pa=static_pressure, h2o_mole_fraction=gas.h2o_mole_fraction,
                    co2_mole_fraction=gas.co2_mole_fraction, allow_pressure_extrapolation=True)
            except DomainError:
                radiation_skipped += 1
                return 0.0
            return radiation.gas_to_black_wall_heat_flux(
                gas_emissivity=emission.value, gas_absorptivity=absorption.value,
                gas_temperature_k=gas_temperature, wall_temperature_k=wall_temperature)

        station_source = sources.get(index, 0.0)

        def balance(wall_temperature: float, *, station=station, temperature=temperature,
                    land=land, circumference_per_channel=circumference_per_channel,
                    coolant_coefficient=coolant_coefficient, radiative=radiative,
                    width=width, height=height,
                    thickness=thickness,
                    source=station_source) -> tuple[float, float, float, float, float, float]:
            h_gas, recovery = bartz_coefficient(gas, throat_radius_m=contour.throat_radius_m,
                                                area_ratio=station.area_ratio, mach=station.mach,
                                                wall_temperature_k=wall_temperature)
            q_rad = radiative(wall_temperature)
            q_gas = h_gas * (recovery - wall_temperature) + q_rad
            # Taylor's T_w is the coolant-side wall, one conduction drop below the hot wall;
            # bisection probes far below the root are floored at the bulk temperature.
            h_cool = coolant_coefficient(max(wall_temperature - (q_gas + source) * thickness / conductivity,
                                             temperature))
            fin = conduction.straight_fin(
                thermal_conductivity_w_m_k=conductivity, heat_transfer_coefficient_w_m2_k=h_cool,
                perimeter_m=2.0, cross_section_area_m2=land, length_m=height,
                base_temperature_k=temperature + 1.0, fluid_temperature_k=temperature)
            wetted = width + 2.0 * fin.efficiency * height
            resistance = thickness / conductivity + circumference_per_channel / (h_cool * wetted)
            conducted = (wall_temperature - temperature) / resistance
            return q_gas + source - conducted, q_gas, h_gas, recovery, h_cool, fin.efficiency

        recovery_guess = bartz_coefficient(gas, throat_radius_m=contour.throat_radius_m,
                                           area_ratio=station.area_ratio, mach=station.mach,
                                           wall_temperature_k=temperature)[1]
        low, high = temperature, recovery_guess
        f_low, f_high = balance(low)[0], balance(high)[0]
        # heat conducted in axially can lift the wall above the recovery temperature, and heat
        # conducted away can pull it toward the coolant: widen the bracket in those cases
        for _ in range(40):
            if f_high <= 0.0:
                break
            high += 0.25 * (high - temperature) + 1.0
            f_high = balance(high)[0]
        for _ in range(40):
            if f_low >= 0.0 or low <= 0.5 * temperature:
                break
            low = max(0.5 * temperature, low - 0.05 * temperature)
            f_low = balance(low)[0]
        if f_low < 0.0 or f_high > 0.0:
            raise ConvergenceError("Hot-wall temperature is not bracketed.", code="wall_temperature_bracket",
                                   details={"x_m": station.x_m})
        for _ in range(200):
            middle = 0.5 * (low + high)
            if balance(middle)[0] > 0.0:
                low = middle
            else:
                high = middle
            if high - low < 1.0e-9 * middle:
                break
        hot_wall = 0.5 * (low + high)
        local_state = fluid.state_tp(temperature, pressure)
        coolant_velocity = mass_flux / local_state.density_kg_m3
        recording = True
        _, q_gas, h_gas, recovery, h_cool, efficiency = balance(hot_wall)
        source = station_source
        cold_wall = hot_wall - (q_gas + source) * thickness / conductivity
        hot_wall_mean = hot_wall
        if wall_model == "2d":
            # fixed point on the mean wall temperatures; the fin solution is the start
            mean_hot, mean_cold = hot_wall, cold_wall
            for _ in range(60):
                h_gas, recovery = bartz_coefficient(gas, throat_radius_m=contour.throat_radius_m,
                                                    area_ratio=station.area_ratio, mach=station.mach,
                                                    wall_temperature_k=mean_hot)
                q_rad = radiative(mean_hot)
                h_cool = coolant_coefficient(max(mean_cold, temperature))
                section = solve_channel_wall_section(
                    pitch_m=circumference_per_channel, channel_width_m=width, channel_height_m=height,
                    hot_wall_thickness_m=thickness, thermal_conductivity_w_m_k=conductivity,
                    gas_coefficient_w_m2_k=h_gas, recovery_temperature_k=recovery,
                    coolant_coefficient_w_m2_k=h_cool, coolant_temperature_k=temperature,
                    radiative_heat_flux_w_m2=q_rad, imposed_heat_flux_w_m2=source,
                    closeout_thickness_m=closeout, closeout_conductivity_w_m_k=closeout_conductivity_w_m_k)
                settled = (abs(section.mean_hot_wall_temperature_k - mean_hot) <= 1.0e-9 * mean_hot
                           and abs(section.mean_coolant_wall_temperature_k - mean_cold) <= 1.0e-9 * mean_cold)
                mean_hot = section.mean_hot_wall_temperature_k
                mean_cold = section.mean_coolant_wall_temperature_k
                if settled:
                    break
            else:
                raise ConvergenceError("2-D wall temperatures did not settle.", code="wall_2d_not_converged",
                                       details={"x_m": station.x_m})
            q_gas = section.mean_heat_flux_w_m2 - source
            hot_wall, hot_wall_mean, cold_wall = section.peak_hot_wall_temperature_k, mean_hot, mean_cold
        recording = False
        solved[index] = ThermalStation(
            x_m=station.x_m, radius_m=station.radius_m, area_ratio=station.area_ratio, mach=station.mach,
            adiabatic_wall_temperature_k=recovery, gas_side_coefficient_w_m2_k=h_gas,
            radiative_heat_flux_w_m2=q_gas - h_gas * (recovery - hot_wall_mean), heat_flux_w_m2=q_gas,
            hot_wall_temperature_k=hot_wall, cold_wall_temperature_k=cold_wall,
            coolant_temperature_k=temperature, coolant_pressure_pa=pressure,
            coolant_side_coefficient_w_m2_k=h_cool, coolant_reynolds=reynolds,
            land_fin_efficiency=efficiency, wall_limit_margin_k=None if limit is None else limit - hot_wall,
            coolant_velocity_m_s=coolant_velocity, coolant_mach=coolant_velocity / local_state.speed_of_sound_m_s,
            channel_width_m=width, channel_height_m=height, hot_wall_thickness_m=thickness,
            hot_wall_mean_temperature_k=hot_wall_mean, axial_conduction_flux_w_m2=source)
        if step == len(order) - 1:
            break
        following_index = order[step + 1]
        following = stations[following_index]
        segment = sqrt((following.x_m - station.x_m) ** 2 + (following.radius_m - station.radius_m) ** 2)
        path = segment / cos_helix
        heat = (q_gas + source) * 2.0 * pi * station.radius_m * segment
        total_heat += heat
        distance += path
        enthalpy += heat / mass_flow
        density = local_state.density_kg_m3
        next_hydraulic, next_mass_flux = passage(following_index)
        junction = 0.0
        if counts[following_index] != counts[index]:
            # bifurcation or merge: sudden change of the total flow area (GFSSP Eqs. (38)-(41))
            total_area = mass_flow / mass_flux
            next_total_area = mass_flow / next_mass_flux
            coefficient = sudden_area_change_coefficient(total_area / next_total_area, reynolds, friction)
            junction = coefficient * mass_flux**2 / (2.0 * density)
            junctions += 1
            junction_loss += junction
        next_temperature, next_pressure = temperature, pressure
        for _ in range(200):
            # fixed-point iteration on the downstream state: Wadel Eq. (14) and Euler's rho v dv = -dp
            next_temperature = _temperature_from_enthalpy(coolant.fluid, enthalpy, next_pressure, next_temperature)
            next_density = fluid.state_tp(next_temperature, next_pressure).density_kg_m3
            mean_density = 0.5 * (density + next_density)
            velocity, next_velocity = mass_flux / density, next_mass_flux / next_density
            mean_velocity = 0.5 * (velocity + next_velocity)
            viscous = friction * path / (0.5 * (hydraulic + next_hydraulic)) * mean_density * mean_velocity**2 / 2.0
            # integral of rho v dv with the trapezoidal mean of rho v = G
            momentum = 0.5 * (mass_flux + next_mass_flux) * (next_velocity - velocity)
            updated = pressure - viscous - momentum - junction
            if updated <= 0.0:
                raise DomainError("Coolant pressure drop exceeds the inlet pressure.",
                                  code="coolant_pressure_exhausted")
            converged = abs(updated - next_pressure) <= 1.0e-12 * pressure
            next_pressure = updated
            if converged:
                break
        else:
            # the iteration contracts unless the coolant is close to choking
            raise ConvergenceError("Coolant pressure did not converge over a segment; the flow is near choking.",
                                   code="coolant_pressure_not_converged", details={"x_m": following.x_m})
        pressure = next_pressure
        temperature = _temperature_from_enthalpy(coolant.fluid, enthalpy, pressure, next_temperature)

    outlet_enthalpy = fluid.state_tp(temperature, pressure).enthalpy_j_kg
    residual = mass_flow * (outlet_enthalpy - inlet_state.enthalpy_j_kg) - total_heat
    ordered = tuple(solved[index] for index in range(len(stations)))
    hottest = max(item.hot_wall_temperature_k for item in ordered)
    warnings.append(WarningMessage(
        "rte_class_model",
        "One-dimensional gas and coolant march with "
        + ("two-dimensional wall conduction across channel and land" if wall_model == "2d"
           else "a fin model of the land")
        + "; no axial wall conduction, no curvature or entrance enhancement, frozen chamber gas "
        "properties; not validated against hardware data.",
        severity="info"))
    for parameter, (stated, indices) in extrapolated.items():
        warnings.append(WarningMessage(
            "coolant_correlation_extrapolated",
            f"{parameter} lies outside the stated range of {stated} at {len(indices)} of "
            f"{len(stations)} stations; the correlation was extrapolated there.",
            severity="warning", field="coolant"))
    if gas_radiation == "wsgg":
        warnings.append(WarningMessage(
            "wsgg_pressure_extrapolated",
            f"Gas radiation uses 1 atm WSGG fits extrapolated in pressure; {radiation_skipped} station "
            "evaluations were outside the 300-2400 K fit range and got no radiation.",
            severity="critical" if radiation_skipped else "warning"))
    if junctions:
        warnings.append(WarningMessage(
            "channel_count_change",
            f"{junctions} channel-count change(s) add {junction_loss:.4g} Pa of sudden-area-change loss "
            "(GFSSP square expansion/reduction); flow-splitting effects beyond the area change are not modelled.",
            severity="info", field="channels"))
    fastest = max(ordered, key=lambda item: item.coolant_mach)
    if fastest.coolant_mach > COOLANT_MACH_WARNING:
        warnings.append(WarningMessage(
            "coolant_high_mach",
            f"Coolant Mach number reaches {fastest.coolant_mach:.2f} ({fastest.coolant_velocity_m_s:.0f} m/s) at "
            f"x = {fastest.x_m:.4g} m; the flow approaches choking and the momentum pressure drop dominates. "
            "Enlarge the channels or raise the inlet pressure.", severity="critical", field="channels"))
    if limit is not None and hottest > limit:
        warnings.append(WarningMessage(
            "wall_limit_exceeded",
            f"Hot-wall temperature {hottest:.0f} K exceeds the {wall_material} limit {limit:.0f} K.",
            severity="critical"))
    return RegenChamberResult(
        stations=ordered,
        total_heat_load_w=total_heat,
        coolant_outlet_temperature_k=temperature,
        coolant_outlet_pressure_pa=pressure,
        coolant_pressure_drop_pa=coolant.inlet_pressure_pa - pressure,
        maximum_hot_wall_temperature_k=hottest,
        wall_limit_temperature_k=limit,
        energy_balance_residual_w=residual,
        warnings=tuple(warnings),
        wall_model=wall_model,
    )


def hot_gas_from_reactants(*, reactants: Mapping[str, float], chamber_pressure_pa: float,
                           throat_radius_of_curvature_m: float) -> HotGas:
    """Build :class:`HotGas` for any reactant mixture from the CEA-gated solver [1]_, [2]_.

    ``reactants`` maps reactant keys (gas species such as ``"H2"`` at
    298.15 K, or condensed reactants such as ``"O2(L)"``) to moles. Frozen
    ``cp`` and ``gamma`` come from the equilibrium chamber composition;
    viscosity and Prandtl number from NASA CEA's transport data at the
    chamber state.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``. https://github.com/nasa/cea
    """

    chamber_pressure_pa = _positive(real(chamber_pressure_pa, 'chamber_pressure_pa'),
                                     'chamber_pressure_pa')
    throat_radius_of_curvature_m = _positive(real(throat_radius_of_curvature_m,
                                                 'throat_radius_of_curvature_m'),
                                             'throat_radius_of_curvature_m')
    from rocket_propulsion.thermochemistry import (
        EquilibriumProblem,
        LocalEquilibriumProvider,
        solve_rocket_performance,
    )
    from rocket_propulsion.thermochemistry.transport import equilibrium_transport_properties

    mixture = dict(reactants)
    performance = solve_rocket_performance(reactants=mixture, chamber_pressure_pa=chamber_pressure_pa,
                                           area_ratios=(2.0,))
    chamber = LocalEquilibriumProvider().solve(EquilibriumProblem.hp(mixture, pressure_pa=chamber_pressure_pa))
    transport = equilibrium_transport_properties(chamber)
    fractions = dict(performance.chamber.mole_fractions)
    return HotGas(
        chamber_pressure_pa=chamber_pressure_pa,
        chamber_temperature_k=chamber.temperature_k,
        characteristic_velocity_m_s=performance.characteristic_velocity_m_s,
        specific_heat_j_kg_k=chamber.cp_j_kg_k,
        viscosity_pa_s=transport.viscosity_pa_s,
        prandtl_number=transport.frozen_prandtl,
        gamma=chamber.gamma,
        throat_radius_of_curvature_m=throat_radius_of_curvature_m,
        h2o_mole_fraction=fractions.get("H2O", 0.0),
        co2_mole_fraction=fractions.get("CO2", 0.0),
    )


def hot_gas_from_propellant(*, propellant_pair: str, chamber_pressure_pa: float,
                            oxidizer_fuel_ratio: float | None = None,
                            throat_radius_of_curvature_m: float) -> HotGas:
    """Build :class:`HotGas` for a propellant preset; see :func:`hot_gas_from_reactants` [1]_, [2]_.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``. https://github.com/nasa/cea
    """

    from rocket_propulsion.thermochemistry import (
        get_propellant_pair,
        reactants_for_mass_ratio,
    )

    pair = get_propellant_pair(propellant_pair)
    ratio = pair.reference_oxidizer_fuel_ratio if oxidizer_fuel_ratio is None else oxidizer_fuel_ratio
    return hot_gas_from_reactants(reactants=reactants_for_mass_ratio(pair, ratio),
                                  chamber_pressure_pa=chamber_pressure_pa,
                                  throat_radius_of_curvature_m=throat_radius_of_curvature_m)


__all__: Sequence[str] = (
    "Coolant",
    "CoolingChannels",
    "HotGas",
    "RegenChamberResult",
    "ThermalStation",
    "bartz_coefficient",
    "hot_gas_from_propellant",
    "hot_gas_from_reactants",
    "solve_regen_chamber",
)
