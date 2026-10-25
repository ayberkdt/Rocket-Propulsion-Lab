"""Prescribed-state Bartz gas-side heat-transfer profiles.

The gas-side coefficient is the Bartz correlation as printed in NASA
SP-125 [1]_, evaluated with user-supplied constant stagnation properties,
prescribed station Mach numbers and area ratios, and a uniform prescribed
wall temperature. The recovery-temperature convention and the signed heat
flux follow Linne et al. [2]_. There is no coupled wall-conduction or coolant
solve, no radiation, reacting boundary-layer, film-cooling or deposition
model, and no hardware validation; results carry
``physical_model_qualified=False``.

Inputs and results are owned, immutable snapshots: every result replays its
station values, integrals, peak, warnings and policy fields from its own
inputs on construction. Replay is an internal-consistency check, not
authentication of measurements or hardware.

References
----------
.. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
       Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
       equations 1, 2 and 4,
       https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
from functools import lru_cache
from math import inf, isfinite, nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import _ordered_items, real, text
from rocket_propulsion.geometry import NozzleContour, NozzleStation

from .bartz_numerics import _BARTZ_POLICY, _profile_exports

_INTEGRATION='piecewise_linear_heat_flux_and_radius_on_straight_meridional_segments'
_RESIDUAL='integral_ledger_identity_not_measured_energy_closure'
_CORRELATION='Bartz prescribed constant-property gas-side correlation'
_PROPERTIES='user-supplied constant stagnation bulk gas properties'
_SOURCE='NASA SP-125 equations 4-13/4-14; Linne et al. 2012 recovery and signed flux equations 1/4'
_URL='https://ntrs.nasa.gov/citations/19710019929'
_VALIDITY=('prescribed gas state and uniform wall temperature; no conjugate wall/coolant solve',
           'empirical operating-domain applicability is unqualified',
           'no radiation, reacting boundary-layer, film cooling or deposition model',
           'linear sampled-profile reconstruction; no continuous-profile quadrature error bound')


def _owned(value,kind,label):
    """Return a fresh ``kind`` instance copied from a mapping or a ``kind`` record.

    Copying through the constructor re-runs its validation, so callers own
    an independent, checked snapshot. Malformed mappings raise
    :class:`DomainError` naming ``label``.
    """
    if isinstance(value,Mapping):
        try:return kind(**dict(value))
        except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
            raise DomainError('Malformed retained mapping.',field=label) from error
    if isinstance(value,kind):
        return kind(**{f.name:getattr(value,f.name) for f in fields(kind)})
    raise DomainError('Typed retained record required.',field=label)


@dataclass(frozen=True,slots=True)
class BartzContourStation:
    """One prescribed station of the thermal contour.

    Attributes
    ----------
    x_m : float
        Axial position, m.
    radius_m : float
        Local wall radius, m; must be positive.
    area_ratio : float
        Local area ratio ``A/At``; must be positive.
    mach : float
        Prescribed local Mach number; must be nonnegative.

    Raises
    ------
    DomainError
        If a value is not a finite real or violates the sign conditions.

    Notes
    -----
    Mach number and area ratio are prescribed by the caller; this record
    does not check isentropic consistency between them.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    x_m:float
    radius_m:float
    area_ratio:float
    mach:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.radius_m<=0 or self.area_ratio<=0 or self.mach<0:
            raise DomainError('Radius/area ratio must be positive; Mach nonnegative.',field='Thermal contour')


def _rounding_cell(value):
    """Return the exact endpoints of the binary64 rounding cell around ``value``.

    The cell is the set of real numbers that round to ``value`` under
    round-to-nearest. It describes floating-point representation only, not
    measurement uncertainty.
    """
    center=Fraction(value);before=Fraction(nextafter(value,-inf));after=nextafter(value,inf)
    return (center+before)/2, (center+Fraction(after))/2 if isfinite(after) else center+(center-before)/2


@dataclass(frozen=True,slots=True)
class BartzContour:
    """Axially ordered stations sampled on a nozzle wall, with the throat radius.

    Attributes
    ----------
    throat_radius_m : float
        Throat radius, m; the Bartz diameter is twice this value.
    stations : tuple of BartzContourStation
        At least two stations with strictly increasing ``x_m``.

    Raises
    ------
    DomainError
        If the throat radius is not positive, fewer than two stations are
        given, ``x_m`` is not strictly increasing, or a station's area ratio
        disagrees with ``(radius / throat_radius)^2`` beyond the binary64
        rounding cells of the three values.

    Notes
    -----
    The area-ratio check is a representation-coherence guard, not a
    geometry calibration; values are never clamped.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    throat_radius_m:float
    stations:tuple[BartzContourStation,...]

    def __post_init__(self):
        radius=real(self.throat_radius_m,'Throat radius')
        if radius<=0:raise DomainError('Throat radius must be positive.',field='Throat radius')
        object.__setattr__(self,'throat_radius_m',radius)
        points=tuple(_owned(s,BartzContourStation,'Contour station') for s in _ordered_items(self.stations,'Contour stations'))
        if len(points)<2:raise DomainError('A surface profile needs at least two stations.',field='Contour stations')
        rlo,rhi=_rounding_cell(radius)
        for index,point in enumerate(points):
            if index and point.x_m<=points[index-1].x_m:
                raise DomainError('Contour x must increase strictly.',field='Contour stations')
            lo,hi=_rounding_cell(point.radius_m);alo,ahi=_rounding_cell(point.area_ratio)
            if hi<rlo or ahi<(lo/rhi)**2 or alo>(hi/rlo)**2:
                raise DomainError('Area ratio and radius disagree beyond binary64 representation cells.',field='Contour geometry')
        object.__setattr__(self,'stations',points)


def _contour(value):
    """Convert a :class:`NozzleContour`, mapping or record into a :class:`BartzContour`.

    Only the throat radius and each station's ``x_m``, ``radius_m``,
    ``area_ratio`` and ``mach`` are retained; other contour metadata is not
    used.
    """
    if isinstance(value,NozzleContour):
        points=[]
        for s in _ordered_items(value.stations,'Contour stations'):
            if not isinstance(s,NozzleStation):raise DomainError('Nozzle station required.',field='Contour stations')
            points.append(BartzContourStation(s.x_m,s.radius_m,s.area_ratio,s.mach))
        return BartzContour(value.throat_radius_m,tuple(points))
    return _owned(value,BartzContour,'Thermal contour')


@dataclass(frozen=True,slots=True)
class BartzInputs:
    """Owned, validated inputs of a prescribed Bartz profile.

    Attributes
    ----------
    contour : BartzContour
        Thermal contour; a :class:`~rocket_propulsion.geometry.NozzleContour`
        or mapping is converted on construction.
    chamber_pressure_pa : float
        Chamber stagnation pressure ``Pc``, Pa.
    chamber_temperature_k : float
        Chamber stagnation temperature ``Tc``, K.
    wall_temperature_k : float
        Uniform prescribed gas-side wall temperature ``Tw``, K.
    characteristic_velocity_m_s : float
        Characteristic velocity ``c*``, m/s.
    dynamic_viscosity_pa_s : float
        Constant gas dynamic viscosity ``mu``, Pa s.
    specific_heat_j_kg_k : float
        Constant gas specific heat ``cp``, J/(kg K).
    prandtl_number : float
        Constant Prandtl number ``Pr``.
    gamma : float
        Ratio of specific heats; must exceed one.
    throat_radius_of_curvature_m : float
        Throat wall radius of curvature ``Rc``, m.

    Raises
    ------
    DomainError
        If a scalar is not a finite positive real, ``gamma <= 1``, or
        ``Tw >= Tc``.

    Notes
    -----
    ``Tw < Tc`` is a restriction of this prescribed model, not a physical
    law; it does not imply ``Tw`` below the local adiabatic wall
    temperature, so negative local heat flux remains possible.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    contour:BartzContour
    chamber_pressure_pa:float
    chamber_temperature_k:float
    wall_temperature_k:float
    characteristic_velocity_m_s:float
    dynamic_viscosity_pa_s:float
    specific_heat_j_kg_k:float
    prandtl_number:float
    gamma:float
    throat_radius_of_curvature_m:float

    def __post_init__(self):
        object.__setattr__(self,'contour',_contour(self.contour))
        for f in fields(self):
            if f.name=='contour':continue
            value=real(getattr(self,f.name),f.name)
            if value<=0:raise DomainError('Thermal quantities must be positive.',field=f.name)
            object.__setattr__(self,f.name,value)
        if self.gamma<=1:raise DomainError('Gamma must exceed one.',field='gamma')
        if self.wall_temperature_k>=self.chamber_temperature_k:
            raise DomainError('This prescribed model requires wall temperature below chamber temperature.',field='wall_temperature_k')


@dataclass(frozen=True,slots=True)
class HeatFluxStation:
    """Finite thermal results at one contour station.

    Attributes
    ----------
    x_m, radius_m, mach : float
        Station position, radius and prescribed Mach number.
    adiabatic_wall_temperature_k : float
        Recovery (adiabatic wall) temperature ``Taw``, K.
    gas_side_coefficient_w_m2_k : float
        Bartz gas-side heat-transfer coefficient ``h``, W/(m^2 K).
    heat_flux_w_m2 : float
        Signed heat flux ``h (Taw - Tw)``, positive into the wall, W/m^2.
    cumulative_heat_w : float
        Signed heat integrated over the wall from the first station, W.

    Raises
    ------
    DomainError
        If a value is not a finite real, the radius, ``Taw`` or ``h`` is not
        positive, or the Mach number is negative.

    Notes
    -----
    This record checks finiteness and signs only. Replay of the equations
    and of the cumulative integral belongs to the enclosing
    :class:`BartzHeatTransferResult`.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """
    x_m:float
    radius_m:float
    mach:float
    adiabatic_wall_temperature_k:float
    gas_side_coefficient_w_m2_k:float
    heat_flux_w_m2:float
    cumulative_heat_w:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.radius_m<=0 or self.mach<0 or self.adiabatic_wall_temperature_k<=0 or self.gas_side_coefficient_w_m2_k<=0:
            raise DomainError('Invalid finite thermal station.',field='Heat-flux station')


def _warnings(value):
    """Normalize warning records; at most five are accepted.

    The bound covers two base warnings, two screening-threshold warnings and
    one signed-flux warning.
    """
    from collections.abc import Set as AbstractSet
    from itertools import islice
    if isinstance(value,(str,bytes,bytearray,memoryview,Mapping,AbstractSet)):
        raise DomainError('Warnings must be ordered records.',field='Thermal warnings')
    try:items=tuple(islice(value,6))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Malformed thermal warnings.',field='Thermal warnings') from error
    if len(items)>5:raise DomainError('Too many thermal warnings.',field='Thermal warnings')
    result=[]
    for item in items:
        item=_owned(item,WarningMessage,'Thermal warning')
        severity=text(item.severity,'Warning severity')
        if severity not in ('info','warning','critical'):raise DomainError('Unknown warning severity.',field='Warning severity')
        result.append(WarningMessage(text(item.code,'Warning code'),text(item.message,'Warning message'),
            severity,None if item.field is None else text(item.field,'Warning field')))
    return tuple(result)


@dataclass(frozen=True,slots=True)
class BartzHeatTransferResult:
    """Replayed prescribed Bartz profile with integrals, peak and qualification flags.

    Attributes
    ----------
    stations : tuple of HeatFluxStation
        One record per contour station.
    peak_heat_flux_w_m2 : float
        Largest signed station heat flux, W/m^2.
    peak_location_m : float
        Axial position of the first station that attains the peak, m.
    total_heat_load_w : float
        Signed heat integrated over the sampled wall, W.
    energy_balance_residual_w : float
        Always zero; labelled by ``residual_kind`` as an integral ledger
        identity, not a measured energy closure.
    correlation, property_model, source, source_url : str
        Model and source labels.
    validity : tuple of str
        Declared scope limits of the prescribed model.
    warnings : tuple of WarningMessage
        Model, ledger and screening-threshold warnings.
    inputs : BartzInputs
        Owned inputs from which the result is replayed.
    contains_negative_heat_flux : bool
        True when some station transfers heat toward the gas.
    physical_model_qualified : bool
        Always False for this model.
    quadrature_error_bound_w : float or None
        ``None``: no error bound against the unsampled continuous profile is
        claimed.
    integration_model, residual_kind, numerical_policy, pi_basis : str
        Identifiers of the integration rule, residual meaning, numerical
        evaluation policy and the binary64 value of pi used for geometry.

    Raises
    ------
    DomainError
        If a field is malformed, a flag is not an actual ``bool``, or any
        value differs from the replay of ``inputs``.

    Notes
    -----
    With ``D = 2 throat_radius``, ``Ar = A/At``, ``C = 1 + (gamma - 1) M^2 / 2``
    and ``B = (Tw/Tc C + 1)/2``, the coefficient is [1]_

        h = 0.026 mu^(1/5) cp Pr^(-3/5) D^(-1/5) (Pc/c*)^(4/5) (D/Rc)^(1/10)
            Ar^(-9/10) B^(-17/25) C^(-3/25)

    The turbulent recovery convention [2]_ uses ``r = Pr^(1/3)`` and
    ``Taw = Tc (1 + r (C - 1)) / C``; the signed flux is
    ``q = h (Taw - Tw)``. On each straight meridional segment of length
    ``L`` both ``q`` and the radius are reconstructed linearly, and the exact
    integral of that reconstruction is
    ``pi L / 3 [(2 r0 + r1) q0 + (r0 + 2 r1) q1]``. This bounds nothing about
    the unsampled true profile. The printed decimal coefficient and exponents
    are evaluated as exact rationals; that is a numerical convention, not a
    claim that they are exact physical constants.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """
    stations:tuple[HeatFluxStation,...]
    peak_heat_flux_w_m2:float
    peak_location_m:float
    total_heat_load_w:float
    energy_balance_residual_w:float
    correlation:str
    property_model:str
    validity:tuple[str,...]
    warnings:tuple[WarningMessage,...]
    source:str
    source_url:str
    inputs:BartzInputs|None=None
    contains_negative_heat_flux:bool=False
    physical_model_qualified:bool=False
    quadrature_error_bound_w:float|None=None
    integration_model:str=_INTEGRATION
    residual_kind:str=_RESIDUAL
    numerical_policy:str=_BARTZ_POLICY
    pi_basis:str='fixed_binary64_math_pi'

    def __post_init__(self):
        source=_owned(self.inputs,BartzInputs,'Bartz inputs');object.__setattr__(self,'inputs',source)
        for name in ('peak_heat_flux_w_m2','peak_location_m','total_heat_load_w','energy_balance_residual_w'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        for name in ('contains_negative_heat_flux','physical_model_qualified'):
            if type(getattr(self,name)) is not bool:raise DomainError('Thermal flags must be actual booleans.',field=name)
        if self.quadrature_error_bound_w is not None:
            object.__setattr__(self,'quadrature_error_bound_w',real(self.quadrature_error_bound_w,'Quadrature error bound'))
        for name in ('correlation','property_model','source','source_url','integration_model','residual_kind','numerical_policy','pi_basis'):
            object.__setattr__(self,name,text(getattr(self,name),name))
        object.__setattr__(self,'validity',tuple(text(s,'Thermal scope') for s in _ordered_items(self.validity,'Thermal scope',size=len(_VALIDITY))))
        object.__setattr__(self,'warnings',_warnings(self.warnings))
        stations=tuple(_owned(s,HeatFluxStation,'Heat-flux station') for s in _ordered_items(self.stations,'Heat-flux stations',size=len(source.contour.stations)))
        object.__setattr__(self,'stations',stations)
        expected=_solve(source)
        if any(getattr(self,name)!=value for name,value in expected.items() if name!='station_values') or any(
                any(getattr(s,name)!=value for name,value in wanted.items()) for s,wanted in zip(stations,expected['station_values'])):
            raise DomainError('Thermal result does not replay from owned inputs.',field='Bartz result')


@lru_cache(maxsize=64)
def _solve(inputs):
    """Return every result field of the prescribed profile for ``inputs``.

    The diameter and pressure screening thresholds raise warnings only; they
    are not an experimentally established validity domain of the
    correlation.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """
    values=_profile_exports(inputs)
    warnings=[WarningMessage('bartz_prescribed_model_unqualified',
        'Prescribed constant gas properties and wall temperature; physical applicability, wall/coolant closure and hardware are unqualified.'),
        WarningMessage('thermal_integral_ledger_identity','Zero residual is an integral ledger identity, not measured energy closure.','info')]
    diameter=2*Fraction(inputs.contour.throat_radius_m)
    if not Fraction(.01)<=diameter<=1:
        warnings.append(WarningMessage('bartz_diameter_extrapolation','Outside legacy diameter screening thresholds; these are not an independently validated correlation domain.'))
    if inputs.chamber_pressure_pa>30e6:
        warnings.append(WarningMessage('bartz_pressure_extrapolation','Above the legacy pressure screening threshold; it does not establish physical applicability.'))
    if values['contains_negative_heat_flux']:
        warnings.append(WarningMessage('negative_prescribed_heat_flux','Some prescribed wall locations transfer heat toward the gas; signed flux and integrated load are retained.'))
    values.update(energy_balance_residual_w=0.,correlation=_CORRELATION,property_model=_PROPERTIES,
        validity=_VALIDITY,warnings=tuple(warnings),source=_SOURCE,source_url=_URL,
        physical_model_qualified=False,quadrature_error_bound_w=None,integration_model=_INTEGRATION,
        residual_kind=_RESIDUAL,numerical_policy=_BARTZ_POLICY,pi_basis='fixed_binary64_math_pi')
    return values


def calculate_bartz_heat_transfer(*,contour,chamber_pressure_pa,chamber_temperature_k,wall_temperature_k,
        characteristic_velocity_m_s,dynamic_viscosity_pa_s,specific_heat_j_kg_k,prandtl_number,gamma,throat_radius_of_curvature_m):
    """Evaluate a prescribed Bartz gas-side heat-flux profile along a nozzle contour.

    Parameters
    ----------
    contour : BartzContour, NozzleContour or mapping
        Wall stations with prescribed Mach number and area ratio.
    chamber_pressure_pa : float
        Chamber stagnation pressure, Pa.
    chamber_temperature_k : float
        Chamber stagnation temperature, K.
    wall_temperature_k : float
        Uniform prescribed gas-side wall temperature, K; must be below the
        chamber temperature.
    characteristic_velocity_m_s : float
        Characteristic velocity, m/s.
    dynamic_viscosity_pa_s : float
        Constant gas dynamic viscosity, Pa s.
    specific_heat_j_kg_k : float
        Constant gas specific heat, J/(kg K).
    prandtl_number : float
        Constant Prandtl number.
    gamma : float
        Ratio of specific heats; must exceed one.
    throat_radius_of_curvature_m : float
        Throat wall radius of curvature, m.

    Returns
    -------
    BartzHeatTransferResult
        Station values, peak, signed total heat load, warnings and
        qualification flags. Heat flux is positive into the wall.

    Raises
    ------
    DomainError
        If an input violates the domain checks of :class:`BartzInputs` or
        :class:`BartzContour`, or an exported value is outside the nonzero
        finite binary64 range.
    ConvergenceError
        If a rounding, sign or peak-ordering decision cannot be resolved
        within the numerical precision budget.

    Notes
    -----
    See :class:`BartzHeatTransferResult` for the governing relations. Gas
    properties are user-supplied constants rather than an equilibrium
    calculation, and wall temperature is prescribed rather than solved.
    Negative local flux is retained, not clipped to zero. Diameters outside
    0.01-1 m and chamber pressures above 30 MPa add screening warnings.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """
    inputs=BartzInputs(contour,chamber_pressure_pa,chamber_temperature_k,wall_temperature_k,
        characteristic_velocity_m_s,dynamic_viscosity_pa_s,specific_heat_j_kg_k,prandtl_number,gamma,throat_radius_of_curvature_m)
    values=dict(_solve(inputs));stations=tuple(HeatFluxStation(**s) for s in values.pop('station_values'))
    return BartzHeatTransferResult(**values,stations=stations,inputs=inputs)
