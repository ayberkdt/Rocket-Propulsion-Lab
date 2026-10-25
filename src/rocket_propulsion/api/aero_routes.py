"""HTTP routes for real-fluid thermodynamics and external aerodynamics.

Handlers translate JSON to library calls and sample plot series; every number
comes from :mod:`rocket_propulsion.thermodynamics.real_fluids`,
:mod:`rocket_propulsion.aerodynamics` or the advanced
:mod:`rocket_propulsion.compressible` modules, which document their equations
and sources. Non-finite values are returned as ``null`` so responses are
strict JSON.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson, "Fundamentals of Aerodynamics",
       6th ed., McGraw-Hill, 2017.
.. [2] USSA-1976. NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976".
       https://ntrs.nasa.gov/citations/19770009539
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from math import cos, isfinite, pi, radians, sin, sqrt
from typing import Any, TypedDict

from rocket_propulsion.aerodynamics import (
    airfoil_geometry,
    boundary_layer,
    entry,
    launch_vehicle,
    lifting_line,
    panel_methods,
    slender_body,
    thin_airfoil,
)
from rocket_propulsion.compressible import cone_flow, linearized_flow, shock_expansion
from rocket_propulsion.compressible.off_design_nozzle import standard_atmosphere_pressure_pa
from rocket_propulsion.core.constants import AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import ModelMetadata
from rocket_propulsion.thermodynamics.real_fluids import (
    REFERENCE_FLUID_KEYS,
    cubic_counterpart,
    reference_fluid,
)

AIR_GAMMA = 1.4  # calorically perfect air, Anderson Ch. 7
TAUBER_STAGNATION_CONSTANT = 1.83e-4  # TAUBER-1989-TP2914 Eq. (40), SI, hemisphere, cold wall


def _clean(value: Any) -> Any:
    """Return JSON-safe data: dataclasses to dicts, non-finite floats to ``None``."""
    if is_dataclass(value) and not isinstance(value, type):
        value = asdict(value)
    if isinstance(value, dict):
        return {str(key): _clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(item) for item in value]
    if isinstance(value, float) and not isfinite(value):
        return None
    if hasattr(value, "value") and type(value).__module__ != "builtins" and not isinstance(value, (int, float, str)):
        return value.value
    return value


def _number(payload: dict[str, Any], key: str, default: float) -> float:
    raw = payload.get(key, default)
    if raw in (None, ""):
        raw = default
    if isinstance(raw, bool):
        raise InputError(f"Field {key} must be a number.", field=key)
    try:
        value = float(raw)
    except (TypeError, ValueError, OverflowError) as error:
        raise InputError(f"Field {key} must be a number.", field=key) from error
    if not isfinite(value):
        raise InputError(f"Field {key} must be finite.", field=key)
    return value


def _positive(payload: dict[str, Any], key: str, default: float) -> float:
    """Return ``_number`` for a length, area or flux that must be strictly positive."""
    value = _number(payload, key, default)
    if value <= 0.0:
        raise InputError(f"Field {key} must be greater than zero.", field=key)
    return value


def _integer(payload: dict[str, Any], key: str, default: int, minimum: int, maximum: int) -> int:
    value = _number(payload, key, float(default))
    if not value.is_integer() or not minimum <= value <= maximum:
        raise InputError(f"Field {key} must be an integer from {minimum} to {maximum}.", field=key)
    return int(value)


def _text(payload: dict[str, Any], key: str, default: str) -> str:
    raw = payload.get(key, default)
    if not isinstance(raw, str) or not raw.strip():
        raise InputError(f"Field {key} must be a non-empty string.", field=key)
    return raw.strip()


def _linspace(start: float, stop: float, count: int) -> list[float]:
    return [start + (stop - start) * index / (count - 1) for index in range(count)]


def _atmosphere(geometric_altitude_m: float) -> dict[str, float]:
    """U.S. Standard Atmosphere 1976 state at a geometric altitude below 86 km."""
    pressure = standard_atmosphere_pressure_pa(entry.ussa1976_geopotential_altitude(geometric_altitude_m))
    density = entry.ussa1976_density(geometric_altitude_m)
    temperature = pressure / (density * AIR_SPECIFIC_GAS_CONSTANT_J_KG_K)
    return {
        "geometric_altitude_m": geometric_altitude_m,
        "pressure_pa": pressure,
        "density_kg_m3": density,
        "temperature_k": temperature,
        "speed_of_sound_m_s": sqrt(AIR_GAMMA * AIR_SPECIFIC_GAS_CONSTANT_J_KG_K * temperature),
    }


def real_fluid_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Reference-EOS state, isobar and saturation curve of a propellant [1]_.

    Parameters
    ----------
    payload : dict
        fluid : str, default "oxygen"
            One of the bundled reference fluids.
        temperature_k, pressure_pa : float, default 90 K, 1 MPa
            State to evaluate.
        point_count : int, default 61
            Isobar sample count (11-201).

    Returns
    -------
    dict
        ``state``, the Peng-Robinson ``cubic_state`` at the same point when it
        solves with its relative density difference, ``isobar`` rows (failed
        points omitted) and ``saturation`` rows.

    Raises
    ------
    InputError, DomainError
        For unknown fluids or states outside the equation's range.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID-RF-2026. NIST Chemistry WebBook fluid tables used to
           verify the bundled equations of state.
    """

    key = _text(payload, "fluid", "oxygen")
    if key not in REFERENCE_FLUID_KEYS:
        raise InputError("Unknown fluid.", field="fluid", details={"supported": list(REFERENCE_FLUID_KEYS)})
    fluid = reference_fluid(key)
    model = fluid.model
    temperature = _number(payload, "temperature_k", 90.0)
    pressure = _number(payload, "pressure_pa", 1.0e6)
    count = _integer(payload, "point_count", 61, 11, 201)
    state = fluid.state_tp(temperature, pressure)
    cubic_state = None
    try:
        cubic_state = cubic_counterpart(key, "peng_robinson").state_tp(temperature, pressure)
    except (DomainError, InputError, ArithmeticError):
        cubic_state = None
    isobar = []
    upper = min(model.t_max_k, max(4.0 * model.t_crit_k, temperature * 2.0))
    for value in _linspace(model.t_min_k * 1.001, upper, count):
        try:
            point = fluid.state_tp(value, pressure)
        except (DomainError, ArithmeticError):
            continue
        isobar.append({"temperature_k": value, "phase": point.phase, "density_kg_m3": point.density_kg_m3,
                       "cp_j_kg_k": point.cp_j_kg_k, "speed_of_sound_m_s": point.speed_of_sound_m_s,
                       "compressibility_factor": point.compressibility_factor})
    saturation = []
    for value in _linspace(model.t_min_k * 1.001, model.t_crit_k * 0.995, 30):
        try:
            sat = fluid.saturation_t(value)
        except (DomainError, ArithmeticError):
            continue
        saturation.append({"temperature_k": value, "pressure_pa": sat.pressure_pa,
                           "liquid_density_kg_m3": sat.liquid.density_kg_m3,
                           "vapor_density_kg_m3": sat.vapor.density_kg_m3,
                           "enthalpy_of_vaporization_j_kg": sat.enthalpy_of_vaporization_j_kg})
    return _clean({
        "fluid": key,
        "critical_temperature_k": model.t_crit_k,
        "critical_pressure_pa": model.p_crit_pa,
        "validity": {"t_min_k": model.t_min_k, "t_max_k": model.t_max_k, "p_max_pa": model.p_max_pa},
        "state": state,
        "cubic_state": cubic_state,
        "cubic_relative_density_difference": (
            None if cubic_state is None else cubic_state.density_kg_m3 / state.density_kg_m3 - 1.0
        ),
        "isobar": isobar,
        "saturation": saturation,
        "warnings": [warning for warning in state.warnings],
    })


def airfoil_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Thin-airfoil and lifting-panel solutions for a NACA section [1]_.

    Parameters
    ----------
    payload : dict
        designation : str, default "2412"
        angle_of_attack_deg : float, default 4
        panels_per_surface : int, default 50 (20-120)

    Returns
    -------
    dict
        ``thin`` and ``panel`` results at the angle, ``cp`` surface rows,
        ``geometry`` rows and a lift curve ``sweep`` for both methods.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson, "Fundamentals of Aerodynamics",
           6th ed., 2017, Ch. 4 (thin airfoil theory, panel methods).
    """

    section = airfoil_geometry.naca_section(_text(payload, "designation", "2412"))
    alpha = _number(payload, "angle_of_attack_deg", 4.0)
    per_surface = _integer(payload, "panels_per_surface", 50, 20, 120)
    nodes = section.panel_nodes(per_surface)
    thin = thin_airfoil.thin_airfoil_naca(section, radians(alpha))
    panel = panel_methods.solve_lifting_panels(nodes, radians(alpha))
    cp = [{"x_over_c": (start[0] + end[0]) / 2.0, "cp": value, "upper": (start[1] + end[1]) / 2.0 >= 0.0}
          for start, end, value in zip(nodes, nodes[1:], panel.pressure_coefficients)]
    sweep = []
    for angle in range(-4, 13, 2):
        sweep.append({
            "angle_of_attack_deg": float(angle),
            "thin_lift_coefficient": thin_airfoil.thin_airfoil_naca(section, radians(angle)).lift_coefficient,
            "panel_lift_coefficient": panel_methods.solve_lifting_panels(nodes, radians(angle)).lift_coefficient,
        })
    return _clean({
        "section": section,
        "angle_of_attack_deg": alpha,
        "thin": {key: value for key, value in asdict(thin).items() if key != "fourier_coefficients"},
        "panel": {"lift_coefficient": panel.lift_coefficient,
                  "pressure_lift_coefficient": panel.pressure_lift_coefficient,
                  "pressure_drag_coefficient": panel.pressure_drag_coefficient,
                  "moment_coefficient": panel.moment_coefficient, "panel_count": len(nodes) - 1},
        "cp": cp,
        "geometry": [{"x": x, "y": y} for x, y in nodes],
        "sweep": sweep,
    })


def wing_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Prandtl lifting-line solution of a trapezoidal wing [1]_.

    Parameters
    ----------
    payload : dict
        aspect_ratio : float, default 8
        taper_ratio : float, default 0.4
        angle_of_attack_deg : float, default 5
        tip_twist_deg : float, default 0
        section_lift_slope_per_rad : float, default 2π (thin airfoil)

    Returns
    -------
    dict
        Wing coefficients, spanwise rows and a taper-ratio ``sweep`` of δ.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson, "Fundamentals of Aerodynamics",
           6th ed., 2017, Ch. 5 (lifting-line theory).
    """

    aspect = _number(payload, "aspect_ratio", 8.0)
    taper = _number(payload, "taper_ratio", 0.4)
    alpha = radians(_number(payload, "angle_of_attack_deg", 5.0))
    twist = radians(_number(payload, "tip_twist_deg", 0.0))
    slope = _number(payload, "section_lift_slope_per_rad", 2.0 * pi)
    if aspect <= 0.0 or not 0.0 < taper <= 1.0:
        raise DomainError("Aspect ratio must be positive and taper in (0, 1].", field="taper_ratio")
    root = 1.0
    span = aspect * root * (1.0 + taper) / 2.0
    wing = lifting_line.trapezoidal_wing(span, root, taper, tip_twist_rad=twist, section_lift_slope=slope)
    result = lifting_line.solve_lifting_line(wing, alpha)
    spanwise = [{"y_over_semispan": y / (span / 2.0), "section_lift_coefficient": cl, "circulation": gamma}
                for y, cl, gamma in zip(result.stations_y, result.section_lift_coefficients, result.circulation)]
    sweep = []
    for value in _linspace(0.1, 1.0, 10):
        taper_wing = lifting_line.trapezoidal_wing(aspect * (1.0 + value) / 2.0, 1.0, value, section_lift_slope=slope)
        solution = lifting_line.solve_lifting_line(taper_wing, alpha)
        sweep.append({"taper_ratio": value, "induced_drag_factor": solution.induced_drag_factor,
                      "span_efficiency": solution.span_efficiency})
    return _clean({
        "aspect_ratio": result.aspect_ratio, "lift_coefficient": result.lift_coefficient,
        "induced_drag_coefficient": result.induced_drag_coefficient,
        "induced_drag_factor": result.induced_drag_factor, "span_efficiency": result.span_efficiency,
        "spanwise": spanwise, "sweep": sweep,
    })


def supersonic_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Taylor-Maccoll cone flow and supersonic diamond-airfoil coefficients [1]_.

    Parameters
    ----------
    payload : dict
        mach : float, default 2
        cone_half_angle_deg : float, default 10
        airfoil_half_angle_deg : float, default 5
        angle_of_attack_deg : float, default 2

    Returns
    -------
    dict
        ``cone`` result, diamond-airfoil ``shock_expansion``, ``ackeret`` and
        ``busemann`` coefficients and an angle-of-attack ``sweep``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson, "Fundamentals of Aerodynamics",
           6th ed., 2017, Ch. 9, 12 and 13.
    """

    mach = _number(payload, "mach", 2.0)
    cone = cone_flow.solve_cone_flow(mach, _number(payload, "cone_half_angle_deg", 10.0))
    section = shock_expansion.diamond_airfoil(_number(payload, "airfoil_half_angle_deg", 5.0))
    alpha = _number(payload, "angle_of_attack_deg", 2.0)

    def coefficients(angle: float) -> dict[str, Any]:
        exact = shock_expansion.shock_expansion_coefficients(section, mach, angle)
        first = linearized_flow.linearized_section_coefficients(section, mach, angle, order=1)
        second = linearized_flow.linearized_section_coefficients(section, mach, angle, order=2)
        return {"angle_of_attack_deg": angle,
                "shock_expansion_lift": exact.lift_coefficient, "shock_expansion_drag": exact.drag_coefficient,
                "ackeret_lift": first.lift_coefficient, "ackeret_drag": first.drag_coefficient,
                "busemann_lift": second.lift_coefficient, "busemann_drag": second.drag_coefficient}

    sweep = []
    for angle in _linspace(0.0, 10.0, 11):
        try:
            sweep.append(coefficients(angle))
        except DomainError:
            break
    return _clean({"mach": mach, "cone": cone, "at_angle": coefficients(alpha), "sweep": sweep})


class _FlatPlateEdge(TypedDict):
    """Edge state shared by every flat-plate reference-temperature call."""

    edge_temperature_k: float
    edge_pressure_pa: float
    edge_velocity_m_s: float
    gas_constant_j_kg_k: float
    wall_temperature_k: float | None


def flat_plate_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Compressible flat-plate skin friction and heating at a flight condition [1]_.

    Parameters
    ----------
    payload : dict
        mach : float, default 2
        altitude_m : float, default 15000 (geometric, below 86 km)
        length_m : float, default 5
        wall_temperature_k : float, optional (adiabatic wall when absent)

    Returns
    -------
    dict
        ``atmosphere``, laminar and turbulent reference-temperature results,
        Van Driest II at the plate end, and an ``x`` sweep of local cf.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson, "Fundamentals of Aerodynamics",
           6th ed., 2017, Ch. 18-19.
    .. [2] HOPKINS-1972-TND6945. E. J. Hopkins and M. Inouye, NASA TN D-6945
           (Van Driest II). https://ntrs.nasa.gov/citations/19730001588
    """

    mach = _number(payload, "mach", 2.0)
    altitude = _number(payload, "altitude_m", 15_000.0)
    length = _number(payload, "length_m", 5.0)
    raw_wall = payload.get("wall_temperature_k")
    wall = None if raw_wall in (None, "") else _number(payload, "wall_temperature_k", 0.0)
    air = _atmosphere(altitude)
    velocity = mach * air["speed_of_sound_m_s"]
    common: _FlatPlateEdge = {
        "edge_temperature_k": air["temperature_k"], "edge_pressure_pa": air["pressure_pa"],
        "edge_velocity_m_s": velocity, "gas_constant_j_kg_k": AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
        "wall_temperature_k": wall}
    laminar = boundary_layer.compressible_flat_plate_reference_temperature(x_m=length, regime="laminar", **common)
    turbulent = boundary_layer.compressible_flat_plate_reference_temperature(x_m=length, regime="turbulent", **common)
    reynolds = air["density_kg_m3"] * velocity * length / boundary_layer.SutherlandViscosity().viscosity(air["temperature_k"])
    van_driest = boundary_layer.van_driest_ii(
        mach=mach, edge_temperature_k=air["temperature_k"], reynolds_x=reynolds, wall_temperature_k=wall,
        wall_to_adiabatic_temperature=None if wall is not None else 1.0)
    sweep = []
    for x in _linspace(length / 40.0, length, 40):
        sweep.append({
            "x_m": x,
            "laminar_cf": boundary_layer.compressible_flat_plate_reference_temperature(
                x_m=x, regime="laminar", **common).local_skin_friction,
            "turbulent_cf": boundary_layer.compressible_flat_plate_reference_temperature(
                x_m=x, regime="turbulent", **common).local_skin_friction,
        })
    return _clean({"atmosphere": air, "velocity_m_s": velocity, "reynolds_length": reynolds,
                   "laminar": laminar, "turbulent": turbulent, "van_driest_ii": van_driest, "sweep": sweep})


def entry_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Planar point-mass atmospheric entry with stagnation heating [1]_, [2]_.

    Parameters
    ----------
    payload : dict
        entry_velocity_m_s : float, default 7000
        entry_flight_path_angle_deg : float, default -10
        mass_kg, drag_coefficient, reference_area_m2 : float, default 1500, 1.0, 4.0
        nose_radius_m : float, default 1.0
        lift_to_drag : float, default 0

    Returns
    -------
    dict
        Peak deceleration and heating, heat load and a decimated trajectory
        (altitude, velocity, deceleration, heating vs time).

    References
    ----------
    .. [1] CHAPMAN-1958-TN4276. D. R. Chapman, NACA TN 4276 (planar entry equations).
           https://ntrs.nasa.gov/citations/19930085059
    .. [2] TAUBER-1989-TP2914. M. E. Tauber, NASA TP-2914, Eq. (40).
           https://ntrs.nasa.gov/citations/19890017745
    """

    trajectory = entry.integrate_planar_entry(
        entry_altitude_m=_number(payload, "entry_altitude_m", 85_000.0),
        entry_velocity_m_s=_number(payload, "entry_velocity_m_s", 7_000.0),
        entry_flight_path_angle_deg=_number(payload, "entry_flight_path_angle_deg", -10.0),
        mass_kg=_number(payload, "mass_kg", 1_500.0),
        drag_coefficient=_number(payload, "drag_coefficient", 1.0),
        reference_area_m2=_number(payload, "reference_area_m2", 4.0),
        lift_to_drag=_number(payload, "lift_to_drag", 0.0),
        nose_radius_m=_number(payload, "nose_radius_m", 1.0),
        stagnation_heating_constant=TAUBER_STAGNATION_CONSTANT,
        time_step_s=0.1,
        final_velocity_m_s=100.0,
    )
    samples = trajectory.samples
    step = max(1, len(samples) // 300)
    rows = [asdict(sample) for sample in samples[::step]]
    if rows[-1]["time_s"] != samples[-1].time_s:
        rows.append(asdict(samples[-1]))
    return _clean({
        "termination": trajectory.termination,
        "peak_deceleration": trajectory.peak_deceleration,
        "peak_heating": trajectory.peak_heating,
        "integrated_heat_load_j_m2": trajectory.integrated_heat_load_j_m2,
        "trajectory": rows,
        "warnings": trajectory.warnings,
    })


def launch_vehicle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Slender-body normal force, viscous crossflow and drag build-up of an ascent body [1]_.

    Parameters
    ----------
    payload : dict
        nose_shape : {"tangent_ogive", "cone"}, default "tangent_ogive"
        nose_length_m, body_length_m, diameter_m : float, default 6, 40, 3.7
        mach : float, default 0.8
        altitude_m : float, default 11000
        angle_of_attack_deg : float, default 4
        crossflow_drag_coefficient : float, default 1.2

    Returns
    -------
    dict
        Geometry, ``normal_force`` (CNα, centre of pressure), ``crossflow`` at
        the angle, ``drag`` components, dynamic pressure and an angle sweep.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, "The Practical Calculation of the
           Aerodynamic Characteristics of Slender Finned Vehicles", 1967.
    .. [2] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951.
           https://ntrs.nasa.gov/citations/19930090962
    """

    shape = _text(payload, "nose_shape", "tangent_ogive")
    if shape not in {"tangent_ogive", "cone"}:
        raise InputError("Nose shape must be tangent_ogive or cone.", field="nose_shape")
    nose_length = _number(payload, "nose_length_m", 6.0)
    body_length = _number(payload, "body_length_m", 40.0)
    radius = _number(payload, "diameter_m", 3.7) / 2.0
    mach = _number(payload, "mach", 0.8)
    altitude = _number(payload, "altitude_m", 11_000.0)
    alpha = _number(payload, "angle_of_attack_deg", 4.0)
    crossflow_cd = _number(payload, "crossflow_drag_coefficient", 1.2)
    nose = (slender_body.tangent_ogive_nose if shape == "tangent_ogive" else slender_body.cone_nose)(nose_length, radius)
    body = slender_body.AxisymmetricBody((nose, slender_body.cylinder_section(body_length, radius)))
    geometry = slender_body.body_geometry(body)
    area = pi * radius * radius
    normal = slender_body.combine_normal_force_components(
        slender_body.slender_body_normal_force(body, reference_area_m2=area))
    crossflow = slender_body.viscous_crossflow_coefficients(
        body, angle_of_attack_deg=alpha, crossflow_drag_coefficient=crossflow_cd,
        reference_area_m2=area, reference_length_m=2.0 * radius)
    air = _atmosphere(altitude)
    velocity = mach * air["speed_of_sound_m_s"]
    total_length = geometry.length_m
    reynolds = air["density_kg_m3"] * velocity * total_length / boundary_layer.SutherlandViscosity().viscosity(air["temperature_k"])
    friction = launch_vehicle.barrowman_skin_friction(reynolds_number=reynolds, mach=mach)
    sonic = launch_vehicle.barrowman_skin_friction(reynolds_number=reynolds, mach=1.0)
    friction_drag = launch_vehicle.body_friction_drag_coefficient(
        compressible_skin_friction=friction.compressible_coefficient, body_fineness_ratio=geometry.fineness_ratio,
        wetted_area_m2=geometry.wetted_area_m2, reference_area_m2=area)
    base = launch_vehicle.barrowman_base_drag(
        mach=mach, base_area_m2=geometry.base_area_m2, reference_area_m2=area,
        wetted_area_m2=geometry.wetted_area_m2, compressible_skin_friction=friction.compressible_coefficient,
        sonic_skin_friction=sonic.compressible_coefficient)
    forebody = None
    if mach < 1.0:
        forebody = launch_vehicle.subsonic_forebody_pressure_drag_coefficient(
            mach=mach, nose_shape=shape, nose_fineness_ratio=nose_length / (2.0 * radius),
            body_fineness_ratio=geometry.fineness_ratio, wetted_area_m2=geometry.wetted_area_m2,
            reference_area_m2=area, compressible_skin_friction=friction.compressible_coefficient,
            sonic_skin_friction=sonic.compressible_coefficient)
    sweep = []
    for angle in _linspace(0.0, 12.0, 13):
        point = slender_body.viscous_crossflow_coefficients(
            body, angle_of_attack_deg=angle, crossflow_drag_coefficient=crossflow_cd,
            reference_area_m2=area, reference_length_m=2.0 * radius)
        sweep.append({"angle_of_attack_deg": angle, "normal_force_coefficient": point.normal_force_coefficient,
                      "potential_normal_force_coefficient": point.potential_normal_force_coefficient})
    dynamic_pressure = launch_vehicle.dynamic_pressure_at_altitude(
        geopotential_altitude_m=entry.ussa1976_geopotential_altitude(altitude), mach=mach)
    normal_force = crossflow.normal_force_coefficient * dynamic_pressure * area
    return _clean({
        "geometry": {key: value for key, value in asdict(geometry).items() if key != "segments"},
        "reference_area_m2": area,
        "normal_force": normal,
        "crossflow": crossflow,
        "drag": {"skin_friction": friction, "friction_drag_coefficient": friction_drag,
                 "base_drag_coefficient": base.base_drag_coefficient,
                 "forebody_pressure_drag_coefficient": forebody, "base": base},
        "atmosphere": air,
        "velocity_m_s": velocity,
        "dynamic_pressure_pa": dynamic_pressure,
        "normal_force_n": normal_force,
        "axial_load_direction": {"cos_alpha": cos(radians(alpha)), "sin_alpha": sin(radians(alpha))},
        "sweep": sweep,
    })


AERO_ROUTES = {
    "/api/v1/real-fluid": real_fluid_route,
    "/api/v1/airfoil": airfoil_route,
    "/api/v1/wing": wing_route,
    "/api/v1/supersonic-body": supersonic_route,
    "/api/v1/flat-plate": flat_plate_route,
    "/api/v1/entry": entry_route,
    "/api/v1/launch-vehicle-aero": launch_vehicle_route,
}

_REFERENCES = ("Anderson, Fundamentals of Aerodynamics, 6th ed., 2017",)
AERO_ROUTE_METADATA = {
    "/api/v1/real-fluid": ModelMetadata(
        model="reference Helmholtz equations of state", version="1.0",
        assumptions=("Single-component reference EOS within the published T/p range",
                     "Verified against NIST WebBook within 11 ppm"),
        units={"temperature": "K", "pressure": "Pa", "density": "kg/m^3"},
        references=("Leachman 2009; Stewart 1991; Span 2000; Setzmann-Wagner 1991; NIST IR 8474",)),
    "/api/v1/airfoil": ModelMetadata(
        model="thin airfoil theory and lifting panel method", version="1.0",
        assumptions=("Inviscid incompressible flow", "Kutta condition at the trailing edge"),
        units={"angle": "deg", "coefficients": "1"}, references=_REFERENCES + ("NACA Report 824",)),
    "/api/v1/wing": ModelMetadata(
        model="Prandtl lifting line (Glauert series)", version="1.0",
        assumptions=("Unswept high-aspect-ratio wing", "Linear section lift"),
        units={"angle": "deg", "coefficients": "1"}, references=_REFERENCES),
    "/api/v1/supersonic-body": ModelMetadata(
        model="Taylor-Maccoll cone flow and shock-expansion theory", version="1.0",
        assumptions=("Calorically perfect air, gamma 1.4", "Attached shocks"),
        units={"angle": "deg", "coefficients": "1"}, references=_REFERENCES + ("NASA SP-3004",)),
    "/api/v1/flat-plate": ModelMetadata(
        model="reference-temperature and Van Driest II flat plate", version="1.0",
        assumptions=("Zero pressure gradient", "Transition not predicted; laminar and turbulent shown"),
        units={"skin_friction": "1", "heat_flux": "W/m^2"}, references=_REFERENCES + ("NASA TN D-6945",)),
    "/api/v1/entry": ModelMetadata(
        model="planar point-mass entry with Tauber stagnation heating", version="1.0",
        assumptions=("Non-rotating spherical Earth, USSA-1976 below 86 km", "Constant CD, cold-wall heating"),
        units={"velocity": "m/s", "heating": "W/m^2"}, references=("NACA TN 4276", "NASA TP-2914")),
    "/api/v1/launch-vehicle-aero": ModelMetadata(
        model="slender-body, viscous crossflow and Barrowman drag build-up", version="1.0",
        assumptions=("Axisymmetric nose-cylinder body without fins", "Small sideslip; no plume or protuberances"),
        units={"force": "N", "pressure": "Pa", "coefficients": "1"},
        references=("Barrowman 1967", "NACA Report 1048", "NASA SP-8035")),
}
