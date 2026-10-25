"""Two-dimensional conduction in the wall cross-section of a cooled channel.

The fin model of :mod:`.regen_chamber` treats the land between two channels
as a straight fin on a hot wall of uniform temperature. Real walls are hotter
above the channel centre than above the land, and the peak hot-wall
temperature (what limits life) sits there. This module solves steady 2-D
conduction in one symmetric half pitch of the wall cross-section by the
finite-volume method:

* hot wall ``0 <= z <= t`` across the half pitch ``0 <= y <= W/2`` (``y = 0`` on
  the channel centreline, ``y = W/2`` on the land centreline);
* land ``t <= z <= t + h`` over ``w/2 <= y <= W/2`` (the channel occupies
  ``0 <= y < w/2``);
* hot face ``z = 0``: ``q = h_g (T_aw - T) + q_rad`` (Robin plus an imposed flux);
* channel floor ``z = t, y < w/2`` and side ``y = w/2, t < z < t + h``:
  ``q = h_c (T - T_b)``;
* symmetry planes and the land top (bonded to the closeout) are adiabatic.

Fourier conduction with constant ``k`` and the convective (Robin) boundary
condition are those of [1]_. Boundary conductances combine the half-cell
conduction with the film: ``G = A / (dz/(2k) + 1/h)``. The symmetric
positive-definite system is solved by banded Cholesky factorization, rows
ordered so the bandwidth is the number of cells across the pitch.

References
----------
.. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
       Transfer Textbook", v6.00, 2024, Ch. 1-2 (Fourier's law, convective
       boundary condition, conduction equation).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import fsum, isfinite, sqrt

from rocket_propulsion.core.errors import DomainError, InputError

__all__ = ["WallSection", "solve_channel_wall_section"]

#: Largest accepted grid (library choice: a pure-Python banded solve in well under a second).
MAXIMUM_CELLS = 20_000


@dataclass(frozen=True, slots=True)
class WallSection:
    """Solved temperature field of one half pitch of the channel wall.

    Attributes
    ----------
    peak_hot_wall_temperature_k, minimum_hot_wall_temperature_k : float
        Extremes of the hot-surface temperature ``T_aw - q/h_g`` (radiation
        excluded) across the pitch.
    mean_hot_wall_temperature_k : float
        Area-weighted mean hot-surface temperature.
    mean_coolant_wall_temperature_k : float
        Area-weighted mean temperature of the wetted surface (floor and side).
    heat_per_length_w_m : float
        Heat entering the hot face per unit channel length, per half pitch.
    mean_heat_flux_w_m2 : float
        ``heat_per_length / (W/2)``, the flux per unit hot-wall area.
    energy_residual_w_m : float
        Heat in minus heat to the coolant; round-off for a converged solve.
    peak_location_y_m : float
        Distance of the hottest hot-surface cell centre from the channel
        centreline.
    hot_surface : tuple of (float, float)
        ``(y, T_surface)`` for each hot-face cell.
    """

    peak_hot_wall_temperature_k: float
    minimum_hot_wall_temperature_k: float
    mean_hot_wall_temperature_k: float
    mean_coolant_wall_temperature_k: float
    heat_per_length_w_m: float
    mean_heat_flux_w_m2: float
    energy_residual_w_m: float
    peak_location_y_m: float
    hot_surface: tuple[tuple[float, float], ...]


def _positive(value: float, name: str) -> float:
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise DomainError(f"{name} must be finite and positive.", field=name)
    return number


def _cells(count: object, name: str) -> int:
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise InputError(f"{name} must be a positive integer.", field=name)
    return count


def _banded_cholesky_solve(diagonal: list[float], bands: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve ``A x = b`` for symmetric positive-definite banded ``A``.

    ``bands[k][i]`` holds ``A[i, i + k + 1]``. Factorization ``A = L L^T`` in
    band storage, then forward and back substitution.
    """
    size, width = len(diagonal), len(bands)
    # lower[i][k] = L[i, i - k - 1]
    lower = [[0.0] * width for _ in range(size)]
    pivot = [0.0] * size
    for i in range(size):
        for k in range(min(width, i) - 1, -1, -1):
            j = i - k - 1
            total = bands[k][j]
            for m in range(k + 1, min(width, j + k + 1)):
                # L[i, j - (m - k)] * L[j, j - (m - k)]
                total -= lower[i][m] * lower[j][m - k - 1]
            lower[i][k] = total / pivot[j]
        value = diagonal[i] - fsum(entry * entry for entry in lower[i])
        if not value > 0.0:
            raise DomainError("Wall conduction matrix is not positive definite.", code="wall_matrix_singular")
        pivot[i] = sqrt(value)
    forward = [0.0] * size
    for i in range(size):
        total = rhs[i] - fsum(lower[i][k] * forward[i - k - 1] for k in range(min(width, i)))
        forward[i] = total / pivot[i]
    solution = [0.0] * size
    for i in range(size - 1, -1, -1):
        total = forward[i] - fsum(lower[i + k + 1][k] * solution[i + k + 1]
                                  for k in range(min(width, size - 1 - i)))
        solution[i] = total / pivot[i]
    return solution


def solve_channel_wall_section(*, pitch_m: float, channel_width_m: float, channel_height_m: float,
                               hot_wall_thickness_m: float, thermal_conductivity_w_m_k: float,
                               gas_coefficient_w_m2_k: float, recovery_temperature_k: float,
                               coolant_coefficient_w_m2_k: float, coolant_temperature_k: float,
                               radiative_heat_flux_w_m2: float = 0.0, cells_channel: int = 12,
                               cells_land: int = 8, cells_wall: int = 8, cells_height: int = 16) -> WallSection:
    """Solve steady 2-D conduction in a half pitch of a channel wall [1]_.

    Parameters
    ----------
    pitch_m : float
        Channel pitch ``W`` at the hot wall (channel plus land).
    channel_width_m, channel_height_m : float
        Channel width ``w < W`` and depth ``h``.
    hot_wall_thickness_m : float
        Hot-wall thickness ``t`` under the channel.
    thermal_conductivity_w_m_k : float
        Wall conductivity ``k``.
    gas_coefficient_w_m2_k, recovery_temperature_k : float
        Hot-gas film coefficient and adiabatic-wall temperature.
    coolant_coefficient_w_m2_k, coolant_temperature_k : float
        Coolant film coefficient and bulk temperature.
    radiative_heat_flux_w_m2 : float, default 0
        Gas radiation absorbed by the hot face (imposed, non-negative).
    cells_channel, cells_land, cells_wall, cells_height : int
        Cells across the half channel, the half land, the hot wall and the
        land height.

    Returns
    -------
    WallSection
        Surface temperatures, mean and peak, heat per length and residual.

    Raises
    ------
    DomainError
        For non-positive or non-finite inputs, ``w >= W``, a negative
        radiative flux or a coolant hotter than the recovery temperature.
    InputError
        For cell counts that are not positive integers or a grid above
        20 000 cells.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard, Ch. 1-2.
    """

    pitch = _positive(pitch_m, "pitch_m")
    width = _positive(channel_width_m, "channel_width_m")
    height = _positive(channel_height_m, "channel_height_m")
    thickness = _positive(hot_wall_thickness_m, "hot_wall_thickness_m")
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    h_gas = _positive(gas_coefficient_w_m2_k, "gas_coefficient_w_m2_k")
    h_cool = _positive(coolant_coefficient_w_m2_k, "coolant_coefficient_w_m2_k")
    t_aw = _positive(recovery_temperature_k, "recovery_temperature_k")
    t_b = _positive(coolant_temperature_k, "coolant_temperature_k")
    q_rad = float(radiative_heat_flux_w_m2)
    if not isfinite(q_rad) or q_rad < 0.0:
        raise DomainError("Radiative heat flux must be finite and non-negative.", field="radiative_heat_flux_w_m2")
    if width >= pitch:
        raise DomainError("The channel must be narrower than the pitch (a land is required).",
                          field="channel_width_m")
    if t_b > t_aw:
        raise DomainError("Coolant temperature exceeds the recovery temperature.", field="coolant_temperature_k")
    n_ch, n_land = _cells(cells_channel, "cells_channel"), _cells(cells_land, "cells_land")
    n_wall, n_height = _cells(cells_wall, "cells_wall"), _cells(cells_height, "cells_height")
    n_y = n_ch + n_land
    if n_y * (n_wall + n_height) > MAXIMUM_CELLS:
        raise InputError(f"The wall grid exceeds {MAXIMUM_CELLS} cells.", field="cells_height")

    half_channel, half_land = 0.5 * width, 0.5 * (pitch - width)
    dy = [half_channel / n_ch] * n_ch + [half_land / n_land] * n_land
    y_centre = []
    position = 0.0
    for size in dy:
        y_centre.append(position + 0.5 * size)
        position += size
    dz = [thickness / n_wall] * n_wall + [height / n_height] * n_height

    # cell (row j in z, column i in y) exists in the wall layers, and in the land columns above
    index: dict[tuple[int, int], int] = {}
    for j in range(n_wall + n_height):
        for i in range(n_y):
            if j < n_wall or i >= n_ch:
                index[j, i] = len(index)
    size = len(index)
    diagonal = [0.0] * size
    bands = [[0.0] * size for _ in range(n_y)]  # neighbour within one row of cells across the pitch
    rhs = [0.0] * size

    def couple(a: int, b: int, conductance: float) -> None:
        low, high = min(a, b), max(a, b)
        diagonal[low] += conductance
        diagonal[high] += conductance
        bands[high - low - 1][low] -= conductance

    hot_faces: list[tuple[int, float, float]] = []  # (cell, conductance to T_aw, dy)
    coolant_faces: list[tuple[int, float, float]] = []  # (cell, conductance to T_b, area)
    for (j, i), cell in index.items():
        if (j, i + 1) in index:
            couple(cell, index[j, i + 1], k * dz[j] / (0.5 * (dy[i] + dy[i + 1])))
        if (j + 1, i) in index:
            couple(cell, index[j + 1, i], k * dy[i] / (0.5 * (dz[j] + dz[j + 1])))
        if j == 0:
            conductance = dy[i] / (0.5 * dz[j] / k + 1.0 / h_gas)
            diagonal[cell] += conductance
            rhs[cell] += conductance * t_aw + q_rad * dy[i]
            hot_faces.append((cell, conductance, dy[i]))
        if j == n_wall - 1 and i < n_ch:  # channel floor
            conductance = dy[i] / (0.5 * dz[j] / k + 1.0 / h_cool)
            diagonal[cell] += conductance
            rhs[cell] += conductance * t_b
            coolant_faces.append((cell, conductance, dy[i]))
        if j >= n_wall and i == n_ch:  # channel side wall
            conductance = dz[j] / (0.5 * dy[i] / k + 1.0 / h_cool)
            diagonal[cell] += conductance
            rhs[cell] += conductance * t_b
            coolant_faces.append((cell, conductance, dz[j]))

    temperatures = _banded_cholesky_solve(diagonal, bands, rhs)
    gained = [conductance * (t_aw - temperatures[cell]) + q_rad * width_y
              for cell, conductance, width_y in hot_faces]
    heat_in = fsum(gained)
    heat_out = fsum(conductance * (temperatures[cell] - t_b) for cell, conductance, _ in coolant_faces)
    surface = []
    for (cell, conductance, width_y), y in zip(hot_faces, y_centre, strict=True):
        convective = conductance * (t_aw - temperatures[cell]) / width_y
        surface.append((y, t_aw - convective / h_gas))
    peak = max(surface, key=lambda item: item[1])
    wetted = fsum(area for _, _, area in coolant_faces)
    # wetted-surface temperature T_b + q/h_c, area-weighted
    coolant_wall = fsum((t_b + (temperatures[cell] - t_b) * conductance / area / h_cool) * area
                        for cell, conductance, area in coolant_faces) / wetted
    return WallSection(
        peak_hot_wall_temperature_k=peak[1],
        minimum_hot_wall_temperature_k=min(t for _, t in surface),
        mean_hot_wall_temperature_k=fsum(t * d for (_, t), d in zip(surface, dy, strict=True)) / (0.5 * pitch),
        mean_coolant_wall_temperature_k=coolant_wall,
        heat_per_length_w_m=heat_in,
        mean_heat_flux_w_m2=heat_in / (0.5 * pitch),
        energy_residual_w_m=heat_in - heat_out,
        peak_location_y_m=peak[0],
        hot_surface=tuple(surface),
    )
