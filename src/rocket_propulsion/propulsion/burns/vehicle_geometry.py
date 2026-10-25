"""Body-frame thruster geometry and propellant-dependent mass properties.

Translational propagation only needs the net force along a commanded axis.
Attitude-coupled (6-DOF) simulation, momentum unloading and thrust-vector
alignment also need *where* each thruster pushes and *where* the centre of
mass is while propellant drains.  This module supplies both in the vehicle
body frame:

``ThrusterLayout``
    Per-unit-thrust influence matrix (force and torque about the centre of
    mass), resultant force and torque for given thrust levels, cosine (axial)
    efficiency, and the torque derivative with respect to the centre of mass.
``PropellantMassModel``
    Total mass, centre of mass and inertia about the centre of mass from a dry
    rigid body plus point-mass propellant stations whose masses are the
    propagated tank states, with analytic derivatives of the centre of mass
    and inertia with respect to every tank mass.

Attitude dynamics and frame rotation stay in the consuming framework.  All
vectors here are body-frame and SI.  Thruster torque and rigid-body mass
properties follow [1]_; centre-of-mass moment arms and thrust-vector
misalignment during thrusting maneuvers are discussed in [2]_.  Thrust
allocation uses a bounded-variable least-squares active-set method [3]_ [4]_.

References
----------
.. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
.. [2] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
       Maneuvers", https://ntrs.nasa.gov/citations/19710016722.
.. [3] C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems",
       Prentice-Hall, 1974 (NNLS algorithm; reprinted by SIAM, 1995).
.. [4] P. B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: an
       Algorithm and Applications", Computational Statistics 10, 1995.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import sqrt

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import (
    real,
    real_mapping,
    real_matrix,
    records,
    text,
    vector3,
)

Vector = tuple[float, float, float]
Matrix3 = tuple[Vector, Vector, Vector]

_ZERO: Vector = (0.0, 0.0, 0.0)


def _vector(value: tuple[float, ...], label: str) -> Vector:
    """Return a finite three-component vector or raise ``DomainError``."""
    return vector3(value, label)


def _cross(a: Vector, b: Vector) -> Vector:
    """Return the cross product ``a x b``."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vector, b: Vector) -> float:
    """Return the dot product ``a . b``."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub(a: Vector, b: Vector) -> Vector:
    """Return the difference ``a - b``."""
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _skew(a: Vector) -> Matrix3:
    """Return ``[a]x`` such that ``[a]x b = a x b``."""

    return ((0.0, -a[2], a[1]), (a[2], 0.0, -a[0]), (-a[1], a[0], 0.0))


def _shift_operator(r: Vector) -> Matrix3:
    """Return the parallel-axis operator per unit mass, ``(r.r) I - r r^T``.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    rr = _dot(r, r)
    return tuple(
        tuple((rr if i == j else 0.0) - r[i] * r[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def _matrix_add(*terms: tuple[float, Matrix3]) -> Matrix3:
    """Return the weighted sum ``sum(weight * matrix)`` of 3x3 matrices."""
    return tuple(
        tuple(sum(weight * matrix[i][j] for weight, matrix in terms) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class ThrusterMount:
    """One thruster's body-frame application point and force direction.

    Parameters
    ----------
    thruster_id : str
        Non-empty thruster identifier.
    position_m : tuple of float
        Body-frame application point, in m.
    thrust_axis : tuple of float
        Direction of the force *on the vehicle* (opposite to the exhaust);
        non-zero, normalized on construction.

    Raises
    ------
    DomainError
        If the identifier is blank, a vector is not finite and
        three-dimensional, or the axis is zero.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    thruster_id: str
    position_m: Vector
    thrust_axis: Vector

    def __post_init__(self) -> None:
        text(self.thruster_id, "Thruster identifier")
        object.__setattr__(self, "position_m", _vector(self.position_m, "Thruster position"))
        axis = _vector(self.thrust_axis, "Thrust axis")
        norm = sqrt(_dot(axis, axis))
        if norm <= 0.0:
            raise DomainError("Thrust axis must be non-zero.")
        object.__setattr__(self, "thrust_axis", tuple(value / norm for value in axis))


@dataclass(frozen=True, slots=True)
class ForceTorque:
    """Resultant body-frame force and torque about the centre of mass.

    Parameters
    ----------
    force_n : tuple of float
        Net force ``F = sum T_i u_i``, in N.
    torque_n_m : tuple of float
        Net torque ``sum (p_i - c) x T_i u_i`` about the centre of mass, N m.
    total_thrust_n : float
        Scalar sum of thruster thrusts, in N.
    torque_wrt_center_of_mass : tuple of tuple of float
        Jacobian ``d tau / d c = [F]x``, in N.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    force_n: Vector
    torque_n_m: Vector
    total_thrust_n: float
    torque_wrt_center_of_mass: Matrix3


def _least_squares(columns: list[list[float]], target: list[float]) -> list[float] | None:
    """Solve ``min ||A s - b||`` by Householder QR; ``None`` if rank deficient.

    A column is treated as dependent when its remaining norm is at most
    ``1e-12`` times the largest absolute matrix entry (at least one).

    References
    ----------
    .. [1] C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems",
           Prentice-Hall, 1974 (NNLS algorithm; reprinted by SIAM, 1995).
    """

    rows = len(target)
    count = len(columns)
    matrix = [[columns[j][i] for j in range(count)] for i in range(rows)]
    rhs = list(target)
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    for k in range(count):
        norm = sqrt(sum(matrix[i][k] ** 2 for i in range(k, rows)))
        if norm <= 1e-12 * max(scale, 1.0):
            return None
        alpha = -norm if matrix[k][k] >= 0.0 else norm
        v = [0.0] * rows
        v[k] = matrix[k][k] - alpha
        for i in range(k + 1, rows):
            v[i] = matrix[i][k]
        v_norm2 = sum(value * value for value in v[k:])
        if v_norm2 > 0.0:
            for j in range(k, count):
                factor = 2.0 * sum(v[i] * matrix[i][j] for i in range(k, rows)) / v_norm2
                for i in range(k, rows):
                    matrix[i][j] -= factor * v[i]
            factor = 2.0 * sum(v[i] * rhs[i] for i in range(k, rows)) / v_norm2
            for i in range(k, rows):
                rhs[i] -= factor * v[i]
    solution = [0.0] * count
    for k in range(count - 1, -1, -1):
        total = rhs[k] - sum(matrix[k][j] * solution[j] for j in range(k + 1, count))
        solution[k] = total / matrix[k][k]
    return solution


def _bvls(
    columns: list[list[float]],
    target: list[float],
    upper: list[float],
    *,
    tolerance: float,
    max_iterations: int,
) -> tuple[list[float], int]:
    """Bounded-variable least squares ``min ||A x - b||``, ``0 <= x <= upper``.

    Active-set method: Lawson-Hanson NNLS [1]_ generalized to upper bounds in
    the manner of Stark and Parker [2]_.  With every upper bound infinite it
    reduces to NNLS.  ``blocked`` holds variables that could not be freed
    without a rank deficiency or without moving the solution; it is cleared
    whenever the solution changes.

    Returns
    -------
    x : list of float
        Bounded solution.
    iterations : int
        Number of inner least-squares iterations.

    Raises
    ------
    ConvergenceError
        If the inner iteration count exceeds ``max_iterations``.

    References
    ----------
    .. [1] C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems",
           Prentice-Hall, 1974 (NNLS algorithm; reprinted by SIAM, 1995).
    .. [2] P. B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: an
           Algorithm and Applications", Computational Statistics 10, 1995.
    """

    count = len(columns)
    rows = len(target)
    x = [0.0] * count
    free: list[int] = []
    at_upper: set[int] = set()
    blocked: set[int] = set()

    def gradient(values: list[float]) -> list[float]:
        residual = [
            target[i] - sum(columns[j][i] * values[j] for j in range(count)) for i in range(rows)
        ]
        return [sum(columns[j][i] * residual[i] for i in range(rows)) for j in range(count)]

    def bound_tolerance(j: int) -> float:
        return 1e-12 * max(1.0, upper[j] if upper[j] != float("inf") else 1.0)

    iterations = 0
    w = gradient(x)
    while True:
        candidates = [
            j
            for j in range(count)
            if j not in free
            and j not in blocked
            and ((j not in at_upper and w[j] > tolerance) or (j in at_upper and w[j] < -tolerance))
        ]
        if not candidates:
            return x, iterations
        entering = max(candidates, key=lambda j: abs(w[j]))
        entered_from_upper = entering in at_upper
        at_upper.discard(entering)
        free.append(entering)
        previous = list(x)
        while True:
            iterations += 1
            if iterations > max_iterations:
                raise ConvergenceError("Thrust allocation did not converge.")
            fixed = [j for j in range(count) if j not in free]
            rhs = [
                target[i] - sum(columns[j][i] * x[j] for j in fixed) for i in range(rows)
            ]
            solved = _least_squares([columns[j] for j in free], rhs)
            if solved is None:
                free.remove(entering)
                if entered_from_upper:
                    at_upper.add(entering)
                blocked.add(entering)
                break
            trial = list(x)
            for position, j in enumerate(free):
                trial[j] = solved[position]
            if all(0.0 < trial[j] < upper[j] for j in free):
                x = trial
                break
            steps = []
            for j in free:
                if trial[j] <= 0.0 and x[j] > trial[j]:
                    steps.append(x[j] / (x[j] - trial[j]))
                elif trial[j] >= upper[j] and trial[j] > x[j]:
                    steps.append((upper[j] - x[j]) / (trial[j] - x[j]))
            alpha = max(0.0, min(steps)) if steps else 0.0
            x = [x[j] + alpha * (trial[j] - x[j]) if j in free else x[j] for j in range(count)]
            still_free = []
            for j in free:
                if x[j] <= bound_tolerance(j):
                    x[j] = 0.0
                elif x[j] >= upper[j] - bound_tolerance(j):
                    x[j] = upper[j]
                    at_upper.add(j)
                else:
                    still_free.append(j)
            free = still_free
            if not free:
                break
        if x == previous:
            blocked.add(entering)
        else:
            blocked.clear()
        w = gradient(x)


@dataclass(frozen=True, slots=True)
class ThrustAllocation:
    """Non-negative thruster commands that best reproduce a force/torque demand.

    Parameters
    ----------
    thrusts_n : tuple of (str, float)
        Thrust per thruster identifier, in N; each in ``[0, maximum]``.
    force_n, torque_n_m : tuple of float
        Achieved force (N) and torque about the centre of mass (N m).
    force_residual_n, torque_residual_n_m : tuple of float
        Achieved minus demanded force and torque.
    weighted_residual_norm : float
        Residual norm with torque rows divided by the torque length scale.
    exact : bool
        Whether the weighted residual is within the requested tolerance,
        i.e. the demand is achievable within the limits.
    thrusters_at_maximum : tuple of str
        Thrusters at their upper bound.
    iterations : int
        Active-set iterations used.

    Notes
    -----
    When ``exact`` is false the allocation is the best bounded approximation
    [1]_ [2]_ and ``thrusters_at_maximum`` shows which limits bind, so the
    caller can scale the demand or re-plan.

    References
    ----------
    .. [1] C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems",
           Prentice-Hall, 1974 (NNLS algorithm; reprinted by SIAM, 1995).
    .. [2] P. B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: an
           Algorithm and Applications", Computational Statistics 10, 1995.
    """

    thrusts_n: tuple[tuple[str, float], ...]
    force_n: Vector
    torque_n_m: Vector
    force_residual_n: Vector
    torque_residual_n_m: Vector
    weighted_residual_norm: float
    exact: bool
    thrusters_at_maximum: tuple[str, ...]
    iterations: int


@dataclass(frozen=True, slots=True)
class ThrusterLayout:
    """Immutable set of body-mounted thrusters.

    Parameters
    ----------
    mounts : tuple of ThrusterMount
        At least one mount; identifiers must be unique.

    Raises
    ------
    DomainError
        If the layout is empty or identifiers are duplicated.

    Notes
    -----
    For thruster ``i`` at ``p_i`` with axis ``u_i`` and thrust ``T_i``, about
    the centre of mass ``c`` [1]_ [2]_::

        F = sum T_i u_i,     tau = sum (p_i - c) x T_i u_i
        d tau / d c = [F]x

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    .. [2] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", https://ntrs.nasa.gov/citations/19710016722.
    """

    mounts: tuple[ThrusterMount, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "mounts", records(self.mounts, ThrusterMount, "Thruster mounts"))
        if not self.mounts:
            raise DomainError("A thruster layout requires at least one thruster.")
        identifiers = [mount.thruster_id for mount in self.mounts]
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Thruster identifiers must be unique.")

    @property
    def thruster_ids(self) -> tuple[str, ...]:
        """Thruster identifiers in mount order."""
        return tuple(mount.thruster_id for mount in self.mounts)

    def influence_matrix(self, center_of_mass_m: Vector) -> tuple[tuple[float, ...], ...]:
        """Return the 6xN force/torque produced by one newton on each thruster.

        Parameters
        ----------
        center_of_mass_m : tuple of float
            Body-frame centre of mass, in m.

        Returns
        -------
        tuple of tuple of float
            Rows 0-2 are force (N/N), rows 3-5 torque about
            ``center_of_mass_m`` (N m/N); columns follow
            :attr:`thruster_ids`.  A thrust allocator inverts this matrix
            subject to non-negative, bounded thrust.

        Raises
        ------
        DomainError
            If the centre of mass is not a finite three-vector.

        References
        ----------
        .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
        """

        center = _vector(center_of_mass_m, "Centre of mass")
        columns = [
            (*mount.thrust_axis, *_cross(_sub(mount.position_m, center), mount.thrust_axis))
            for mount in self.mounts
        ]
        return tuple(tuple(column[row] for column in columns) for row in range(6))

    def resultant(
        self, thrusts_n: Mapping[str, float], center_of_mass_m: Vector
    ) -> ForceTorque:
        """Return net force, torque about the centre of mass and ``d tau/d c``.

        Parameters
        ----------
        thrusts_n : Mapping of str to float
            Thrust per thruster identifier, in N; finite and non-negative.
            Thrusters not listed are off.
        center_of_mass_m : tuple of float
            Body-frame centre of mass, in m.

        Returns
        -------
        ForceTorque
            Resultant force, torque and torque Jacobian.

        Raises
        ------
        DomainError
            If a thruster is unknown or a thrust is negative or non-finite.

        References
        ----------
        .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
        .. [2] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
               Maneuvers", https://ntrs.nasa.gov/citations/19710016722.
        """

        center = _vector(center_of_mass_m, "Centre of mass")
        thrusts_n = real_mapping(thrusts_n, "Thruster thrusts")
        unknown = set(thrusts_n) - set(self.thruster_ids)
        if unknown:
            raise DomainError("Unknown thrusters: " + ", ".join(sorted(unknown)))
        force = [0.0, 0.0, 0.0]
        torque = [0.0, 0.0, 0.0]
        total = 0.0
        for mount in self.mounts:
            thrust = thrusts_n.get(mount.thruster_id, 0.0)
            if thrust < 0.0:
                raise DomainError("Thruster thrust must be finite and non-negative.")
            if thrust == 0.0:
                continue
            total += thrust
            arm = _cross(_sub(mount.position_m, center), mount.thrust_axis)
            for index in range(3):
                force[index] += thrust * mount.thrust_axis[index]
                torque[index] += thrust * arm[index]
        net: Vector = (force[0], force[1], force[2])
        return ForceTorque(
            force_n=net,
            torque_n_m=(torque[0], torque[1], torque[2]),
            total_thrust_n=total,
            torque_wrt_center_of_mass=_skew(net),
        )

    def allocate(
        self,
        *,
        force_n: Vector,
        torque_n_m: Vector,
        center_of_mass_m: Vector,
        torque_length_scale_m: float = 1.0,
        maximum_thrust_n: float | Mapping[str, float] | None = None,
        tolerance: float = 1e-9,
        balance_weight: float = 0.0,
    ) -> ThrustAllocation:
        """Find bounded, non-negative thrusts minimizing the force/torque mismatch.

        Parameters
        ----------
        force_n : tuple of float
            Demanded body-frame force, in N.
        torque_n_m : tuple of float
            Demanded torque about the centre of mass, in N m.
        center_of_mass_m : tuple of float
            Body-frame centre of mass, in m.
        torque_length_scale_m : float, optional
            Positive length dividing torque rows so newtons and newton-metres
            are weighed consistently, in m.
        maximum_thrust_n : float or Mapping of str to float or None, optional
            Common or per-thruster upper bound, in N; ``None`` is unbounded.
            A mapping must name every thruster.  Bounds must be positive.
        tolerance : float, optional
            Positive relative tolerance on the weighted residual for
            ``exact``.
        balance_weight : float, optional
            Non-negative regularization weight ``lambda``.

        Returns
        -------
        ThrustAllocation
            Thrusts, achieved force and torque, residuals and binding limits.

        Raises
        ------
        DomainError
            If an input is outside its domain or a bound mapping is
            incomplete.
        ConvergenceError
            If the active-set iteration does not converge.

        Notes
        -----
        Solves ``min ||W (A T - b)||`` subject to ``0 <= T <= T_max`` with a
        bounded-variable active-set method (Lawson-Hanson NNLS [1]_ extended
        to upper bounds [2]_), where ``A`` is :meth:`influence_matrix` and
        ``W`` divides torque rows by ``torque_length_scale_m``.  If several
        allocations meet the demand exactly and ``balance_weight`` is zero,
        one of them is returned.  A positive ``balance_weight`` (lambda,
        dimensionless: weighted residual per newton of thrust) solves the
        strictly convex problem ``min ||W (A T - b)||^2 + lambda^2 ||T||^2``
        instead: the unique minimum-effort allocation, at the cost of a
        residual of order ``lambda^2``.  ``exact`` compares the weighted
        residual with ``tolerance`` times ``max(1, ||W b||)``.

        References
        ----------
        .. [1] C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems",
               Prentice-Hall, 1974 (NNLS algorithm; reprinted by SIAM, 1995).
        .. [2] P. B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: an
               Algorithm and Applications", Computational Statistics 10, 1995.
        """

        desired_force = _vector(force_n, "Desired force")
        desired_torque = _vector(torque_n_m, "Desired torque")
        length = real(torque_length_scale_m, "Torque length scale")
        if length <= 0.0:
            raise DomainError("Torque length scale must be greater than zero.")
        tolerance = real(tolerance, "Allocation tolerance")
        if tolerance <= 0.0:
            raise DomainError("Allocation tolerance must be greater than zero.")
        if maximum_thrust_n is None:
            maxima = {identifier: float("inf") for identifier in self.thruster_ids}
        elif isinstance(maximum_thrust_n, Mapping):
            maxima = real_mapping(maximum_thrust_n, "Maximum thrust")
            if set(maxima) != set(self.thruster_ids):
                raise DomainError("Maximum thrust must be given for every thruster.")
        else:
            limit = real(maximum_thrust_n, "Maximum thrust")
            maxima = {identifier: limit for identifier in self.thruster_ids}
        if any(value <= 0.0 for value in maxima.values()):
            raise DomainError("Maximum thrust must be greater than zero.")
        matrix = self.influence_matrix(center_of_mass_m)
        weights = (1.0, 1.0, 1.0, 1.0 / length, 1.0 / length, 1.0 / length)
        balance = real(balance_weight, "Balance weight")
        if balance < 0.0:
            raise DomainError("Balance weight must be non-negative.")
        count = len(self.mounts)
        columns = [
            [weights[row] * matrix[row][column] for row in range(6)]
            + ([balance if other == column else 0.0 for other in range(count)] if balance else [])
            for column in range(count)
        ]
        target = [weights[row] * value for row, value in enumerate((*desired_force, *desired_torque))]
        if balance:
            target += [0.0] * count
        demand_scale = max(1.0, sqrt(sum(value * value for value in target)))
        thrusts, iterations = _bvls(
            columns,
            target,
            [maxima[identifier] for identifier in self.thruster_ids],
            tolerance=1e-14 * demand_scale,
            max_iterations=60 * max(1, len(columns)),
        )
        allocation = dict(zip(self.thruster_ids, thrusts, strict=True))
        result = self.resultant(allocation, center_of_mass_m)
        force_residual = tuple(a - b for a, b in zip(result.force_n, desired_force, strict=True))
        torque_residual = tuple(
            a - b for a, b in zip(result.torque_n_m, desired_torque, strict=True)
        )
        norm = sqrt(
            sum(value * value for value in force_residual)
            + sum((value / length) ** 2 for value in torque_residual)
        )
        return ThrustAllocation(
            thrusts_n=tuple(allocation.items()),
            force_n=result.force_n,
            torque_n_m=result.torque_n_m,
            force_residual_n=force_residual,  # type: ignore[arg-type]
            torque_residual_n_m=torque_residual,  # type: ignore[arg-type]
            weighted_residual_norm=norm,
            exact=norm <= tolerance * demand_scale,
            thrusters_at_maximum=tuple(
                identifier
                for identifier, value in allocation.items()
                if value >= maxima[identifier] * (1.0 - 1e-12)
            ),
            iterations=iterations,
        )

    def axial_efficiency(self, thrusts_n: Mapping[str, float], axis: Vector) -> float:
        """Return the axial (cosine) efficiency ``(F . a) / (|a| sum T_i)``.

        Parameters
        ----------
        thrusts_n : Mapping of str to float
            Thrust per thruster identifier, in N.
        axis : tuple of float
            Non-zero reference axis; need not be a unit vector.

        Returns
        -------
        float
            Fraction of the summed thrust delivered along ``axis``; less than
            one for a canted cluster.

        Raises
        ------
        DomainError
            If the axis is zero, a thrust is invalid, or no thruster fires.

        Notes
        -----
        The ratio quantifies the cosine loss of canted thrusters and the axial
        component of a misaligned thrust vector [1]_ [2]_.

        References
        ----------
        .. [1] NASA-SP-125. NASA, "Design of Liquid Propellant Rocket Engines",
               2nd ed., https://ntrs.nasa.gov/citations/19710019929.
        .. [2] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
               Maneuvers", https://ntrs.nasa.gov/citations/19710016722.
        """

        direction = _vector(axis, "Axis")
        norm = sqrt(_dot(direction, direction))
        if norm <= 0.0:
            raise DomainError("Axis must be non-zero.")
        result = self.resultant(thrusts_n, _ZERO)
        if result.total_thrust_n <= 0.0:
            raise DomainError("Axial efficiency requires at least one firing thruster.")
        return _dot(result.force_n, direction) / (norm * result.total_thrust_n)


@dataclass(frozen=True, slots=True)
class PropellantStation:
    """Point-mass location of one propagated tank state in the body frame.

    Parameters
    ----------
    state_name : str
        Name of the propagated tank-mass state.
    position_m : tuple of float
        Body-frame propellant centroid, in m.

    Raises
    ------
    DomainError
        If the name is blank or the position is not a finite three-vector.

    Notes
    -----
    The propellant is lumped at a fixed centroid [1]_; slosh and centroid
    travel inside the tank are not represented.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    state_name: str
    position_m: Vector

    def __post_init__(self) -> None:
        text(self.state_name, "Propellant station state name")
        object.__setattr__(self, "position_m", _vector(self.position_m, "Station position"))


@dataclass(frozen=True, slots=True)
class MassPropertiesEvaluation:
    """Mass, centre of mass, inertia and their tank-mass derivatives.

    Parameters
    ----------
    mass_kg : float
        Total vehicle mass, in kg.
    center_of_mass_m : tuple of float
        Body-frame centre of mass, in m.
    inertia_kg_m2 : tuple of tuple of float
        Inertia about the current centre of mass in body axes, kg m**2.
    center_of_mass_wrt_state : tuple of (str, tuple of float)
        ``dc/dm_k`` per tank state, in m/kg.
    inertia_wrt_state : tuple of (str, 3x3 tuple)
        ``dI_c/dm_k`` per tank state, in m**2.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    mass_kg: float
    center_of_mass_m: Vector
    inertia_kg_m2: Matrix3
    center_of_mass_wrt_state: tuple[tuple[str, Vector], ...]
    inertia_wrt_state: tuple[tuple[str, Matrix3], ...]


@dataclass(frozen=True, slots=True)
class PropellantMassModel:
    """Dry rigid body plus point-mass propellant stations.

    Parameters
    ----------
    dry_mass_kg : float
        Dry mass ``m_d``, in kg.  Finite and positive.
    dry_center_of_mass_m : tuple of float
        Dry centre of mass ``d``, in m.
    dry_inertia_kg_m2 : 3x3 tuple of float
        Dry inertia ``J_d`` about ``d``, in kg m**2; symmetric with
        non-negative diagonal satisfying the triangle inequality.
    stations : tuple of PropellantStation
        One station per tank state.

    Raises
    ------
    DomainError
        If the dry mass is non-positive, the inertia is not a symmetric 3x3
        matrix with admissible principal values, or a tank state has more
        than one station.

    Notes
    -----
    With tank masses ``m_k`` at ``p_k``, ``M = m_d + sum m_k`` and
    ``c = (m_d d + sum m_k p_k) / M``.  Using the parallel-axis operator
    ``S(r) = (r.r) I - r r^T`` [1]_::

        I_c = J_d + m_d S(d) + sum m_k S(p_k) - M S(c)
        dc/dm_k   = (p_k - c) / M
        dI_c/dm_k = S(p_k) - S(c) - [2 (c.e) I - (e c^T + c e^T)],  e = p_k - c

    The symmetry check uses a relative tolerance of ``1e-9`` and the
    triangle-inequality check an absolute tolerance of ``1e-9`` kg m**2.

    References
    ----------
    .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
    """

    dry_mass_kg: float
    dry_center_of_mass_m: Vector
    dry_inertia_kg_m2: Matrix3
    stations: tuple[PropellantStation, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "dry_mass_kg", real(self.dry_mass_kg, "Dry mass"))
        if self.dry_mass_kg <= 0.0:
            raise DomainError("Dry mass must be finite and greater than zero.")
        object.__setattr__(
            self, "dry_center_of_mass_m", _vector(self.dry_center_of_mass_m, "Dry centre of mass")
        )
        inertia = real_matrix(self.dry_inertia_kg_m2, "Dry inertia")
        if len(inertia) != 3 or any(len(row) != 3 for row in inertia):
            raise DomainError("Dry inertia must be a 3x3 matrix.")
        object.__setattr__(self, "dry_inertia_kg_m2", inertia)
        for i in range(3):
            for j in range(i):
                if abs(inertia[i][j] - inertia[j][i]) > 1e-9 * max(1.0, abs(inertia[i][j])):
                    raise DomainError("Dry inertia must be symmetric.")
        principal = (inertia[0][0], inertia[1][1], inertia[2][2])
        if any(value < 0.0 for value in principal) or any(
            principal[i] > principal[(i + 1) % 3] + principal[(i + 2) % 3] + 1e-9
            for i in range(3)
        ):
            raise DomainError("Dry inertia diagonal violates the triangle inequality.")
        object.__setattr__(
            self, "stations", records(self.stations, PropellantStation, "Propellant stations")
        )
        names = [station.state_name for station in self.stations]
        if len(names) != len(set(names)):
            raise DomainError("Each tank state may have only one propellant station.")

    def evaluate(self, additional_states: Mapping[str, float]) -> MassPropertiesEvaluation:
        """Return mass properties for the current tank-state masses.

        Parameters
        ----------
        additional_states : Mapping of str to float
            Tank-state masses, in kg; must contain every station's state and
            be finite and non-negative.

        Returns
        -------
        MassPropertiesEvaluation
            Mass, centre of mass, inertia about it, and their analytic
            derivatives with respect to every tank mass (see class Notes).

        Raises
        ------
        DomainError
            If a station's state is missing or a mass is negative or
            non-finite.

        References
        ----------
        .. [1] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
        """

        states = real_mapping(additional_states, "Tank states")
        masses: list[float] = []
        for station in self.stations:
            if station.state_name not in states:
                raise DomainError(f"Missing tank state {station.state_name!r}.")
            value = states[station.state_name]
            if value < 0.0:
                raise DomainError("Tank-state masses must be finite and non-negative.")
            masses.append(value)
        total = self.dry_mass_kg + sum(masses)
        moment = [self.dry_mass_kg * value for value in self.dry_center_of_mass_m]
        for mass, station in zip(masses, self.stations, strict=True):
            for index in range(3):
                moment[index] += mass * station.position_m[index]
        center: Vector = (moment[0] / total, moment[1] / total, moment[2] / total)
        terms: list[tuple[float, Matrix3]] = [
            (1.0, self.dry_inertia_kg_m2),
            (self.dry_mass_kg, _shift_operator(self.dry_center_of_mass_m)),
            (-total, _shift_operator(center)),
        ]
        terms.extend(
            (mass, _shift_operator(station.position_m))
            for mass, station in zip(masses, self.stations, strict=True)
        )
        inertia = _matrix_add(*terms)
        center_partials: list[tuple[str, Vector]] = []
        inertia_partials: list[tuple[str, Matrix3]] = []
        shift_center = _shift_operator(center)
        for station in self.stations:
            offset = _sub(station.position_m, center)
            center_partials.append(
                (station.state_name, tuple(value / total for value in offset))  # type: ignore[arg-type]
            )
            correction = tuple(
                tuple(
                    (2.0 * _dot(center, offset) if i == j else 0.0)
                    - (offset[i] * center[j] + center[i] * offset[j])
                    for j in range(3)
                )
                for i in range(3)
            )
            inertia_partials.append(
                (
                    station.state_name,
                    _matrix_add(
                        (1.0, _shift_operator(station.position_m)),
                        (-1.0, shift_center),
                        (-1.0, correction),  # type: ignore[arg-type]
                    ),
                )
            )
        return MassPropertiesEvaluation(
            mass_kg=total,
            center_of_mass_m=center,
            inertia_kg_m2=inertia,
            center_of_mass_wrt_state=tuple(center_partials),
            inertia_wrt_state=tuple(inertia_partials),
        )
