# Thruster geometry and mass properties

`propulsion/burns/vehicle_geometry.py` provides the body-frame quantities that
attitude-coupled (6-DOF) simulation, thrust-vector alignment checks and
momentum-unloading analysis need. Attitude propagation and frame rotation stay
with the consumer.

## Thruster layout

A `ThrusterMount` has a body-frame application point \(\mathbf p_i\) and the
unit direction \(\hat u_i\) of the force on the vehicle (opposite to the
exhaust). For thrusts \(T_i\) and centre of mass \(\mathbf c\):

\[
  \mathbf F=\sum_i T_i\hat u_i,\qquad
  \boldsymbol\tau=\sum_i(\mathbf p_i-\mathbf c)\times T_i\hat u_i,\qquad
  \frac{\partial\boldsymbol\tau}{\partial\mathbf c}=[\mathbf F]_\times .
\]

`ThrusterLayout` provides:

- `resultant(thrusts, center_of_mass)` → force, torque, total thrust and
  \(\partial\boldsymbol\tau/\partial\mathbf c\);
- `influence_matrix(center_of_mass)` → the 6×N force/torque per newton, the
  matrix a thrust allocator inverts under non-negative bounded thrust (the
  allocator itself is out of scope);
- `axial_efficiency(thrusts, axis)` → \((\mathbf F\cdot\hat a)/\sum T_i\), the
  cosine loss of a canted cluster. This is the value to fold into the axial
  thrust of a tabulated artifact; the bridge does not apply cant twice.

### Thrust allocation

`allocate(force_n=..., torque_n_m=..., center_of_mass_m=..., maximum_thrust_n=...)`
returns the thrusts that best reproduce a force/torque demand within the
thruster limits:

\[
  \min_{0\le T\le T_\max}\ \lVert W(AT-b)\rVert,\qquad
  W=\operatorname{diag}(1,1,1,1/L,1/L,1/L),
\]

where \(A\) is the influence matrix and \(L\) (`torque_length_scale_m`)
puts newtons and newton-metres on a common scale. It is solved with a
bounded-variable active-set method: Lawson–Hanson NNLS extended to upper
bounds in the manner of Stark and Parker. It is dependency-free, with
Householder QR for the free-set least-squares subproblems. Thrusters can only
push, so non-negativity is a hard constraint, and the upper limits are imposed
inside the optimization.

Imposing the limits inside the optimization is important. A sparse
non-negative solution can need more than the limit on two thrusters even when
spreading the load over four meets the demand exactly; clipping such a
solution afterwards would silently change the force. The bounded solver finds
the feasible split. When the demand is not achievable within the limits, the
result is the best bounded approximation: `exact` is false and
`thrusters_at_maximum` shows which limits bind.

When several allocations meet the demand exactly, one of them is returned. A
positive `balance_weight` \(\lambda\) instead solves the strictly convex
problem

\[
  \min_{0\le T\le T_\max}\ \lVert W(AT-b)\rVert^2+\lambda^2\lVert T\rVert^2 ,
\]

the unique minimum-effort allocation. It spreads thrust over equivalent
thrusters, at the cost of a residual of order \(\lambda^2\). It is solved
exactly by the same bounded solver with \(\lambda I\) appended to the
matrix. The tests check the closed-form symmetric split
\(T_i = 60/(4+\lambda^2)\) and certify optimality on random bounded,
weighted problems through the KKT conditions. The tests compare results with the
exact optimum found by enumerating every lower/upper/free assignment of the
thrusters, on random layouts with random limits.

Momentum unloading with unbalanced thrusters produces a small net
\(\mathbf F\); its \(\Delta\mathbf v=\mathbf F\,\Delta t/m\) is a classic
orbit-determination error source and follows directly from `resultant`.

## Propellant-dependent mass properties

`PropellantMassModel` combines a dry rigid body (mass \(m_d\), centre of mass
\(\mathbf d\), inertia \(J_d\) about \(\mathbf d\)) with point-mass propellant
stations whose masses \(m_k\) are the propagated tank states. With
\(S(\mathbf r)=(\mathbf r\cdot\mathbf r)I-\mathbf r\mathbf r^{\mathsf T}\):

\[
  M=m_d+\sum m_k,\qquad
  \mathbf c=\frac{m_d\mathbf d+\sum m_k\mathbf p_k}{M},\qquad
  I_c=J_d+m_dS(\mathbf d)+\sum m_kS(\mathbf p_k)-MS(\mathbf c).
\]

Analytic derivatives with respect to each tank mass:

\[
  \frac{\partial\mathbf c}{\partial m_k}=\frac{\mathbf p_k-\mathbf c}{M},\qquad
  \frac{\partial I_c}{\partial m_k}=S(\mathbf p_k)-S(\mathbf c)
  -\big[2(\mathbf c\cdot\mathbf e)I-(\mathbf e\mathbf c^{\mathsf T}+\mathbf c\mathbf e^{\mathsf T})\big],
  \quad \mathbf e=\mathbf p_k-\mathbf c .
\]

Combined with the bridge's tank drains \(\dot m_k\), these give
\(\dot{\mathbf c}\) and \(\dot I_c\) for the rotational equations, and with
`ThrusterLayout` the torque change as the centre of mass migrates. Both
partials are checked against central finite differences, and the inertia is
checked against an independent parallel-axis sum.

## Limits

- Propellant is a point mass at a fixed centroid: no slosh, no centroid travel
  as a tank empties, no bladder/diaphragm mechanics.
- Minimum impulse bits and pulse quantization are not modelled.
- Gimbal dynamics and plume impingement are not modelled.

## References

- C. L. Lawson and R. J. Hanson, *Solving Least Squares Problems*,
  Prentice-Hall, 1974; reprinted by SIAM, 1995 (NNLS algorithm; cited
  bibliographically).
- P. B. Stark and R. L. Parker, *Bounded-Variable Least-Squares: an Algorithm
  and Applications*, Computational Statistics 10, 1995 (cited
  bibliographically).

- B. Wie, *Space Vehicle Dynamics and Control*, 2nd ed., AIAA, 2008
  (cited bibliographically).
- [NASA SP-125, Design of Liquid Propellant Rocket Engines](https://ntrs.nasa.gov/citations/19710019929)
  (the source this repository already cites for engine-cluster cant
  efficiency in `cluster.py`).
