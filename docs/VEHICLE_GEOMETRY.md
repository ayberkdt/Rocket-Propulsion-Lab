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
- Thrust allocation, gimbal dynamics and plume impingement are not modelled.

## References

- B. Wie, *Space Vehicle Dynamics and Control*, 2nd ed., AIAA, 2008
  (cited bibliographically).
- [NASA SP-125, Design of Liquid Propellant Rocket Engines](https://ntrs.nasa.gov/citations/19710019929)
  (the source this repository already cites for engine-cluster cant
  efficiency in `cluster.py`).
