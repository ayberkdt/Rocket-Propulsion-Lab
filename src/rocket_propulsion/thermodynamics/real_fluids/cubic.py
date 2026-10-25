"""Peng-Robinson and Soave-Redlich-Kwong cubic equations as a fast fallback.

Both equations share the form [1]_ (Eq. 9)

    p = R T / (v - b) - a(T) / ((v + Delta1 b) (v + Delta2 b)),

with ``a(T) = a0 (1 + m (1 - sqrt(T/Tc)))^2`` [1]_ (Eq. 13). The constants
are those of Table 1 of Bell and Jaeger [1]_, which cites Peng and Robinson
[2]_ (and Robinson and Peng [4]_ for ``omega > 0.491``) and Soave [3]_:

=====  =================  ===============  =========================================
EOS    a0                 b                m(omega)
=====  =================  ===============  =========================================
SRK    0.42747 R2Tc2/pc   0.08664 RTc/pc   0.48 + 1.574 w - 0.176 w2
PR     0.45724 R2Tc2/pc   0.07780 RTc/pc   0.37464 + 1.54226 w - 0.26992 w2
                                           (w > 0.491: 0.379642 + 1.48503 w
                                           - 0.164423 w2 + 0.016666 w3)
=====  =================  ===============  =========================================

``(Delta1, Delta2)`` is ``(1, 0)`` for SRK and ``(1 + sqrt 2, 1 - sqrt 2)``
for PR. The pressure-explicit equation is converted to the residual reduced
Helmholtz energy [1]_ (Eqs. 25-29)

    alphar = -ln(1 - b rho_r delta)
             - a(T)/(R T) ln((1 + Delta1 b rho_r delta)/(1 + Delta2 b rho_r delta))
               / (b (Delta1 - Delta2)),

so the cubic models plug into the same :class:`RealFluid` machinery as the
reference equations; departure functions follow from
:meth:`RealFluid.departures`. Only ``T_c``, ``p_c`` and the acentric factor
are needed, and the accuracy is correspondingly lower: every result carries a
``cubic_eos_lower_accuracy`` warning. The original papers [2]_ [3]_ were not
opened for this work; the equations and constants are those printed by [1]_.

References
----------
.. [1] BELL-JAGER-2016-CUBIC. I. H. Bell and A. Jaeger, "Helmholtz Energy
       Transformations of Common Cubic Equations of State for Use with Pure
       Fluids and Mixtures", J. Res. NIST 121, 238-263, 2016, Eq. (9) p. 241,
       Table 1 and Eq. (13) p. 242, Eqs. (25)-(29) p. 244,
       https://doi.org/10.6028/jres.121.011
.. [2] PENG-ROBINSON-1976. D.-Y. Peng and D. B. Robinson, "A New Two-Constant
       Equation of State", Ind. Eng. Chem. Fundam. 15(1), 59-64, 1976,
       https://doi.org/10.1021/i160057a011
.. [3] SOAVE-1972. G. Soave, "Equilibrium Constants from a Modified
       Redlich-Kwong Equation of State", Chem. Eng. Sci. 27(6), 1197-1203,
       1972, https://doi.org/10.1016/0009-2509(72)80096-4
.. [4] ROBINSON-PENG-1978. D. B. Robinson and D.-Y. Peng, "The
       Characterization of the Heptanes and Heavier Fractions for the GPA
       Peng-Robinson Programs", GPA Research Report RR-28, 1978, as cited
       in [1]_ (reference 23), https://doi.org/10.6028/jres.121.011
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import log, sqrt
from typing import Literal

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import real

from .helmholtz import FluidModel, IdealHelmholtz, RealFluid, ResidualDerivatives
from .reference_fluids import reference_fluid_info, reference_fluid_model

CubicKind = Literal["peng_robinson", "soave_redlich_kwong"]
_SQRT2 = sqrt(2.0)
_CONSTANTS = {
    # kind: (Omega_a, Omega_b, Delta1, Delta2) -- Bell and Jaeger (2016) Table 1
    "peng_robinson": (0.45724, 0.07780, 1.0 + _SQRT2, 1.0 - _SQRT2),
    "soave_redlich_kwong": (0.42747, 0.08664, 1.0, 0.0),
}


def _m_parameter(kind: str, omega: float) -> float:
    if kind == "soave_redlich_kwong":
        return 0.48 + 1.574 * omega - 0.176 * omega * omega
    if omega <= 0.491:
        return 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    return 0.379642 + 1.48503 * omega - 0.164423 * omega ** 2 + 0.016666 * omega ** 3


@dataclass(frozen=True, slots=True)
class CubicResidual:
    """Residual Helmholtz energy of a Peng-Robinson or SRK equation.

    Parameters
    ----------
    kind : {"peng_robinson", "soave_redlich_kwong"}
        Equation.
    t_crit_k, p_crit_pa : float
        Critical temperature and pressure.
    acentric_factor : float
        Pitzer acentric factor.
    gas_constant_j_mol_k : float
        Molar gas constant.
    t_reduce_k, rho_reduce_mol_m3 : float
        Reducing parameters of ``tau`` and ``delta``; ``t_reduce_k`` must
        equal ``t_crit_k`` [1]_ (Sec. 2.2.1).

    Attributes
    ----------
    a0_pa_m6_mol2, b_m3_mol, m : float
        Attraction constant, covolume and alpha-function slope.
    reduced_density_limit : float
        ``1 / (b rho_r)``, the covolume limit.

    Raises
    ------
    InputError
        If the kind is unknown or a constant is not positive.

    References
    ----------
    .. [1] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238
           (2016), Table 1, Eqs. (13) and (25)-(29),
           https://doi.org/10.6028/jres.121.011
    """

    kind: CubicKind
    t_crit_k: float
    p_crit_pa: float
    acentric_factor: float
    gas_constant_j_mol_k: float
    t_reduce_k: float
    rho_reduce_mol_m3: float
    a0_pa_m6_mol2: float = field(default=0.0, init=False)
    b_m3_mol: float = field(default=0.0, init=False)
    m: float = field(default=0.0, init=False)
    reduced_density_limit: float = field(default=0.0, init=False)

    def __post_init__(self) -> None:
        if self.kind not in _CONSTANTS:
            raise InputError(f"Unknown cubic equation {self.kind!r}.", field="kind")
        for name in ("t_crit_k", "p_crit_pa", "gas_constant_j_mol_k", "t_reduce_k",
                     "rho_reduce_mol_m3"):
            value = real(getattr(self, name), name)
            if value <= 0.0:
                raise InputError(f"{name} must be positive.", field=name)
            object.__setattr__(self, name, value)
        if self.t_reduce_k != self.t_crit_k:
            raise InputError("t_reduce_k must equal t_crit_k for a pure cubic fluid.",
                             field="t_reduce_k")
        omega = real(self.acentric_factor, "acentric_factor")
        object.__setattr__(self, "acentric_factor", omega)
        omega_a, omega_b, _, _ = _CONSTANTS[self.kind]
        r, tc, pc = self.gas_constant_j_mol_k, self.t_crit_k, self.p_crit_pa
        object.__setattr__(self, "a0_pa_m6_mol2", omega_a * r * r * tc * tc / pc)
        object.__setattr__(self, "b_m3_mol", omega_b * r * tc / pc)
        object.__setattr__(self, "m", _m_parameter(self.kind, omega))
        object.__setattr__(self, "reduced_density_limit",
                           1.0 / (self.b_m3_mol * self.rho_reduce_mol_m3))

    def attraction(self, temperature_k: float) -> tuple[float, float, float]:
        """Return ``a(T)``, ``da/dT`` and ``d2a/dT2``.

        Parameters
        ----------
        temperature_k : float
            Positive temperature.

        Returns
        -------
        tuple of float
            Attraction parameter (Pa m6/mol2) and its derivatives.

        Raises
        ------
        DomainError
            If ``1 + m (1 - sqrt(T/Tc)) <= 0``, where the alpha function
            stops decreasing and is no longer physical.

        References
        ----------
        .. [1] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238
               (2016), Eq. (13), p. 242, https://doi.org/10.6028/jres.121.011
        """

        root = sqrt(temperature_k * self.t_crit_k)
        f = 1.0 + self.m * (1.0 - sqrt(temperature_k / self.t_crit_k))
        if f <= 0.0:
            raise DomainError(
                f"Cubic alpha function is non-physical at T = {temperature_k} K "
                f"(1 + m (1 - sqrt(T/Tc)) = {f}).", field="temperature_k")
        f1 = -self.m / (2.0 * root)
        f2 = self.m / (4.0 * temperature_k * root)
        a0 = self.a0_pa_m6_mol2
        return a0 * f * f, 2.0 * a0 * f * f1, 2.0 * a0 * (f1 * f1 + f * f2)

    def derivatives(self, tau: float, delta: float) -> ResidualDerivatives:
        """Return ``alphar`` and its first and second derivatives.

        Parameters
        ----------
        tau, delta : float
            Inverse reduced temperature and reduced density.

        Returns
        -------
        tuple of float
            ``(alphar, alphar_delta, alphar_delta_delta, alphar_tau,
            alphar_tau_tau, alphar_delta_tau)``.

        Raises
        ------
        DomainError
            At or beyond the covolume limit, or where the alpha function is
            non-physical.

        Notes
        -----
        With ``x = b rho_r delta``, ``A(tau) = a(T)/(R T)`` and
        ``psi = ln((1 + Delta1 x)/(1 + Delta2 x)) / (b (Delta1 - Delta2))``,
        ``alphar = -ln(1 - x) - A psi`` [1]_ (Eqs. 25-29);
        ``dA/dtau = (a - T a')/(R T_r)`` and
        ``d2A/dtau2 = T^3 a''/(R T_r^2)``.

        References
        ----------
        .. [1] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238
               (2016), Eqs. (25)-(29), p. 244,
               https://doi.org/10.6028/jres.121.011
        """

        _, _, d1, d2 = _CONSTANTS[self.kind]
        b, rr = self.b_m3_mol, self.rho_reduce_mol_m3
        x = b * rr * delta
        if x >= 1.0:
            raise DomainError("Density at or beyond the cubic covolume limit.",
                              field="density_kg_m3")
        temperature = self.t_reduce_k / tau
        a, a1, a2 = self.attraction(temperature)
        r = self.gas_constant_j_mol_k
        big_a = a / (r * temperature)
        big_a_t = (a - temperature * a1) / (r * self.t_reduce_k)
        big_a_tt = temperature ** 3 * a2 / (r * self.t_reduce_k ** 2)
        e1, e2 = 1.0 + d1 * x, 1.0 + d2 * x
        psi = log(e1 / e2) / (b * (d1 - d2))
        psi_d = rr / (e1 * e2)
        psi_dd = -rr * b * rr * (d1 * e2 + d2 * e1) / (e1 * e2) ** 2
        one = 1.0 - x
        m_d = b * rr / one
        return (-log(one) - big_a * psi, m_d - big_a * psi_d,
                m_d * m_d - big_a * psi_dd, -big_a_t * psi, -big_a_tt * psi,
                -big_a_t * psi_d)


def constant_cp_ideal(cp_over_r: float) -> IdealHelmholtz:
    """Return an ideal-gas part with constant ``c_p0 / R``.

    Parameters
    ----------
    cp_over_r : float
        Ideal-gas isobaric heat capacity divided by ``R``; greater than 1.

    Returns
    -------
    IdealHelmholtz
        ``alpha0 = ln(delta) + (c_p0/R - 1) ln(tau)``; the enthalpy and
        entropy reference is arbitrary.

    Raises
    ------
    InputError
        If ``cp_over_r <= 1``.

    References
    ----------
    .. [1] ORTIZVEGA-2023-HE. Ortiz Vega et al., NIST IR 8474 (2023), Eq. (3),
           p. 3 (the monatomic case c_p0/R = 2.5),
           https://doi.org/10.6028/NIST.IR.8474
    """

    value = real(cp_over_r, "cp_over_r")
    if value <= 1.0:
        raise InputError("cp_over_r must exceed 1.", field="cp_over_r")
    return IdealHelmholtz(a1=0.0, a2=0.0, log_tau=value - 1.0)


def cubic_model(kind: CubicKind, *, name: str, t_crit_k: float, p_crit_pa: float,
                acentric_factor: float, molar_mass_kg_mol: float, ideal: IdealHelmholtz,
                t_min_k: float, t_max_k: float, p_max_pa: float,
                gas_constant_j_mol_k: float = MOLAR_GAS_CONSTANT_J_MOL_K,
                rho_reduce_mol_m3: float | None = None,
                references: tuple[str, ...] = ()) -> FluidModel:
    """Build a cubic :class:`FluidModel` for an arbitrary fluid.

    Parameters
    ----------
    kind : {"peng_robinson", "soave_redlich_kwong"}
        Equation.
    name : str
        Fluid name.
    t_crit_k, p_crit_pa, acentric_factor : float
        Sourced critical constants and acentric factor of the fluid.
    molar_mass_kg_mol : float
        Molar mass.
    ideal : IdealHelmholtz
        Ideal-gas part on the same reducing parameters (for example
        :func:`constant_cp_ideal`).
    t_min_k, t_max_k, p_max_pa : float
        Range the caller accepts; there is no intrinsic validity range for a
        cubic equation, so it must be stated explicitly.
    gas_constant_j_mol_k : float, optional
        Molar gas constant; default is the SI 2019 value.
    rho_reduce_mol_m3 : float, optional
        Reducing density; default ``p_c / (R T_c)``.
    references : tuple of str, optional
        Sources of the constants.

    Returns
    -------
    FluidModel
        Model with ``model_kind="cubic"``.

    Raises
    ------
    InputError
        If a parameter is invalid.

    References
    ----------
    .. [1] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238
           (2016), Table 1, https://doi.org/10.6028/jres.121.011
    """

    if not isinstance(ideal, IdealHelmholtz):
        raise InputError("ideal must be an IdealHelmholtz.", field="ideal")
    r = real(gas_constant_j_mol_k, "gas_constant_j_mol_k")
    tc = real(t_crit_k, "t_crit_k")
    pc = real(p_crit_pa, "p_crit_pa")
    if tc <= 0.0 or pc <= 0.0 or r <= 0.0:
        raise InputError("Critical constants and R must be positive.", field="t_crit_k")
    rho_r = pc / (r * tc) if rho_reduce_mol_m3 is None else real(rho_reduce_mol_m3,
                                                                  "rho_reduce_mol_m3")
    residual = CubicResidual(kind, tc, pc, acentric_factor, r, tc, rho_r)
    label = "Peng-Robinson" if kind == "peng_robinson" else "Soave-Redlich-Kwong"
    return FluidModel(
        name=f"{name} ({label} cubic)", key=f"{name}:{kind}", model_kind="cubic",
        molar_mass_kg_mol=molar_mass_kg_mol, gas_constant_j_mol_k=r, t_reduce_k=tc,
        rho_reduce_mol_m3=rho_r, t_crit_k=tc, p_crit_pa=pc, t_min_k=t_min_k,
        t_max_k=t_max_k, p_max_pa=p_max_pa, ideal=ideal, residual=residual,
        references=("BELL-JAGER-2016-CUBIC Table 1", *references),
        notes=("Cubic equation of state: lower accuracy fallback.",))


def cubic_counterpart(reference_key: str, kind: CubicKind = "peng_robinson") -> RealFluid:
    """Return a cubic calculator built from a reference fluid's constants.

    Parameters
    ----------
    reference_key : str
        Key of a reference fluid (see ``REFERENCE_FLUID_KEYS``).
    kind : {"peng_robinson", "soave_redlich_kwong"}, optional
        Equation; default Peng-Robinson.

    Returns
    -------
    RealFluid
        Calculator using the reference equation's ``T_c``, ``p_c``, gas
        constant, molar mass, ideal-gas part and validity range, the
        WebBook acentric factor, and ``rho_r`` equal to the reference
        critical density so that the reference ideal-gas part applies
        unchanged [1]_ (Sec. 2.2.1).

    Raises
    ------
    InputError
        If the key or kind is unknown.

    References
    ----------
    .. [1] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238
           (2016), Sec. 2.2.1, p. 243, https://doi.org/10.6028/jres.121.011
    .. [2] NIST-WEBBOOK-FLUID. NIST Chemistry WebBook SRD 69 (acentric
           factors), https://webbook.nist.gov/chemistry/fluid/
    """

    ref = reference_fluid_model(reference_key)
    info = reference_fluid_info(reference_key)
    model = cubic_model(
        kind, name=ref.name, t_crit_k=ref.t_crit_k, p_crit_pa=ref.p_crit_pa,
        acentric_factor=info.acentric_factor, molar_mass_kg_mol=ref.molar_mass_kg_mol,
        ideal=ref.ideal, t_min_k=ref.t_min_k, t_max_k=ref.t_max_k, p_max_pa=ref.p_max_pa,
        gas_constant_j_mol_k=ref.gas_constant_j_mol_k,
        rho_reduce_mol_m3=ref.rho_reduce_mol_m3,
        references=(*ref.references, "NIST-WEBBOOK-FLUID acentric factor"))
    return RealFluid(model)
