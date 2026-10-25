"""Reference viscosity and thermal conductivity of cryogenic propellants.

The correlations are the ones the NIST Chemistry WebBook (SRD 69) uses for
these fluids [6]_:

==================  ==================================  ==================================
key                 viscosity                           thermal conductivity
==================  ==================================  ==================================
``nitrogen``        Lemmon and Jacobsen 2004 [1]_       Lemmon and Jacobsen 2004 [1]_
``oxygen``          Lemmon and Jacobsen 2004 [1]_       Lemmon and Jacobsen 2004 [1]_
``parahydrogen``    Muzny et al. 2013 via CoolProp [7]_ Assael et al. 2011 [2]_
``normal_hydrogen`` Muzny et al. 2013 via CoolProp [7]_ Assael et al. 2011 [2]_
``methane``         Friend, Ely and Ingham 1989 [3]_    Friend, Ely and Ingham 1989 [3]_
``helium``          Arp et al. 1998 (TN 1334) [5]_      Hands and Arp 1981 via CoolProp [7]_
==================  ==================================  ==================================

Every correlation is a function of temperature and density, written as a
dilute-gas term, a residual (excess) term and, for the thermal conductivity,
a critical enhancement. Nitrogen, oxygen and hydrogen use the simplified
Olchowy-Sengers crossover model as printed in [1]_ (Eqs. 7-11) and [2]_
(Eqs. 4-7); methane uses the enhancement of [3]_ (Eqs. 18-22). The
thermodynamic inputs of the enhancements (``cp``, ``cv``, ``(drho/dp)_T``,
``(dp/dT)_rho``) and the density for ``(T, p)`` inputs come from the
reference equations of state of :mod:`.reference_fluids`, which are the
equations the transport correlations were fitted with.

Provenance. Coefficients of nitrogen, oxygen, hydrogen conductivity,
methane and helium viscosity were read from the primary publications in
this work. Two correlations were NOT read from their primary paper:

* Hydrogen viscosity (Muzny, Huber and Kazakov, J. Chem. Eng. Data 58,
  969-979, 2013, normal hydrogen; the WebBook uses it for parahydrogen too).
  The paper was not openly accessible and was not opened; only its erratum
  [4]_ was. The coefficients are transcribed from CoolProp [7]_
  (``Hydrogen.json``/``ParaHydrogen.json`` and
  ``viscosity_hydrogen_higher_order_hardcoded``), with the erratum's
  Avogadro constant 6.022137e23 /mol, exponent ``-i`` and scale density
  90.909090909 kg/m3. Correctness rests on that transcription, verified
  against NIST WebBook viscosities (16 rows per hydrogen form, largest
  relative deviation 5e-12) and the erratum's three check values.
* Helium conductivity (Hands and Arp, Cryogenics 21, 697-703, 1981, a
  subscription article that was not opened), transcribed from CoolProp's
  ``conductivity_hardcoded_helium``, which follows the REFPROP source
  (including two factors that CoolProp notes are not in the paper).
  Correctness rests on the transcription verified against NIST WebBook
  conductivities (largest relative deviation 1e-6).

Deviations from the WebBook, with reasons:

* Methane viscosity as used by the WebBook (Quinones-Cisneros, Huber and
  Deiters, unpublished, 2011 [6]_). The published Friend-Ely-Ingham
  viscosity is implemented instead; it differs from the WebBook values by
  up to about 4.5 % (see ``docs/REAL_FLUID_TRANSPORT.md``). Because the
  methane conductivity contains the viscosity (Eqs. 13a and 18a of [3]_),
  it differs from the WebBook by up to about 0.7 % (2 % near the critical
  point); it reproduces the paper's own Table A3.
* Eq. (16) of [3]_ prints ``delta*_sigma = 11`` outside the saturated-vapor
  branch; 1 is used because only 1 reproduces Table A3 of the same paper.

Example
-------
>>> from rocket_propulsion.thermodynamics.real_fluids.transport import transport_tp
>>> lox = transport_tp("oxygen", 90.0, 1.0e6)
>>> round(lox.viscosity_pa_s * 1e6, 2)  # NIST WebBook (accessed 2026-10-25): 197.43 uPa s
197.43

References
----------
.. [1] LEMMON-JACOBSEN-2004-TRANSPORT. E. W. Lemmon and R. T Jacobsen,
       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen,
       Argon, and Air", Int. J. Thermophys. 25(1), 21-69, 2004, Eqs. (1)-(11)
       and Tables I-V, pp. 23-28, https://trc.nist.gov/refprop/FAQ/NAO.PDF
.. [2] ASSAEL-2011-H2-CONDUCTIVITY. M. J. Assael, J.-A. M. Assael, M. L. Huber,
       R. A. Perkins and Y. Takata, "Correlation of the Thermal Conductivity
       of Normal and Parahydrogen from the Triple Point to 1000 K and up to
       100 MPa", J. Phys. Chem. Ref. Data 40, 033101, 2011, Eqs. (1)-(7),
       Tables 2, 5 and 7, https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
.. [3] FRIEND-ELY-INGHAM-1989-CH4. D. G. Friend, J. F. Ely and H. Ingham,
       "Thermophysical Properties of Methane", J. Phys. Chem. Ref. Data
       18(2), 583-638, 1989, Eqs. (9)-(22) and Tables 1, 8-10, pp. 590,
       596-598, https://srd.nist.gov/jpcrdreprint/1.555828.pdf
.. [4] MUZNY-2022-H2-VISCOSITY-ERRATUM. C. D. Muzny, M. L. Huber and
       A. F. Kazakov, "Erratum: Correlation for the Viscosity of Normal
       Hydrogen Obtained from Symbolic Regression",
       https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934278
.. [5] ARP-MCCARTY-FRIEND-1998-HE. V. D. Arp, R. D. McCarty and D. G. Friend,
       NIST Technical Note 1334 (revised), 1998, Sec. 4.2, Eqs. (1)-(9),
       pp. 4-5, https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
.. [6] NIST-WEBBOOK-FLUID-TRANSPORT-2026. NIST Chemistry WebBook SRD 69,
       fluid pages ("References and Notes": viscosity and thermal
       conductivity models), https://webbook.nist.gov/chemistry/fluid/
.. [7] COOLPROP-TRANSPORT-C00300E. CoolProp, commit
       c00300e97ad4cdb30032635d7379209d5e45a121: ``dev/fluids/Hydrogen.json``,
       ``dev/fluids/ParaHydrogen.json`` and
       ``src/Backends/Helmholtz/TransportRoutines.cpp``,
       https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from math import atan, exp, log, pi, sqrt

from rocket_propulsion.core.constants import BOLTZMANN_CONSTANT_J_K
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

from .helmholtz import FluidModel, FluidState, RealFluid
from .reference_fluids import reference_fluid

# Boltzmann constant as printed by Lemmon and Jacobsen (2004), p. 27, and by
# Friend, Ely and Ingham (1989), Table 1 (CODATA 1986); the model's own value.
_K_B_CODATA1986 = 1.380658e-23
# Exact SI value (BIPM SI Brochure, 9th ed., 2019; registered in core.constants);
# Assael et al. (2011) name Boltzmann's constant without printing a value.
_K_B_SI2019 = BOLTZMANN_CONSTANT_J_K
_MICRO = 1e-6
_MILLI = 1e-3
_MPA = 1e6
# Universal constants of the simplified Olchowy-Sengers model: R_D, nu, gamma,
# as printed in Lemmon and Jacobsen (2004), p. 27, and Assael et al. (2011), p. 6.
_OS_RD = 1.01
_OS_NU = 0.63
_OS_GAMMA = 1.2415
_CRITICAL_WARN_SHARE = 0.1

#: Fluids with at least one transport correlation.
TRANSPORT_FLUID_KEYS: tuple[str, ...] = (
    "parahydrogen", "normal_hydrogen", "oxygen", "nitrogen", "methane", "helium")


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise InputError(f"{label} must be positive.", field=label)
    return number


def _sum_terms(terms: tuple[tuple[float, float, float, float], ...], tau: float,
               delta: float) -> float:
    """Return ``sum N tau^t delta^d exp(-gamma delta^l)`` with ``gamma = (l != 0)``."""

    total = 0.0
    for n, t, d, el in terms:
        value = n * tau ** t * delta ** d
        if el != 0.0:
            value *= exp(-delta ** el)
        total += value
    return total


# --------------------------------------------------------------------------- data classes
@dataclass(frozen=True, slots=True)
class OlchowySengersParameters:
    """Fluid constants of the simplified Olchowy-Sengers enhancement.

    Attributes
    ----------
    xi0_m : float
        Correlation-length amplitude ``xi_0``.
    gamma_amplitude : float
        Susceptibility amplitude ``Gamma``.
    qd_inverse_m : float
        Effective cutoff length ``q_D^-1``.
    t_ref_k : float
        Reference temperature where the enhancement is negligible.
    rho_crit_kg_m3, p_crit_pa : float
        Critical density and pressure used to reduce the susceptibility.
    boltzmann_j_k : float
        Boltzmann constant used by the source.

    Notes
    -----
    The enhancement is [1]_ (Eqs. 7-11), [2]_ (Eqs. 4-7)::

        lambda_c = rho cp R_D k T / (6 pi eta xi) (Omega - Omega_0)
        Omega    = 2/pi [((cp - cv)/cp) atan(y) + (cv/cp) y],  y = xi / q_D^-1
        Omega_0  = 2/pi [1 - exp(-1 / (1/y + y^2 (rho_c/rho)^2 / 3))]
        xi       = xi_0 [(chi(T) - chi(T_ref) T_ref/T) / Gamma]^(nu/gamma)
        chi      = p_c rho / rho_c^2 (drho/dp)_T

    with ``R_D = 1.01``, ``nu = 0.63``, ``gamma = 1.2415``; ``lambda_c = 0``
    when the bracket of ``xi`` is not positive [1]_ (p. 27).

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (7)-(11), p. 27,
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eqs. (4)-(7), p. 6,
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    """

    xi0_m: float
    gamma_amplitude: float
    qd_inverse_m: float
    t_ref_k: float
    rho_crit_kg_m3: float
    p_crit_pa: float
    boltzmann_j_k: float

    def enhancement(self, *, density_kg_m3: float, temperature_k: float, cp_j_kg_k: float,
                    cv_j_kg_k: float, viscosity_pa_s: float, drho_dp_t: float,
                    drho_dp_t_ref: float) -> float:
        """Return ``lambda_c`` in W/(m K).

        Parameters
        ----------
        density_kg_m3, temperature_k : float
            State.
        cp_j_kg_k, cv_j_kg_k : float
            Heat capacities at the state.
        viscosity_pa_s : float
            Viscosity at the state.
        drho_dp_t, drho_dp_t_ref : float
            ``(drho/dp)_T`` in kg/(m3 Pa) at ``T`` and at ``T_ref`` (same
            density).

        Returns
        -------
        float
            Critical enhancement, zero when the correlation-length bracket is
            not positive.

        References
        ----------
        .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
               Thermophys. 25, 21 (2004), Eqs. (7)-(11), p. 27,
               https://trc.nist.gov/refprop/FAQ/NAO.PDF
        """

        rho = density_kg_m3
        scale = self.p_crit_pa * rho / self.rho_crit_kg_m3 ** 2
        bracket = scale * (drho_dp_t - drho_dp_t_ref * self.t_ref_k / temperature_k)
        bracket /= self.gamma_amplitude
        if bracket <= 0.0:
            return 0.0
        xi = self.xi0_m * bracket ** (_OS_NU / _OS_GAMMA)
        y = xi / self.qd_inverse_m
        ratio = cv_j_kg_k / cp_j_kg_k
        omega = 2.0 / pi * ((1.0 - ratio) * atan(y) + ratio * y)
        rc = self.rho_crit_kg_m3 / rho
        omega0 = 2.0 / pi * (1.0 - exp(-1.0 / (1.0 / y + y * y * rc * rc / 3.0)))
        return (rho * cp_j_kg_k * _OS_RD * self.boltzmann_j_k * temperature_k
                / (6.0 * pi * viscosity_pa_s * xi) * (omega - omega0))


@dataclass(frozen=True, slots=True)
class TransportValidity:
    """Declared range of the transport correlations of one fluid.

    Attributes
    ----------
    fluid : str
        Registry key.
    viscosity_available : bool
        True for every fluid of this module (kept for fluids added later
        without a viscosity correlation).
    viscosity_t_min_k, viscosity_t_max_k, viscosity_p_max_pa : float
        Range of the viscosity correlation.
    conductivity_t_min_k, conductivity_t_max_k, conductivity_p_max_pa : float
        Range of the thermal-conductivity correlation.
    viscosity_model, conductivity_model : str
        Source key and equations.
    notes : tuple of str
        How the range was set.
    """

    fluid: str
    viscosity_available: bool
    viscosity_t_min_k: float
    viscosity_t_max_k: float
    viscosity_p_max_pa: float
    conductivity_t_min_k: float
    conductivity_t_max_k: float
    conductivity_p_max_pa: float
    viscosity_model: str
    conductivity_model: str
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class TransportProperties:
    """Viscosity, thermal conductivity and Prandtl number at one state.

    Attributes
    ----------
    fluid : str
        Registry key.
    phase : {"liquid", "vapor", "supercritical"}
        Phase label of the equation-of-state state.
    temperature_k, pressure_pa, density_kg_m3 : float
        State (pressure and density from the reference equation of state).
    cp_j_kg_k : float
        Isobaric heat capacity from the equation of state.
    viscosity_pa_s : float
        Total viscosity (from the correlation, or caller-supplied).
    viscosity_dilute_pa_s, viscosity_residual_pa_s : float or None
        Dilute-gas and residual (excess) parts; None when supplied.
    thermal_conductivity_w_m_k : float
        Total thermal conductivity.
    conductivity_dilute_w_m_k, conductivity_residual_w_m_k : float
        Dilute-gas and residual (excess) parts.
    conductivity_critical_w_m_k : float
        Critical enhancement.
    prandtl_number : float
        ``cp mu / k``.
    viscosity_model, conductivity_model : str
        Source keys.
    warnings : tuple of WarningMessage
        Accuracy notes.
    """

    fluid: str
    phase: str
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    cp_j_kg_k: float
    viscosity_pa_s: float
    viscosity_dilute_pa_s: float | None
    viscosity_residual_pa_s: float | None
    thermal_conductivity_w_m_k: float
    conductivity_dilute_w_m_k: float
    conductivity_residual_w_m_k: float
    conductivity_critical_w_m_k: float
    prandtl_number: float
    viscosity_model: str
    conductivity_model: str
    warnings: tuple[WarningMessage, ...] = ()


@dataclass(frozen=True, slots=True)
class _Inputs:
    """Equation-of-state quantities the correlations need at ``(T, rho)``."""

    model: FluidModel
    state: FluidState
    drho_dp_t: float
    drho_dp_t_ref: float


# --------------------------------------------------------------------------- Lemmon-Jacobsen
@dataclass(frozen=True, slots=True)
class _LemmonJacobsen:
    """Coefficients of Lemmon and Jacobsen (2004), Tables I-IV."""

    t_crit_k: float
    rho_crit_mol_m3: float
    molar_mass_g_mol: float
    epsilon_k: float
    sigma_nm: float
    viscosity_terms: tuple[tuple[float, float, float, float], ...]
    n1: float
    n2: float
    t2: float
    n3: float
    t3: float
    conductivity_terms: tuple[tuple[float, float, float, float], ...]
    critical: OlchowySengersParameters


# Table II: collision-integral coefficients b_0..b_4 (all four fluids).
_LJ_B = (0.431, -0.4623, 0.08406, 0.005341, -0.00331)


def _lj_viscosity_parts(c: _LemmonJacobsen, temperature: float,
                        rho_mol: float) -> tuple[float, float]:
    """Return ``(eta_0, eta_r)`` in uPa s, Eqs. (1)-(3)."""

    ln_ts = log(temperature / c.epsilon_k)
    omega = exp(sum(b * ln_ts ** i for i, b in enumerate(_LJ_B)))
    eta0 = 0.0266958 * sqrt(c.molar_mass_g_mol * temperature) / (c.sigma_nm ** 2 * omega)
    tau = c.t_crit_k / temperature
    delta = rho_mol / c.rho_crit_mol_m3
    return eta0, _sum_terms(c.viscosity_terms, tau, delta)


_LJ_N2 = _LemmonJacobsen(
    t_crit_k=126.192, rho_crit_mol_m3=11.1839e3, molar_mass_g_mol=28.01348,
    epsilon_k=98.94, sigma_nm=0.3656,
    viscosity_terms=((10.72, 0.1, 2.0, 0.0), (0.03989, 0.25, 10.0, 1.0),
                     (0.001208, 3.2, 12.0, 1.0), (-7.402, 0.9, 2.0, 2.0),
                     (4.620, 0.3, 1.0, 3.0)),
    n1=1.511, n2=2.117, t2=-1.0, n3=-3.332, t3=-0.7,
    conductivity_terms=((8.862, 0.0, 1.0, 0.0), (31.11, 0.03, 2.0, 0.0),
                        (-73.13, 0.2, 3.0, 1.0), (20.03, 0.8, 4.0, 2.0),
                        (-0.7096, 0.6, 8.0, 2.0), (0.2672, 1.9, 10.0, 2.0)),
    critical=OlchowySengersParameters(
        xi0_m=0.17e-9, gamma_amplitude=0.055, qd_inverse_m=0.40e-9, t_ref_k=252.384,
        rho_crit_kg_m3=11.1839e3 * 28.01348e-3, p_crit_pa=3.3958 * _MPA,
        boltzmann_j_k=_K_B_CODATA1986))

_LJ_O2 = _LemmonJacobsen(
    t_crit_k=154.581, rho_crit_mol_m3=13.63e3, molar_mass_g_mol=31.9988,
    epsilon_k=118.5, sigma_nm=0.3428,
    viscosity_terms=((17.67, 0.05, 1.0, 0.0), (0.4042, 0.0, 5.0, 0.0),
                     (0.0001077, 2.10, 12.0, 0.0), (0.3510, 0.0, 8.0, 1.0),
                     (-13.67, 0.5, 1.0, 2.0)),
    n1=1.036, n2=6.283, t2=-0.9, n3=-4.262, t3=-0.6,
    conductivity_terms=((15.31, 0.0, 1.0, 0.0), (8.898, 0.0, 3.0, 0.0),
                        (-0.7336, 0.3, 4.0, 0.0), (6.728, 4.3, 5.0, 2.0),
                        (-4.374, 0.5, 7.0, 2.0), (-0.4747, 1.8, 10.0, 2.0)),
    critical=OlchowySengersParameters(
        xi0_m=0.24e-9, gamma_amplitude=0.055, qd_inverse_m=0.51e-9, t_ref_k=309.162,
        rho_crit_kg_m3=13.63e3 * 31.9988e-3, p_crit_pa=5.043 * _MPA,
        boltzmann_j_k=_K_B_CODATA1986))


# --------------------------------------------------------------------------- Assael 2011
@dataclass(frozen=True, slots=True)
class _Assael:
    """Coefficients of Assael et al. (2011), Tables 2 and 5."""

    t_crit_k: float
    rho_crit_kg_m3: float
    a1: tuple[float, ...]
    a2: tuple[float, ...]
    b1: tuple[float, ...]
    b2: tuple[float, ...]
    critical: OlchowySengersParameters


def _assael_dilute(c: _Assael, temperature: float) -> float:
    tr = temperature / c.t_crit_k
    num = sum(a * tr ** i for i, a in enumerate(c.a1))
    den = sum(a * tr ** i for i, a in enumerate(c.a2))
    return num / den


_ASSAEL_NORMAL = _Assael(
    t_crit_k=33.145, rho_crit_kg_m3=31.262,
    a1=(-3.40976e-1, 4.58820e0, -1.45080e0, 3.26394e-1, 3.16939e-3, 1.90592e-4, -1.13900e-6),
    a2=(1.38497e2, -2.21878e1, 4.57151e0, 1.00000e0),
    b1=(3.63081e-2, -2.07629e-2, 3.14810e-2, -1.43097e-2, 1.74980e-3),
    b2=(1.83370e-3, -8.86716e-3, 1.58260e-2, -1.06283e-2, 2.80673e-3),
    critical=OlchowySengersParameters(
        xi0_m=1.5e-10, gamma_amplitude=0.052, qd_inverse_m=4.0e-10, t_ref_k=49.7175,
        rho_crit_kg_m3=31.262, p_crit_pa=1.2964 * _MPA, boltzmann_j_k=_K_B_SI2019))

_ASSAEL_PARA = _Assael(
    t_crit_k=32.938, rho_crit_kg_m3=31.323,
    a1=(-1.24500e0, 3.10212e2, -3.31004e2, 2.46016e2, -6.57810e1, 1.08260e1, -5.19659e-1,
        1.43979e-2),
    a2=(1.42304e4, -1.93922e4, 1.58379e4, -4.81812e3, 7.28639e2, -3.57365e1, 1.00000e0),
    b1=(2.65975e-2, -1.33826e-3, 1.30219e-2, -5.67678e-3, -9.23380e-5),
    b2=(-1.21727e-3, 3.66663e-3, 3.88715e-3, -9.21055e-3, 4.00723e-3),
    critical=OlchowySengersParameters(
        xi0_m=1.5e-10, gamma_amplitude=0.052, qd_inverse_m=5.0e-10, t_ref_k=49.407,
        rho_crit_kg_m3=31.323, p_crit_pa=1.2858 * _MPA, boltzmann_j_k=_K_B_SI2019))


# --------------------------------------------------------------------------- Friend-Ely-Ingham
# Table 1: critical point and potential parameters; Table 8: Eq. (12) C_1..C_9 and
# Eq. (14) f_1, f_2; Table 9: Eq. (15) (r, s, g) and Eq. (17) (r, s, j); Table 10.
_FEI_TC = 190.551
_FEI_RHOC_MOL_M3 = 10.139e3
_FEI_EPS_K = 174.0
_FEI_R = 8.314510
_FEI_C = (-3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140,
          8.9488430959, -1.8739245042, 0.20966101390, -9.6570437074e-3)
_FEI_F1 = 1.458850
_FEI_F2 = -0.4377162
_FEI_ETA_NUM = ((1.0, 0.0, 0.41250137), (1.0, 1.0, -0.14390912), (2.0, 0.0, 0.10366993),
                (2.0, 1.0, 0.40287464), (2.0, 1.5, -0.24903524), (3.0, 0.0, -0.12953131),
                (3.0, 2.0, 0.06575776), (4.0, 0.0, 0.02566628), (4.0, 1.0, -0.03716526))
_FEI_ETA_DEN = ((1.0, 0.0, -0.38798341), (1.0, 1.0, 0.03533815))
_FEI_LAM = ((1.0, 0.0, 2.4149207), (3.0, 0.0, 0.55166331), (4.0, 0.0, -0.52837734),
            (4.0, 1.0, 0.073809553), (5.0, 0.0, 0.24465507), (5.0, 1.0, -0.047613626))
_FEI_LAM7 = 1.5554612
_FEI_FT = 2.646
_FEI_FRHO = 2.678
_FEI_FA = -0.637
_FEI_P_C = 4.5992 * _MPA
_FEI_SCALED_T = 0.03
_FEI_SCALED_RHO = 0.25


def _fei_dilute_viscosity(temperature: float) -> float:
    """Eqs. (10a) and (12): ``eta_0`` in uPa s."""

    t = temperature / _FEI_EPS_K
    inv_omega = sum(c * t ** (i / 3.0 - 1.0) for i, c in enumerate(_FEI_C))
    return 10.50 * sqrt(t) * inv_omega


def _fei_excess_viscosity(temperature: float, rho_mol: float) -> float:
    """Eq. (15a): ``eta_ex`` in uPa s."""

    tau = _FEI_TC / temperature
    delta = rho_mol / _FEI_RHOC_MOL_M3
    num = sum(g * delta ** r * tau ** s for r, s, g in _FEI_ETA_NUM)
    den = 1.0 + sum(g * delta ** r * tau ** s for r, s, g in _FEI_ETA_DEN)
    return 12.149 * num / den


def _fei_dilute_conductivity(temperature: float, cv0_over_r: float) -> float:
    """Eqs. (13a)-(14): ``lambda_0`` in mW/(m K); ``tau^2 phi_tautau = -cv0/R``."""

    f_int = _FEI_F1 + _FEI_F2 / (temperature / _FEI_EPS_K)
    return 0.51826 * _fei_dilute_viscosity(temperature) * (
        3.75 - f_int * (-cv0_over_r + 1.5))


# --------------------------------------------------------------------------- Muzny 2013 (H2)
# Transcribed from CoolProp commit c00300e (Hydrogen.json and ParaHydrogen.json carry identical
# "viscosity" blocks citing Muzny-JCED-2013; TransportRoutines.cpp,
# viscosity_hydrogen_higher_order_hardcoded); the primary paper was not opened. Avogadro
# constant, the exponent -i and the scale density are those of the erratum.
_MUZNY_EPS_K = 30.41
_MUZNY_SIGMA_NM = 0.297
_MUZNY_M_G_MOL = 2.01588
_MUZNY_C_DILUTE = 2.1357e-8  # Pa s, with M in g/mol, T in K, sigma in nm
_MUZNY_A = (0.20963, -0.455274, 0.143602, -0.0335325, 0.00276981)
_MUZNY_B = (-0.187, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292, 0.5166)
_MUZNY_NA = 6.022137e23  # erratum: value recommended at the time of publication
_MUZNY_TC = 33.145
_MUZNY_RHO_SC = 90.909090909  # kg/m3, erratum item 3
_MUZNY_C = (6.43449673e-6, 4.56334068e-2, 2.32797868e-1, 9.58326120e-1, 1.27941189e-1,
            3.63576595e-1)


def _muzny_viscosity_parts(temperature: float, density_kg_m3: float) -> tuple[float, float]:
    """Return ``(eta_0, eta_1 + eta_h)`` of Muzny et al. (2013) in Pa s.

    ``eta_0 = C sqrt(M T) / (sigma^2 S*)`` with ``ln S* = sum a_i (ln T*)^i``;
    ``eta_1 = eta_0 B_eta rho`` with ``B_eta = N_A sigma^3 sum b_i T*^-i``
    (Rainwater-Friend); ``eta_h = c1 rho_r^2 exp(c2 T_r + c3/T_r +
    c4 rho_r^2/(c5 + T_r) + c6 rho_r^6)``, ``T_r = T/33.145``,
    ``rho_r = rho/90.909090909``.
    """

    t_star = temperature / _MUZNY_EPS_K
    ln_t = log(t_star)
    s_star = exp(sum(a * ln_t ** i for i, a in enumerate(_MUZNY_A)))
    eta0 = _MUZNY_C_DILUTE * sqrt(_MUZNY_M_G_MOL * temperature) / (_MUZNY_SIGMA_NM ** 2 * s_star)
    b_star = sum(b * t_star ** (-i) for i, b in enumerate(_MUZNY_B))
    b_eta = _MUZNY_NA * (_MUZNY_SIGMA_NM * 1e-9) ** 3 * b_star
    eta1 = eta0 * b_eta * density_kg_m3 / (_MUZNY_M_G_MOL * 1e-3)
    tr = temperature / _MUZNY_TC
    rr = density_kg_m3 / _MUZNY_RHO_SC
    c = _MUZNY_C
    eta_h = c[0] * rr * rr * exp(c[1] * tr + c[2] / tr + c[3] * rr * rr / (c[4] + tr)
                                 + c[5] * rr ** 6)
    return eta0, eta1 + eta_h


# --------------------------------------------------------------------------- helium
# Viscosity: Arp, McCarty and Friend, NIST TN 1334 (1998), Sec. 4.2, Eqs. (1)-(9), read from
# the rendered pages (coefficients agree with CoolProp's viscosity_helium_hardcoded except
# that CoolProp truncates the last coefficient of Eq. (2) and omits the 100-110 K blend).
_HE_ETA0P = (-0.135311743, 1.00347841, 1.20654649, -0.149564551, 0.0125208416)  # Eq. (2)
_HE_B = (-47.5295259, 87.6799309, -42.0741589, 8.33128289, -0.589252385)  # Eq. (4)
_HE_C = (547.309267, -904.870586, 431.404928, -81.4504854, 5.37008433)  # Eq. (5)
_HE_D = (-1684.39324, 3331.08630, -1632.19172, 308.804413, -20.2936367)  # Eq. (6)
_UG_CM_S = 1e-7  # 1 ug/(cm s) = 1e-7 Pa s
# Conductivity: Hands and Arp (1981) as implemented in CoolProp commit c00300e,
# TransportRoutines.cpp, conductivity_hardcoded_helium (follows the REFPROP source).
_HE_K_RHOC = 68.0
_HE_K_C = (1.862970530e-4, -7.275964435e-7, -1.427549651e-4, 3.290833592e-5, -5.213335363e-8,
           4.492659933e-8, -5.924416513e-9, 7.087321137e-6, -6.013335678e-6, 8.067145814e-7,
           3.995125013e-7)
# x0, E1, E2, beta, gamma, delta, rho_c (kg/m3), T_c (K), p_c (Pa) of the critical term.
_HE_CRIT = (0.392, 2.8461, 0.27156, 0.3554, 1.1743, 4.304, 69.158, 5.18992, 2.2746e5)


def _poly_x(coef: tuple[float, ...], x: float) -> float:
    """``c0/x + c1 + c2 x + c3 x^2 + c4 x^3`` (TN 1334 Eqs. 2, 4-6)."""

    return coef[0] / x + coef[1] + coef[2] * x + coef[3] * x * x + coef[4] * x ** 3


def _helium_viscosity(temperature: float, density_kg_m3: float) -> float:
    """TN 1334 Eqs. (1)-(9) in Pa s; density in g/cm3 inside the equations.

    ``eta = exp(eta0'(x) + eta_E'(rho, x))`` with ``x = ln T`` for
    ``T <= 100 K``. Above 100 K the dilute part of Steward's equations is
    replaced by Maitland and Smith, Eq. (8), linearly blended in temperature
    between 100 and 110 K (Sec. 4.2.2); above 300 K the excess function is
    evaluated at 300 K (Sec. 4.2.3).
    """

    rho = density_kg_m3 * 1e-3
    x = log(min(temperature, 300.0))
    e0p = _poly_x(_HE_ETA0P, x)
    e_excess = (rho * _poly_x(_HE_B, x) + rho * rho * _poly_x(_HE_C, x)
                + rho ** 3 * _poly_x(_HE_D, x))
    if temperature <= 100.0:
        return exp(e0p + e_excess) * _UG_CM_S
    t = temperature
    maitland = 196.0 * t ** 0.71938 * exp(12.451 / t - 295.67 / (t * t) - 4.1249)
    if t < 110.0:
        w = (t - 100.0) / 10.0
        dilute = (1.0 - w) * exp(e0p) + w * maitland
    else:
        dilute = maitland
    return (exp(e0p + e_excess) - exp(e0p) + dilute) * _UG_CM_S


def _helium_dilute_conductivity(temperature: float) -> float:
    """Hands-Arp dilute-gas conductivity in W/(m K) (CoolProp transcription)."""

    inv = 1.0 / temperature
    summer = (3.739232544 * inv - 2.620316969e1 * inv ** 2 + 5.982252246e1 * inv ** 3
              - 4.926397634e1 * inv ** 4)
    return 2.7870034e-3 * temperature ** 7.034007057e-1 * exp(summer)


def _helium_conductivity_parts(state: FluidState, drho_dp_t: float,
                               viscosity_pa_s: float) -> tuple[float, float, float]:
    """Hands-Arp conductivity ``(lambda_0, lambda_e, lambda_c)`` in W/(m K).

    Transcribed from CoolProp ``conductivity_hardcoded_helium`` (commit
    c00300e), which follows the REFPROP source; the paper was not opened.
    At exactly ``rho = 69.158 kg/m3`` the critical term's density distance
    is floored at 1e-12 to avoid a division by zero.
    """

    t, rho = state.temperature_k, state.density_kg_m3
    lam0 = _helium_dilute_conductivity(t)
    c = _HE_K_C
    t13, t23 = t ** (1.0 / 3.0), t ** (2.0 / 3.0)
    lame = ((c[0] + c[1] * t + c[2] * t13 + c[3] * t23) * rho
            + (c[4] + c[5] * t13 + c[6] * t23) * rho ** 3
            + (c[7] + c[8] * t13 + c[9] * t23 + c[10] / t) * rho * rho * log(rho / _HE_K_RHOC))
    lamc = 0.0
    if 3.5 < t < 12.0:
        x0, e1, e2, beta, gamma, delta_exp, rhoc, tc, pc = _HE_CRIT
        d_t = abs(1.0 - t / tc)
        d_rho = max(abs(1.0 - rho / rhoc), 1e-12)
        k_t = drho_dp_t / rho
        weight = (d_t / 0.2) ** 2 + (d_rho / 0.25) ** 2
        if weight > 1.0:
            k_bar = k_t
        else:
            x = (d_t / d_rho) ** (1.0 / beta)
            u = 1.0 + x / x0
            h = e1 * u * (1.0 + e2 * u ** (2.0 / beta)) ** ((gamma - 1.0) / (2.0 * beta))
            q = u ** (2.0 / beta)
            dhdx = e1 * (e2 * q * (gamma - 1.0) * (e2 * q + 1.0) ** (0.5 * (gamma - 1.0) / beta)
                         + beta * beta * (e2 * q + 1.0) ** (
                             0.5 * (2.0 * beta + gamma - 1.0) / beta)) / (
                beta * beta * x0 * (e2 * q + 1.0))
            rhs = d_rho ** (delta_exp - 1.0) * (delta_exp * h - x / beta * dhdx)
            k_prime = 1.0 / (rhs * (rho / rhoc) ** 2 * pc)
            k_bar = weight * k_t + (1.0 - weight) * k_prime
        lamc = (3.4685233e-17 * 3.726229668 * sqrt(k_bar) * t * t / rho / viscosity_pa_s
                * state.dp_dt_rho_pa_k ** 2 * exp(-18.66 * d_t * d_t - 4.25 * d_rho ** 4))
    return lam0, lame, lamc


# --------------------------------------------------------------------------- registry
_MUZNY_MODEL = ("Muzny, Huber and Kazakov (2013) via COOLPROP-TRANSPORT-C00300E with "
                "MUZNY-2022-H2-VISCOSITY-ERRATUM (paper not opened)")
_MUZNY_RANGE_NOTE = ("Muzny et al. (2013) range as stated in the NIST publication record: "
                     "triple point to 1000 K, pressures up to 200 MPa.")
_VALIDITY = {
    "nitrogen": TransportValidity(
        "nitrogen", True, 63.151, 1000.0, 2200.0 * _MPA, 63.151, 1000.0, 2200.0 * _MPA,
        "LEMMON-JACOBSEN-2004-TRANSPORT Eqs. (1)-(3)",
        "LEMMON-JACOBSEN-2004-TRANSPORT Eqs. (4)-(11)",
        (("The paper states validity over all liquid and vapor states (abstract); the range "
         "of the Span et al. (2000) equation of state is used. Uncertainties above 100 MPa "
         "are stated as unknown (p. 62)."),)),
    "oxygen": TransportValidity(
        "oxygen", True, 54.361, 300.0, 80.0 * _MPA, 54.361, 300.0, 80.0 * _MPA,
        "LEMMON-JACOBSEN-2004-TRANSPORT Eqs. (1)-(3)",
        "LEMMON-JACOBSEN-2004-TRANSPORT Eqs. (4)-(11)",
        (("Limited by the Schmidt-Wagner equation of state used for the density "
         "(triple point to 300 K, 80 MPa)."),)),
    "parahydrogen": TransportValidity(
        "parahydrogen", True, 13.8033, 1000.0, 200.0 * _MPA,
        13.8033, 1000.0, 100.0 * _MPA, _MUZNY_MODEL,
        "ASSAEL-2011-H2-CONDUCTIVITY Eqs. (1)-(7), Table 5",
        ("Title range of Assael et al. (2011): triple point to 1000 K, up to 100 MPa.",
         _MUZNY_RANGE_NOTE)),
    "normal_hydrogen": TransportValidity(
        "normal_hydrogen", True, 13.957, 1000.0, 200.0 * _MPA,
        13.957, 1000.0, 100.0 * _MPA, _MUZNY_MODEL,
        "ASSAEL-2011-H2-CONDUCTIVITY Eqs. (1)-(7), Table 2",
        ("Title range of Assael et al. (2011): triple point to 1000 K, up to 100 MPa.",
         _MUZNY_RANGE_NOTE)),
    "methane": TransportValidity(
        "methane", True, 90.6941, 400.0, 55.0 * _MPA, 90.6941, 625.0, 100.0 * _MPA,
        "FRIEND-ELY-INGHAM-1989-CH4 Eqs. (10a), (12), (15a)",
        "FRIEND-ELY-INGHAM-1989-CH4 Eqs. (13a), (14), (17a), (18a)-(22)",
        (("Abstract of Friend et al. (1989): viscosity about 91-400 K below 55 MPa, "
         "conductivity about 91-700 K below 100 MPa; the upper temperature of the "
         "conductivity is capped at 625 K by the Setzmann-Wagner equation of state."),
         ("Above the viscosity range the conductivity critical term uses the extrapolated "
         "Friend-Ely-Ingham viscosity."))),
    "helium": TransportValidity(
        "helium", True, 2.5, 1500.0, 100.0 * _MPA, 3.5, 1500.0, 100.0 * _MPA,
        "ARP-MCCARTY-FRIEND-1998-HE Eqs. (1)-(9)",
        "Hands and Arp (1981) via COOLPROP-TRANSPORT-C00300E (paper not opened)",
        (("TN 1334, Sec. 1: the helium I code is valid from 2.5 to 1500 K to 100 MPa; "
         "Sec. 4.1: the conductivity validity is unknown below 3.5 K or above 100 MPa."),)),
}


def _key(name: object) -> str:
    if not isinstance(name, str) or name not in _VALIDITY:
        raise InputError(
            f"No transport correlation for {name!r}; choose one of "
            f"{', '.join(TRANSPORT_FLUID_KEYS)}.", field="fluid")
    return name


@cache
def _fluid(key: str) -> RealFluid:
    return reference_fluid(key)


def transport_validity(fluid: str) -> TransportValidity:
    """Return the declared range of the transport correlations of ``fluid``.

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.

    Returns
    -------
    TransportValidity
        Temperature and pressure limits and model keys.

    Raises
    ------
    InputError
        If the fluid has no transport correlation.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), abstract and p. 62,
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), title,
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), abstract,
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    return _VALIDITY[_key(fluid)]


# --------------------------------------------------------------------------- evaluation
def _drho_dp_unchecked(model: FluidModel, density: float, temperature: float) -> float:
    """``(drho/dp)_T`` in kg/(m3 Pa) without range checks (reference temperature)."""

    delta = density / model.molar_mass_kg_mol / model.rho_reduce_mol_m3
    tau = model.t_reduce_k / temperature
    _ar, ar_d, ar_dd, *_ = model.residual.derivatives(tau, delta)
    dp_drho_molar = model.gas_constant_j_mol_k * temperature * (
        1.0 + 2.0 * delta * ar_d + delta * delta * ar_dd)
    return model.molar_mass_kg_mol / dp_drho_molar


def _inputs(key: str, state: FluidState) -> _Inputs:
    model = _fluid(key).model
    t_ref = {"nitrogen": _LJ_N2, "oxygen": _LJ_O2}.get(key)
    if t_ref is not None:
        ref_t = t_ref.critical.t_ref_k
    elif key == "parahydrogen":
        ref_t = _ASSAEL_PARA.critical.t_ref_k
    elif key == "normal_hydrogen":
        ref_t = _ASSAEL_NORMAL.critical.t_ref_k
    else:
        ref_t = state.temperature_k
    return _Inputs(model=model, state=state, drho_dp_t=1.0 / state.dp_drho_t_pa_m3_kg,
                   drho_dp_t_ref=_drho_dp_unchecked(model, state.density_kg_m3, ref_t))


def _check_range(key: str, state: FluidState, which: str) -> None:
    v = _VALIDITY[key]
    if which == "viscosity":
        t_lo, t_hi, p_hi = v.viscosity_t_min_k, v.viscosity_t_max_k, v.viscosity_p_max_pa
    else:
        t_lo, t_hi, p_hi = (v.conductivity_t_min_k, v.conductivity_t_max_k,
                            v.conductivity_p_max_pa)
    t, p = state.temperature_k, state.pressure_pa
    if not (t_lo <= t <= t_hi and p <= p_hi):
        raise DomainError(
            f"{key}: T = {t} K, p = {p} Pa is outside the {which} correlation range "
            f"[{t_lo}, {t_hi}] K, p <= {p_hi} Pa.", field="temperature_k")


def _viscosity_parts(key: str, state: FluidState) -> tuple[float, float]:
    """Return ``(dilute, residual)`` viscosity in Pa s (no range check)."""

    t = state.temperature_k
    if key in ("nitrogen", "oxygen"):
        c = _LJ_N2 if key == "nitrogen" else _LJ_O2
        eta0, etar = _lj_viscosity_parts(c, t, state.molar_density_mol_m3)
    elif key == "methane":
        eta0 = _fei_dilute_viscosity(t)
        etar = _fei_excess_viscosity(t, state.molar_density_mol_m3)
    elif key == "helium":
        dilute = _helium_viscosity(t, 0.0)
        return dilute, _helium_viscosity(t, state.density_kg_m3) - dilute
    else:
        return _muzny_viscosity_parts(t, state.density_kg_m3)
    return eta0 * _MICRO, etar * _MICRO


def _conductivity_parts(key: str, inp: _Inputs,
                        viscosity_pa_s: float) -> tuple[float, float, float]:
    """Return ``(dilute, residual, critical)`` conductivity in W/(m K)."""

    s = inp.state
    t, rho = s.temperature_k, s.density_kg_m3
    if key in ("nitrogen", "oxygen"):
        c = _LJ_N2 if key == "nitrogen" else _LJ_O2
        eta0, _ = _lj_viscosity_parts(c, t, s.molar_density_mol_m3)
        tau = c.t_crit_k / t
        delta = s.molar_density_mol_m3 / c.rho_crit_mol_m3
        lam0 = (c.n1 * eta0 + c.n2 * tau ** c.t2 + c.n3 * tau ** c.t3) * _MILLI
        lamr = _sum_terms(c.conductivity_terms, tau, delta) * _MILLI
        lamc = c.critical.enhancement(
            density_kg_m3=rho, temperature_k=t, cp_j_kg_k=s.cp_j_kg_k, cv_j_kg_k=s.cv_j_kg_k,
            viscosity_pa_s=viscosity_pa_s, drho_dp_t=inp.drho_dp_t,
            drho_dp_t_ref=inp.drho_dp_t_ref)
        return lam0, lamr, lamc
    if key in ("parahydrogen", "normal_hydrogen"):
        c2 = _ASSAEL_PARA if key == "parahydrogen" else _ASSAEL_NORMAL
        tr = t / c2.t_crit_k
        dr = rho / c2.rho_crit_kg_m3
        lam0 = _assael_dilute(c2, t)
        lamr = sum((b1 + b2 * tr) * dr ** (i + 1)
                   for i, (b1, b2) in enumerate(zip(c2.b1, c2.b2, strict=True)))
        lamc = c2.critical.enhancement(
            density_kg_m3=rho, temperature_k=t, cp_j_kg_k=s.cp_j_kg_k, cv_j_kg_k=s.cv_j_kg_k,
            viscosity_pa_s=viscosity_pa_s, drho_dp_t=inp.drho_dp_t,
            drho_dp_t_ref=inp.drho_dp_t_ref)
        return lam0, lamr, lamc
    if key == "helium":
        return _helium_conductivity_parts(s, inp.drho_dp_t, viscosity_pa_s)
    return _fei_conductivity_parts(inp, viscosity_pa_s)


def _fei_conductivity_parts(inp: _Inputs, viscosity_pa_s: float) -> tuple[float, float, float]:
    """Friend-Ely-Ingham ``(lambda_0, lambda_ex, lambda_cr)`` in W/(m K).

    Eq. (16) of the reprint prints ``delta*_sigma = 11`` outside the vapor
    branch; the paper's own Table A3 (p. 632) is reproduced only with
    ``delta*_sigma = 1`` (100 K, 0.1 MPa: 199.7 versus 199.8 mW/(m K); with
    11 the result is 134.8), so 1 is used.
    """

    s, m = inp.state, inp.model
    t = s.temperature_k
    tau_eos = m.t_reduce_k / t
    delta_eos = s.molar_density_mol_m3 / m.rho_reduce_mol_m3
    _a0, _a0_t, a0_tt = m.ideal.evaluate(tau_eos, delta_eos)
    lam0 = _fei_dilute_conductivity(t, -tau_eos * tau_eos * a0_tt)
    tau = _FEI_TC / t
    delta = s.molar_density_mol_m3 / _FEI_RHOC_MOL_M3
    if t < _FEI_TC and delta < 1.0:
        delta_sigma = _fluid("methane").saturation_t(t).vapor.molar_density_mol_m3 / (
            _FEI_RHOC_MOL_M3)
    else:
        delta_sigma = 1.0
    lamex = 6.29638 * (sum(j * delta ** r * tau ** sx for r, sx, j in _FEI_LAM)
                       + _FEI_LAM7 * delta * delta / delta_sigma)
    # Eq. (18a) with the physical forms of its reduced derivatives: the bracket is
    # (dp/dT)_rho / (rho R) and Eq. (19) chi* = p_c / rho_c^2 * rho (drho/dp)_T.
    thermal = s.dp_dt_rho_pa_k / (s.molar_density_mol_m3 * _FEI_R)
    drho_dp_molar = inp.drho_dp_t / m.molar_mass_kg_mol
    chi = _FEI_P_C / _FEI_RHOC_MOL_M3 ** 2 * s.molar_density_mol_m3 * drho_dp_molar
    t_star = 1.0 - t / _FEI_TC
    rho_star = 1.0 - delta
    damping = exp(-(_FEI_FT * sqrt(abs(t_star)) + _FEI_FRHO * rho_star * rho_star
                    + _FEI_FA * rho_star))
    lamcr = (91.855 / (viscosity_pa_s / _MICRO * tau * tau) * thermal * thermal
             * chi ** 0.4681 * damping)
    return lam0 * _MILLI, lamex * _MILLI, lamcr * _MILLI


def _warnings(key: str, state: FluidState, supplied: bool,
              critical_share: float) -> tuple[WarningMessage, ...]:
    found: list[WarningMessage] = []
    if critical_share > _CRITICAL_WARN_SHARE:
        found.append(WarningMessage(
            "transport_near_critical",
            f"The critical enhancement is {100.0 * critical_share:.0f} % of the thermal "
            "conductivity; the correlations are less certain in the critical region."))
    if key == "methane":
        t_star = abs(1.0 - state.temperature_k / _FEI_TC)
        rho_star = abs(1.0 - state.molar_density_mol_m3 / _FEI_RHOC_MOL_M3)
        if t_star < _FEI_SCALED_T and rho_star < _FEI_SCALED_RHO:
            found.append(WarningMessage(
                "methane_scaled_region_not_modelled",
                "Friend et al. (1989) evaluate the susceptibility with a scaled equation of "
                "state for |T*| < 0.03 and |rho*| < 0.25; this implementation (like the NIST "
                "WebBook values it was compared with) uses the Setzmann-Wagner equation of "
                "state there.", severity="info"))
        v = _VALIDITY[key]
        outside = (state.temperature_k > v.viscosity_t_max_k
                   or state.pressure_pa > v.viscosity_p_max_pa)
        if outside and not supplied:
            found.append(WarningMessage(
                "viscosity_extrapolated_in_critical_term",
                "The conductivity critical term uses the Friend-Ely-Ingham viscosity outside "
                "its stated range.", severity="info"))
    if key == "nitrogen" and state.pressure_pa > 100.0 * _MPA:
        found.append(WarningMessage(
            "transport_uncertainty_unknown",
            "Lemmon and Jacobsen (2004) state that uncertainties above 100 MPa are unknown."))
    if supplied:
        found.append(WarningMessage(
            "viscosity_caller_supplied",
            "The viscosity was supplied by the caller; it also enters the conductivity "
            "critical enhancement.", severity="info"))
    return tuple(found)


def _evaluate(key: str, state: FluidState, viscosity_pa_s: object,
              *, need_viscosity: bool) -> TransportProperties:
    supplied = viscosity_pa_s is not None
    _check_range(key, state, "conductivity")
    if supplied:
        eta = _positive(viscosity_pa_s, "viscosity_pa_s")
        eta0 = etar = None
    else:
        if need_viscosity:
            _check_range(key, state, "viscosity")
        eta0, etar = _viscosity_parts(key, state)
        eta = eta0 + etar
    inp = _inputs(key, state)
    lam0, lamr, lamc = _conductivity_parts(key, inp, eta)
    lam = lam0 + lamr + lamc
    v = _VALIDITY[key]
    return TransportProperties(
        fluid=key, phase=state.phase, temperature_k=state.temperature_k,
        pressure_pa=state.pressure_pa, density_kg_m3=state.density_kg_m3,
        cp_j_kg_k=state.cp_j_kg_k, viscosity_pa_s=eta, viscosity_dilute_pa_s=eta0,
        viscosity_residual_pa_s=etar, thermal_conductivity_w_m_k=lam,
        conductivity_dilute_w_m_k=lam0, conductivity_residual_w_m_k=lamr,
        conductivity_critical_w_m_k=lamc, prandtl_number=state.cp_j_kg_k * eta / lam,
        viscosity_model="caller-supplied" if supplied else v.viscosity_model,
        conductivity_model=v.conductivity_model,
        warnings=state.warnings + _warnings(key, state, supplied, lamc / lam))


def transport_from_state(fluid: str, state: FluidState, *,
                         viscosity_pa_s: float | None = None) -> TransportProperties:
    """Return transport properties at a state of the reference equation of state.

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    state : FluidState
        State returned by ``reference_fluid(fluid).state_tp`` or
        ``state_rho_t``.
    viscosity_pa_s : float, optional
        Optional override. When given it replaces the viscosity correlation
        and is used in the conductivity critical enhancement.

    Returns
    -------
    TransportProperties
        Viscosity, thermal conductivity (with its parts) and Prandtl number.

    Raises
    ------
    InputError
        Unknown fluid, a state of another fluid, or a non-positive
        ``viscosity_pa_s``.
    DomainError
        If the state is outside a correlation's declared range.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (1)-(11),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eqs. (1)-(7),
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (9)-(22),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    if not isinstance(state, FluidState) or state.fluid != _fluid(key).model.name \
            or state.model_kind != "reference":
        raise InputError(f"state must be a reference-equation state of {key}.", field="state")
    return _evaluate(key, state, viscosity_pa_s, need_viscosity=True)


def transport_rho_t(fluid: str, density_kg_m3: float, temperature_k: float, *,
                    viscosity_pa_s: float | None = None) -> TransportProperties:
    """Return transport properties at density and temperature.

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    density_kg_m3, temperature_k : float
        Single-phase state (two-phase states are refused by the equation of
        state).
    viscosity_pa_s : float, optional
        Optional viscosity override (also used in the critical term).

    Returns
    -------
    TransportProperties
        Viscosity, thermal conductivity and Prandtl number.

    Raises
    ------
    InputError
        Unknown fluid or a non-positive ``viscosity_pa_s``.
    DomainError
        Outside the equation-of-state or correlation range, or two-phase.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (1)-(11),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eqs. (1)-(7),
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (9)-(22),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    state = _fluid(key).state_rho_t(density_kg_m3, temperature_k)
    return _evaluate(key, state, viscosity_pa_s, need_viscosity=True)


def transport_tp(fluid: str, temperature_k: float, pressure_pa: float, *,
                 viscosity_pa_s: float | None = None) -> TransportProperties:
    """Return viscosity, thermal conductivity and Prandtl number at ``(T, p)``.

    The density, ``cp``, ``cv`` and ``(drho/dp)_T`` come from the reference
    equation of state of the same fluid; see
    :func:`rocket_propulsion.thermodynamics.real_fluids.reference_fluid`.

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    temperature_k, pressure_pa : float
        Single-phase state.
    viscosity_pa_s : float, optional
        Optional viscosity override (also used in the critical term).

    Returns
    -------
    TransportProperties
        Viscosity, thermal conductivity and Prandtl number with the
        equation-of-state density and ``cp``.

    Raises
    ------
    InputError
        Unknown fluid or a non-positive ``viscosity_pa_s``.
    DomainError
        Outside the equation-of-state or correlation range, or on the
        saturation curve.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (1)-(11),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eqs. (1)-(7),
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (9)-(22),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    state = _fluid(key).state_tp(temperature_k, pressure_pa)
    return _evaluate(key, state, viscosity_pa_s, need_viscosity=True)


def thermal_conductivity_rho_t(fluid: str, density_kg_m3: float, temperature_k: float, *,
                               viscosity_pa_s: float | None = None) -> float:
    """Return the thermal conductivity in W/(m K) at density and temperature.

    Unlike :func:`transport_rho_t` this checks only the conductivity range;
    for methane above the viscosity range the critical term uses the
    extrapolated viscosity (with an ``info`` warning in the full result).

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    density_kg_m3, temperature_k : float
        Single-phase state.
    viscosity_pa_s : float, optional
        Optional viscosity override for the critical enhancement.

    Returns
    -------
    float
        Thermal conductivity.

    Raises
    ------
    InputError
        Unknown fluid or a non-positive ``viscosity_pa_s``.
    DomainError
        Outside the equation-of-state or conductivity range.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (4)-(11),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eqs. (1)-(7),
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (13a)-(22),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    state = _fluid(key).state_rho_t(density_kg_m3, temperature_k)
    return _evaluate(key, state, viscosity_pa_s,
                     need_viscosity=False).thermal_conductivity_w_m_k


def viscosity_rho_t(fluid: str, density_kg_m3: float, temperature_k: float) -> float:
    """Return the viscosity in Pa s at density and temperature.

    Checks only the viscosity range, which for hydrogen (200 MPa) is wider
    than the conductivity range (100 MPa).

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    density_kg_m3, temperature_k : float
        Single-phase state.

    Returns
    -------
    float
        Viscosity.

    Raises
    ------
    InputError
        If the fluid is unknown.
    DomainError
        Outside the equation-of-state or viscosity range, or two-phase.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (1)-(3),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (10a), (12), (15a),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [3] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity of Muzny et al. (2013) (primary paper not opened),
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [4] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    state = _fluid(key).state_rho_t(density_kg_m3, temperature_k)
    _check_range(key, state, "viscosity")
    eta0, etar = _viscosity_parts(key, state)
    return eta0 + etar


def dilute_gas_transport(fluid: str, temperature_k: float) -> tuple[float, float]:
    """Return the zero-density viscosity and thermal conductivity.

    Parameters
    ----------
    fluid : str
        One of :data:`TRANSPORT_FLUID_KEYS`.
    temperature_k : float
        Temperature inside the conductivity range of the fluid.

    Returns
    -------
    tuple of (float or None, float)
        ``(eta_0, lambda_0)`` in Pa s and W/(m K).

    Raises
    ------
    InputError
        If the fluid is unknown.
    DomainError
        If the temperature is outside the range.

    References
    ----------
    .. [1] LEMMON-JACOBSEN-2004-TRANSPORT. Lemmon and Jacobsen, Int. J.
           Thermophys. 25, 21 (2004), Eqs. (2) and (5),
           https://trc.nist.gov/refprop/FAQ/NAO.PDF
    .. [2] ASSAEL-2011-H2-CONDUCTIVITY. Assael et al., J. Phys. Chem. Ref.
           Data 40, 033101 (2011), Eq. (2),
           https://srd.nist.gov/jpcrdreprint/1.3606499.pdf
    .. [3] FRIEND-ELY-INGHAM-1989-CH4. Friend, Ely and Ingham, J. Phys. Chem.
           Ref. Data 18, 583 (1989), Eqs. (10a), (12)-(14),
           https://srd.nist.gov/jpcrdreprint/1.555828.pdf
    .. [4] COOLPROP-TRANSPORT-C00300E. CoolProp commit c00300e, hydrogen
           viscosity (Muzny et al. 2013) and helium conductivity (Hands and
           Arp 1981) transcriptions; primary papers not opened,
           https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp
    .. [5] ARP-MCCARTY-FRIEND-1998-HE. Arp, McCarty and Friend, NIST TN 1334
           (1998), Eqs. (1)-(9),
           https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf
    """

    key = _key(fluid)
    t = _positive(temperature_k, "temperature_k")
    v = _VALIDITY[key]
    if not v.conductivity_t_min_k <= t <= v.conductivity_t_max_k:
        raise DomainError(f"{key}: T = {t} K is outside [{v.conductivity_t_min_k}, "
                          f"{v.conductivity_t_max_k}] K.", field="temperature_k")
    if key in ("nitrogen", "oxygen"):
        c = _LJ_N2 if key == "nitrogen" else _LJ_O2
        eta0, _ = _lj_viscosity_parts(c, t, 0.0)
        lam0 = c.n1 * eta0 + c.n2 * (c.t_crit_k / t) ** c.t2 + c.n3 * (c.t_crit_k / t) ** c.t3
        return eta0 * _MICRO, lam0 * _MILLI
    if key == "methane":
        m = _fluid(key).model
        tau = m.t_reduce_k / t
        a0_tt = m.ideal.evaluate(tau, 1.0)[2]
        return (_fei_dilute_viscosity(t) * _MICRO,
                _fei_dilute_conductivity(t, -tau * tau * a0_tt) * _MILLI)
    if key == "helium":
        return _helium_viscosity(t, 0.0), _helium_dilute_conductivity(t)
    c2 = _ASSAEL_PARA if key == "parahydrogen" else _ASSAEL_NORMAL
    return _muzny_viscosity_parts(t, 0.0)[0], _assael_dilute(c2, t)


__all__ = [
    "TRANSPORT_FLUID_KEYS",
    "OlchowySengersParameters",
    "TransportProperties",
    "TransportValidity",
    "dilute_gas_transport",
    "thermal_conductivity_rho_t",
    "transport_from_state",
    "transport_rho_t",
    "transport_tp",
    "transport_validity",
    "viscosity_rho_t",
]
