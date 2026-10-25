"""Reference Helmholtz equations of state for rocket propellants and pressurants.

Six fluids are provided, each with the coefficients of the primary
publication [1]_ - [5]_:

========================  ===================================  ================
key                       equation of state                    validity
========================  ===================================  ================
``parahydrogen``          Leachman et al. 2009 [1]_            13.8033-1000 K,
                                                               <= 2000 MPa
``normal_hydrogen``       Leachman et al. 2009 [1]_            13.957-1000 K,
                                                               <= 2000 MPa
``oxygen``                Schmidt and Wagner 1985 [6]_ as      54.361-300 K,
                          printed by Stewart et al. 1991 [2]_  <= 80 MPa
``nitrogen``              Span et al. 2000 [3]_                63.151-1000 K,
                                                               <= 2200 MPa
``methane``               Setzmann and Wagner 1991 [4]_        90.6941-625 K,
                                                               <= 1000 MPa
``helium``                Ortiz Vega et al. 2023 [5]_          2.1768-1500 K,
                                                               <= 2000 MPa
========================  ===================================  ================

Provenance. The coefficient tables were read from the publisher reprints
hosted by NIST (``srd.nist.gov/jpcrdreprint``) and from NIST IR 8474 [5]_.
For the born-digital papers (hydrogen, nitrogen, helium) the numbers were
extracted from the PDF text layer; for the scanned papers (oxygen, methane)
the OCR layer and the rendered page images were compared. Every residual
coefficient with three or more significant digits was additionally checked
digit-for-digit against the open machine-readable transcription in the
CoolProp fluid files at commit ``c00300e`` [7]_; all of them matched. The
numerical agreement with NIST WebBook [8]_ tables, which are computed from
these same equations, is tested in ``tests/test_real_fluids.py``.

Ideal-gas parts are implemented exactly as printed: hydrogen Eq. (31) with
Table 4 [1]_; oxygen as the Helmholtz integral of the heat-capacity Eq. (11)
[2]_ referenced to the ideal-gas enthalpy 8682 J/mol and entropy
205.037 J/(mol K) at 298.15 K and 0.101325 MPa [2]_ (p. 919); nitrogen
Eq. (53) [3]_; methane Eq. (5.2) with Table 33 [4]_; helium Eq. (3) [5]_.

The hydrogen molar mass (2.01588 g/mol) is not printed in [1]_; it is taken
from the NIST WebBook species record [9]_. Acentric factors, used only by the
cubic fallbacks, are the values displayed on the WebBook fluid pages [8]_.

References
----------
.. [1] LEACHMAN-2009-H2. J. W. Leachman, R. T Jacobsen, S. G. Penoncello and
       E. W. Lemmon, J. Phys. Chem. Ref. Data 38(3), 721-748, 2009: R on
       p. 724, Tables 4-7 on pp. 729-731, Table 8 and Eq. (33) on
       pp. 732-733, https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
.. [2] STEWART-1991-O2. R. B. Stewart, R. T Jacobsen and W. Wagner,
       J. Phys. Chem. Ref. Data 20(5), 917-1021, 1991: constants p. 919,
       Eq. (7) pp. 924-925, Eq. (11) p. 926, Table 4 and Eq. (15) p. 928,
       https://srd.nist.gov/jpcrdreprint/1.555897.pdf
.. [3] SPAN-2000-N2. R. Span et al., J. Phys. Chem. Ref. Data 29(6),
       1361-1433, 2000: constants p. 1364, Eq. (4) p. 1368, Eq. (53)
       p. 1382, Tables 17-18 p. 1385,
       https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
.. [4] SETZMANN-WAGNER-1991-CH4. U. Setzmann and W. Wagner, J. Phys. Chem.
       Ref. Data 20(6), 1061-1155, 1991: constants p. 1064, Eq. (3.2)
       p. 1082, Table 26 p. 1086, Eq. (5.2) and Table 33 p. 1091, Table 35
       p. 1092, https://srd.nist.gov/jpcrdreprint/1.555898.pdf
.. [5] ORTIZVEGA-2023-HE. D. O. Ortiz Vega et al., "An Equation of State for
       the Thermodynamic Properties of Helium", NIST IR 8474, 2023: Table 1
       p. 2, Eq. (3) p. 3, Table 2 p. 4, validity p. 15,
       https://doi.org/10.6028/NIST.IR.8474
.. [6] SCHMIDT-WAGNER-1985-O2. R. Schmidt and W. Wagner, "A New Form of the
       Equation of State for Pure Substances and its Application to Oxygen",
       Fluid Phase Equilib. 19, 175-200, 1985,
       https://doi.org/10.1016/0378-3812(85)87016-3
.. [7] COOLPROP-FLUIDS-C00300E. CoolProp fluid JSON files, commit
       c00300e97ad4cdb30032635d7379209d5e45a121,
       https://github.com/CoolProp/CoolProp/tree/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids
.. [8] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
       "Thermophysical Properties of Fluid Systems", NIST Chemistry WebBook,
       SRD 69, https://webbook.nist.gov/chemistry/fluid/
.. [9] NIST-WEBBOOK-H2-SPECIES. NIST Chemistry WebBook, "Hydrogen" species
       record (molecular weight 2.01588),
       https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from rocket_propulsion.core.constants import USSA1976_SEA_LEVEL_PRESSURE_PA
from rocket_propulsion.core.errors import InputError

from .helmholtz import (
    FluidModel,
    IdealHelmholtz,
    LiquidDensityAncillary,
    RealFluid,
    ResidualHelmholtz,
    ResidualTerm,
    VaporPressureAncillary,
    with_reference_state,
)

_MPA = 1e6


def _pow(n: float, t: float, d: float) -> ResidualTerm:
    return ResidualTerm("power", n, t, d)


def _exp(n: float, t: float, d: float, el: float) -> ResidualTerm:
    return ResidualTerm("exponential", n, t, d, l=el)


def _gauss(n: float, t: float, d: float, eta: float, beta: float, gamma: float,
           epsilon: float) -> ResidualTerm:
    return ResidualTerm("gaussian", n, t, d, eta=eta, beta=beta, gamma=gamma, epsilon=epsilon)


@dataclass(frozen=True, slots=True)
class ReferenceFluidInfo:
    """Constants of a reference fluid needed outside its equation of state.

    Attributes
    ----------
    key : str
        Registry key.
    acentric_factor : float
        Pitzer acentric factor shown on the NIST WebBook fluid page [1]_;
        used only by the cubic fallbacks.
    rho_crit_mol_m3 : float
        Critical molar density of the reference equation.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. NIST Chemistry WebBook SRD 69, fluid pages,
           https://webbook.nist.gov/chemistry/fluid/
    """

    key: str
    acentric_factor: float
    rho_crit_mol_m3: float


# --------------------------------------------------------------------------- hydrogen
# Leachman et al. (2009) Table 5 (N, t, d, p) and Table 6 (phi, beta, gamma, D);
# the Gaussian terms are stored as eta = -phi and beta = -beta (Eq. 32 sign convention).
_PARA_H2_TERMS = (
    _pow(-7.33375, 0.6855, 1),
    _pow(0.01, 1.0, 4),
    _pow(2.60375, 1.0, 1),
    _pow(4.66279, 0.489, 1),
    _pow(0.68239, 0.774, 2),
    _pow(-1.47078, 1.133, 2),
    _pow(0.135801, 1.386, 3),
    _exp(-1.05327, 1.619, 1, 1),
    _exp(0.328239, 1.162, 3, 1),
    _gauss(-0.0577833, 3.96, 2, 1.7437, 0.194, 0.8048, 1.5487),
    _gauss(0.0449743, 5.276, 1, 0.5516, 0.2019, 1.5248, 0.1785),
    _gauss(0.0703464, 0.99, 3, 0.0634, 0.0301, 0.6648, 1.28),
    _gauss(-0.0401766, 6.791, 1, 2.1341, 0.2383, 0.6832, 0.6319),
    _gauss(0.11951, 3.19, 1, 1.777, 0.3253, 1.493, 1.7104),
)
_NORMAL_H2_TERMS = (
    _pow(-6.93643, 0.6844, 1),
    _pow(0.01, 1.0, 4),
    _pow(2.1101, 0.989, 1),
    _pow(4.52059, 0.489, 1),
    _pow(0.732564, 0.803, 2),
    _pow(-1.34086, 1.1444, 2),
    _pow(0.130985, 1.409, 3),
    _exp(-0.777414, 1.754, 1, 1),
    _exp(0.351944, 1.311, 3, 1),
    _gauss(-0.0211716, 4.187, 2, 1.685, 0.171, 0.7164, 1.506),
    _gauss(0.0226312, 5.646, 1, 0.489, 0.2245, 1.3444, 0.156),
    _gauss(0.032187, 0.791, 3, 0.103, 0.1304, 1.4517, 1.736),
    _gauss(-0.0231752, 7.249, 1, 2.506, 0.2785, 0.7204, 0.67),
    _gauss(0.0557346, 2.986, 1, 1.607, 0.3967, 1.5445, 1.662),
)
# Table 4: a_k ln(1 - exp(b_k tau)) is stored as (a_k, theta_k = -b_k).
_PARA_H2_IDEAL = IdealHelmholtz(
    a1=-1.4485891134, a2=1.884521239, log_tau=1.5,
    planck_einstein=((4.30256, 15.1496751472), (13.0289, 25.0925982148),
                     (-47.7365, 29.4735563787), (50.0013, 35.4059141417),
                     (-18.6261, 40.724998482), (0.993973, 163.7925799988),
                     (0.536078, 309.2173173842)))
_NORMAL_H2_IDEAL = IdealHelmholtz(
    a1=-1.4579856475, a2=1.888076782, log_tau=1.5,
    planck_einstein=((1.616, 16.0205159149), (-0.4117, 22.6580178006),
                     (-0.792, 60.0090511389), (0.758, 74.9434303817),
                     (1.217, 206.9392065168)))
_H2_GAS_CONSTANT = 8.314472  # Leachman et al. (2009), nomenclature, p. 724
_H2_MOLAR_MASS = 2.01588e-3  # NIST WebBook hydrogen species record

# --------------------------------------------------------------------------- oxygen
# Stewart, Jacobsen and Wagner (1991) Table 4 (k, i = d, l, j = t, N_k), Eq. (15);
# gamma = 1 for l > 0.
_O2_TERMS = (
    _pow(0.3983768749, 0.0, 1),
    _pow(-1.846157454, 1.5, 1),
    _pow(0.4183473197, 2.5, 1),
    _pow(0.02370620711, -0.5, 2),
    _pow(0.09771730573, 1.5, 2),
    _pow(0.03017891294, 2.0, 2),
    _pow(0.02273353212, 0.0, 3),
    _pow(0.01357254086, 1.0, 3),
    _pow(-0.04052698943, 2.5, 3),
    _pow(0.0005454628515, 0.0, 6),
    _pow(0.0005113182277, 2.0, 7),
    _pow(2.953466883e-07, 5.0, 7),
    _pow(-8.687645072e-05, 2.0, 8),
    _exp(-0.2127082589, 5.0, 1, 2),
    _exp(0.08735941958, 6.0, 1, 2),
    _exp(0.127550919, 3.5, 2, 2),
    _exp(-0.09067701064, 5.5, 2, 2),
    _exp(-0.03540084206, 3.0, 3, 2),
    _exp(-0.03623278059, 7.0, 3, 2),
    _exp(0.0132769929, 6.0, 5, 2),
    _exp(-0.0003254111865, 8.5, 6, 2),
    _exp(-0.008313582932, 4.0, 7, 2),
    _exp(0.002124570559, 6.5, 8, 2),
    _exp(-0.0008325206232, 5.5, 10, 2),
    _exp(-2.626173276e-05, 22.0, 2, 4),
    _exp(0.002599581482, 11.0, 3, 4),
    _exp(0.009984649663, 18.0, 3, 4),
    _exp(0.002199923153, 11.0, 4, 4),
    _exp(-0.02591350486, 23.0, 4, 4),
    _exp(-0.1259630848, 17.0, 5, 4),
    _exp(0.1478355637, 18.0, 5, 4),
    _exp(-0.01011251078, 23.0, 5, 4),
)
_O2_TC = 154.581
_O2_RHOC = 13.63e3
_O2_GAS_CONSTANT = 8.31434


def _oxygen_ideal() -> IdealHelmholtz:
    """Helmholtz integral of Stewart et al. (1991) Eq. (11).

    ``c_p0/R = N1 T^-1.5 + N2 + N3 T^2 + N4 u^2 e^u/(e^u - 1)^2
    + N5 (2/3) eta^2 e^-eta / (1 + (2/3) e^-eta)^2`` with ``u = N6/T`` and
    ``eta = N7/T``. A term ``c T^k`` gives ``-c T_c^k/(k (k + 1)) tau^-k``.
    """

    n1, n2, n3, n4, n5, n6, n7 = (1.06778, 3.50042, 0.166961e-7, 1.01258, 0.944365,
                                  2242.45, 11580.4)
    shape = IdealHelmholtz(
        a1=0.0, a2=0.0, log_tau=n2 - 1.0,
        power=((-n1 * _O2_TC ** -1.5 / ((-1.5) * (-0.5)), 1.5),
               (-n3 * _O2_TC ** 2 / (2.0 * 3.0), -2.0)),
        planck_einstein=((n4, n6 / _O2_TC),),
        two_level=((-n5, 2.0 / 3.0, n7 / _O2_TC),))
    # Reference state of Stewart et al. (1991), p. 919: ideal gas at T0 = 298.15 K and
    # P0 = 0.101325 MPa with H0 = 8682 J/mol and S0 = 205.037 J/(mol K).
    return with_reference_state(
        shape, gas_constant_j_mol_k=_O2_GAS_CONSTANT, t_reduce_k=_O2_TC,
        rho_reduce_mol_m3=_O2_RHOC, temperature_k=298.15,
        pressure_pa=USSA1976_SEA_LEVEL_PRESSURE_PA, enthalpy_j_mol=8682.0,
        entropy_j_mol_k=205.037)


# --------------------------------------------------------------------------- nitrogen
# Span et al. (2000) Table 17 (N_k, i_k = d, j_k = t, l_k) and Table 18 (phi, beta,
# gamma; epsilon = 1 in Eq. 55).
_N2_TERMS = (
    _pow(0.924803575275, 0.25, 1),
    _pow(-0.492448489428, 0.875, 1),
    _pow(0.661883336938, 0.5, 2),
    _pow(-1.92902649201, 0.875, 2),
    _pow(-0.0622469309629, 0.375, 3),
    _pow(0.349943957581, 0.75, 3),
    _exp(0.564857472498, 0.5, 1, 1),
    _exp(-1.61720005987, 0.75, 1, 1),
    _exp(-0.481395031883, 2.0, 1, 1),
    _exp(0.421150636384, 1.25, 3, 1),
    _exp(-0.0161962230825, 3.5, 3, 1),
    _exp(0.172100994165, 1.0, 4, 1),
    _exp(0.00735448924933, 0.5, 6, 1),
    _exp(0.0168077305479, 3.0, 6, 1),
    _exp(-0.00107626664179, 0.0, 7, 1),
    _exp(-0.0137318088513, 2.75, 7, 1),
    _exp(0.000635466899859, 0.75, 8, 1),
    _exp(0.00304432279419, 2.5, 8, 1),
    _exp(-0.0435762336045, 4.0, 1, 2),
    _exp(-0.0723174889316, 6.0, 2, 2),
    _exp(0.0389644315272, 6.0, 3, 2),
    _exp(-0.021220136391, 3.0, 4, 2),
    _exp(0.00408822981509, 3.0, 5, 2),
    _exp(-5.51990017984e-05, 6.0, 8, 2),
    _exp(-0.0462016716479, 16.0, 4, 3),
    _exp(-0.00300311716011, 11.0, 5, 3),
    _exp(0.0368825891208, 15.0, 5, 3),
    _exp(-0.0025585684622, 12.0, 8, 3),
    _exp(0.00896915264558, 12.0, 3, 4),
    _exp(-0.0044151337035, 7.0, 5, 4),
    _exp(0.00133722924858, 4.0, 6, 4),
    _exp(0.000264832491957, 16.0, 9, 4),
    _gauss(19.6688194015, 0.0, 1, 20.0, 325.0, 1.16, 1.0),
    _gauss(-20.911560073, 1.0, 1, 20.0, 325.0, 1.16, 1.0),
    _gauss(0.0167788306989, 2.0, 3, 15.0, 300.0, 1.13, 1.0),
    _gauss(2627.67566274, 3.0, 2, 25.0, 275.0, 1.25, 1.0),
)
# Eq. (53): ln(delta) + a1 ln(tau) + a2 + a3 tau + a4/tau + a5/tau^2 + a6/tau^3
#           + a7 ln(1 - exp(-a8 tau)).
_N2_IDEAL = IdealHelmholtz(
    a1=-12.76952708, a2=-0.00784163, log_tau=2.5,
    power=((-1.934819e-4, -1.0), (-1.247742e-5, -2.0), (6.678326e-8, -3.0)),
    planck_einstein=((1.012941, 26.65788),))

# --------------------------------------------------------------------------- methane
# Setzmann and Wagner (1991) Table 35 (n_i, c_i = l, d_i, t_i; alpha, beta, gamma, Delta).
_CH4_TERMS = (
    _pow(0.04367901028, -0.5, 1),
    _pow(0.6709236199, 0.5, 1),
    _pow(-1.765577859, 1.0, 1),
    _pow(0.8582330241, 0.5, 2),
    _pow(-1.206513052, 1.0, 2),
    _pow(0.512046722, 1.5, 2),
    _pow(-0.0004000010791, 4.5, 2),
    _pow(-0.01247842423, 0.0, 3),
    _pow(0.03100269701, 1.0, 4),
    _pow(0.001754748522, 3.0, 4),
    _pow(-3.171921605e-06, 1.0, 8),
    _pow(-2.24034684e-06, 3.0, 9),
    _pow(2.947056156e-07, 3.0, 10),
    _exp(0.1830487909, 0.0, 1, 1),
    _exp(0.1511883679, 1.0, 1, 1),
    _exp(-0.4289363877, 2.0, 1, 1),
    _exp(0.06894002446, 0.0, 2, 1),
    _exp(-0.01408313996, 0.0, 4, 1),
    _exp(-0.0306305483, 2.0, 5, 1),
    _exp(-0.02969906708, 2.0, 6, 1),
    _exp(-0.01932040831, 5.0, 1, 2),
    _exp(-0.1105739959, 5.0, 2, 2),
    _exp(0.09952548995, 5.0, 3, 2),
    _exp(0.008548437825, 2.0, 4, 2),
    _exp(-0.06150555662, 4.0, 4, 2),
    _exp(-0.04291792423, 12.0, 3, 3),
    _exp(-0.0181320729, 8.0, 5, 3),
    _exp(0.0344590476, 10.0, 5, 3),
    _exp(-0.00238591945, 10.0, 8, 3),
    _exp(-0.01159094939, 10.0, 2, 4),
    _exp(0.06641693602, 14.0, 3, 4),
    _exp(-0.0237154959, 12.0, 4, 4),
    _exp(-0.03961624905, 18.0, 4, 4),
    _exp(-0.01387292044, 22.0, 4, 4),
    _exp(0.03389489599, 18.0, 5, 4),
    _exp(-0.002927378753, 14.0, 6, 4),
    _gauss(9.324799946e-05, 2.0, 2, 20.0, 200.0, 1.07, 1.0),
    _gauss(-6.287171518, 0.0, 0, 40.0, 250.0, 1.11, 1.0),
    _gauss(12.71069467, 1.0, 0, 40.0, 250.0, 1.11, 1.0),
    _gauss(-6.423953466, 2.0, 0, 40.0, 250.0, 1.11, 1.0),
)
# Table 33 / Eq. (5.2): a1 + a2 tau + a3 ln(tau) + sum a_i ln(1 - exp(-theta_i tau)).
_CH4_IDEAL = IdealHelmholtz(
    a1=9.91243972, a2=-6.33270087, log_tau=3.0016,
    planck_einstein=((0.008449, 3.40043240), (4.6942, 10.26951575), (3.4865, 20.43932747),
                     (1.6572, 29.93744884), (1.4115, 79.13351945)))
_CH4_MOLAR_MASS = 16.0428e-3
_CH4_RHOC_MASS = 162.66

# --------------------------------------------------------------------------- helium
# Ortiz Vega et al. (2023) Table 2 (n, t, d, l; eta, beta, gamma, epsilon).
_HE_TERMS = (
    _pow(0.015559018, 1.0, 4),
    _pow(3.0638932, 0.425, 1),
    _pow(-4.2420844, 0.63, 1),
    _pow(0.054418088, 0.69, 2),
    _pow(-0.18971904, 1.83, 2),
    _pow(0.087856262, 0.575, 3),
    _exp(2.2833566, 0.925, 1, 1),
    _exp(-0.53331595, 1.585, 1, 2),
    _exp(-0.53296502, 1.69, 3, 2),
    _exp(0.99444915, 1.51, 2, 1),
    _exp(-0.30078896, 2.9, 2, 2),
    _exp(-1.6432563, 0.8, 1, 1),
    _gauss(0.8029102, 1.26, 2, 1.5497, 0.2471, 3.15, 0.596),
    _gauss(0.026838669, 3.51, 1, 9.245, 0.0983, 2.54505, 0.3423),
    _gauss(0.04687678, 2.785, 2, 4.76323, 0.1556, 1.2513, 0.761),
    _gauss(-0.14832766, 1.0, 1, 6.3826, 2.6782, 1.9416, 0.9747),
    _gauss(0.03016211, 4.22, 1, 8.7023, 2.7077, 0.5984, 0.5868),
    _gauss(-0.019986041, 0.83, 3, 0.255, 0.6621, 2.2282, 0.5627),
    _gauss(0.14283514, 1.575, 2, 0.3523, 0.1775, 1.606, 2.5346),
    _gauss(0.007418269, 3.447, 2, 0.1492, 0.4821, 3.815, 3.6763),
    _gauss(-0.22989793, 0.73, 3, 0.05, 0.3069, 1.61958, 4.5245),
    _gauss(0.79224829, 1.634, 2, 0.1668, 0.1758, 0.6407, 5.039),
    _gauss(-0.049386338, 6.13, 2, 42.2358, 1357.6577, 1.076, 0.959),
)
_HE_IDEAL = IdealHelmholtz(a1=0.1733487932835764, a2=0.4674522201550815, log_tau=1.5)

_PROVENANCE = (("Coefficients read from the primary source and cross-checked digit-for-digit "
               "against CoolProp fluid JSON at commit c00300e."),)


def _para_hydrogen() -> FluidModel:
    tc, pc = 32.938, 1.2858 * _MPA
    return FluidModel(
        name="parahydrogen", key="parahydrogen", model_kind="reference",
        molar_mass_kg_mol=_H2_MOLAR_MASS, gas_constant_j_mol_k=_H2_GAS_CONSTANT,
        t_reduce_k=tc, rho_reduce_mol_m3=15.538e3, t_crit_k=tc, p_crit_pa=pc,
        t_min_k=13.8033, t_max_k=1000.0, p_max_pa=2000.0 * _MPA,
        ideal=_PARA_H2_IDEAL, residual=ResidualHelmholtz(_PARA_H2_TERMS),
        ancillary=VaporPressureAncillary(
            tc, pc, ((-4.87767, 1.0), (1.03359, 1.5), (0.82668, 2.65), (-0.129412, 7.4)),
            "LEACHMAN-2009-H2 Table 8, Eq. (33)"),
        references=("LEACHMAN-2009-H2 Tables 4-7, Eqs. (31)-(32)", "NIST-WEBBOOK-H2-SPECIES"),
        notes=_PROVENANCE)


def _normal_hydrogen() -> FluidModel:
    tc, pc = 33.145, 1.2964 * _MPA
    return FluidModel(
        name="normal hydrogen", key="normal_hydrogen", model_kind="reference",
        molar_mass_kg_mol=_H2_MOLAR_MASS, gas_constant_j_mol_k=_H2_GAS_CONSTANT,
        t_reduce_k=tc, rho_reduce_mol_m3=15.508e3, t_crit_k=tc, p_crit_pa=pc,
        t_min_k=13.957, t_max_k=1000.0, p_max_pa=2000.0 * _MPA,
        ideal=_NORMAL_H2_IDEAL, residual=ResidualHelmholtz(_NORMAL_H2_TERMS),
        ancillary=VaporPressureAncillary(
            tc, pc, ((-4.89789, 1.0), (0.988558, 1.5), (0.349689, 2.0), (0.499356, 2.85)),
            "LEACHMAN-2009-H2 Table 8, Eq. (33)"),
        references=("LEACHMAN-2009-H2 Tables 4-7, Eqs. (31)-(32)", "NIST-WEBBOOK-H2-SPECIES"),
        notes=_PROVENANCE)


def _oxygen() -> FluidModel:
    return FluidModel(
        name="oxygen", key="oxygen", model_kind="reference",
        molar_mass_kg_mol=31.9988e-3, gas_constant_j_mol_k=_O2_GAS_CONSTANT,
        t_reduce_k=_O2_TC, rho_reduce_mol_m3=_O2_RHOC, t_crit_k=_O2_TC,
        p_crit_pa=5.043 * _MPA, t_min_k=54.361, t_max_k=300.0, p_max_pa=80.0 * _MPA,
        ideal=_oxygen_ideal(), residual=ResidualHelmholtz(_O2_TERMS),
        ancillary=VaporPressureAncillary(
            _O2_TC, 5.0430 * _MPA,
            ((-6.043938, 1.0), (1.175627, 1.5), (-0.994086, 3.0), (-3.456781, 7.0),
             (3.361499, 9.0)), "STEWART-1991-O2 Eq. (7)"),
        liquid_ancillary=LiquidDensityAncillary(
            _O2_TC, _O2_RHOC, ((1.507678, 1.0 / 3.0), (0.85810805, 2.0 / 3.0),
                               (0.19035504, 3.0)), "linear", "STEWART-1991-O2 Eq. (10)"),
        references=("SCHMIDT-WAGNER-1985-O2", "STEWART-1991-O2 Table 4, Eqs. (11), (15)"),
        notes=_PROVENANCE + (
            ("Validity is the range of the Stewart et al. (1991) title: triple point to 300 K, "
            "pressures to 80 MPa."),))


def _nitrogen() -> FluidModel:
    tc, pc = 126.192, 3.3958 * _MPA
    return FluidModel(
        name="nitrogen", key="nitrogen", model_kind="reference",
        molar_mass_kg_mol=28.01348e-3, gas_constant_j_mol_k=8.314510,
        t_reduce_k=tc, rho_reduce_mol_m3=11.1839e3, t_crit_k=tc, p_crit_pa=pc,
        t_min_k=63.151, t_max_k=1000.0, p_max_pa=2200.0 * _MPA,
        ideal=_N2_IDEAL, residual=ResidualHelmholtz(_N2_TERMS),
        ancillary=VaporPressureAncillary(
            tc, pc, ((-6.12445284, 1.0), (1.26327220, 1.5), (-0.765910082, 2.5),
                     (-1.77570564, 5.0)), "SPAN-2000-N2 Eq. (4)"),
        liquid_ancillary=LiquidDensityAncillary(
            tc, 11.1839e3, ((1.48654237, 0.3294), (-0.280476066, 2.0 / 3.0),
                            (0.0894143085, 8.0 / 3.0), (-0.119879866, 35.0 / 6.0)), "log",
            "SPAN-2000-N2 Eq. (5)"),
        references=("SPAN-2000-N2 Tables 17-18, Eqs. (53), (55)",), notes=_PROVENANCE)


def _methane() -> FluidModel:
    tc, pc = 190.564, 4.5992 * _MPA
    return FluidModel(
        name="methane", key="methane", model_kind="reference",
        molar_mass_kg_mol=_CH4_MOLAR_MASS, gas_constant_j_mol_k=8.31451,
        t_reduce_k=tc, rho_reduce_mol_m3=_CH4_RHOC_MASS / _CH4_MOLAR_MASS, t_crit_k=tc,
        p_crit_pa=pc, t_min_k=90.6941, t_max_k=625.0, p_max_pa=1000.0 * _MPA,
        ideal=_CH4_IDEAL, residual=ResidualHelmholtz(_CH4_TERMS),
        ancillary=VaporPressureAncillary(
            tc, pc, ((-6.036219, 1.0), (1.409353, 1.5), (-0.4945199, 2.0), (-1.443048, 4.5)),
            "SETZMANN-WAGNER-1991-CH4 Eq. (3.2), Table 26"),
        liquid_ancillary=LiquidDensityAncillary(
            tc, _CH4_RHOC_MASS / _CH4_MOLAR_MASS,
            ((1.9906389, 0.354), (-0.78756197, 0.5), (0.036976723, 2.5)), "log",
            "SETZMANN-WAGNER-1991-CH4 Eq. (3.4), Table 26"),
        references=("SETZMANN-WAGNER-1991-CH4 Tables 33, 35, Eqs. (5.1)-(5.3)",),
        notes=_PROVENANCE + (
            ("The paper's nomenclature (p. 1064) prints p_c = 4.5922 MPa; Eq. (2.3), Table 26 "
            "and p. 1091 give 4.5992 MPa, which is used."),
            ("Validity is 90.6941 K <= T < 625 K, p <= 1000 MPa (p. 1091); the melting "
            "line is not checked.")))


def _helium() -> FluidModel:
    tc = 5.1953
    return FluidModel(
        name="helium", key="helium", model_kind="reference",
        molar_mass_kg_mol=4.002602e-3, gas_constant_j_mol_k=8.314472,
        t_reduce_k=tc, rho_reduce_mol_m3=17.3837e3, t_crit_k=tc, p_crit_pa=0.22832 * _MPA,
        t_min_k=2.1768, t_max_k=1500.0, p_max_pa=2000.0 * _MPA,
        ideal=_HE_IDEAL, residual=ResidualHelmholtz(_HE_TERMS), ancillary=None,
        references=("ORTIZVEGA-2023-HE Tables 1-2, Eqs. (3)-(4)",),
        notes=_PROVENANCE + (
            ("NIST IR 8474 describes this as a preliminary equation; no vapor-pressure "
            "ancillary is given, so saturation starts from the spinodal loop."),
            ("Valid for helium I only (above the lambda line); the lambda line and melting "
            "line are not checked.")))


_BUILDERS = {
    "parahydrogen": _para_hydrogen,
    "normal_hydrogen": _normal_hydrogen,
    "oxygen": _oxygen,
    "nitrogen": _nitrogen,
    "methane": _methane,
    "helium": _helium,
}

_INFO = {
    "parahydrogen": ReferenceFluidInfo("parahydrogen", -0.219, 15.538e3),
    "normal_hydrogen": ReferenceFluidInfo("normal_hydrogen", -0.219, 15.508e3),
    "oxygen": ReferenceFluidInfo("oxygen", 0.0222, _O2_RHOC),
    "nitrogen": ReferenceFluidInfo("nitrogen", 0.0372, 11.1839e3),
    "methane": ReferenceFluidInfo("methane", 0.01142, _CH4_RHOC_MASS / _CH4_MOLAR_MASS),
    "helium": ReferenceFluidInfo("helium", -0.3836, 17.3837e3),
}

REFERENCE_FLUID_KEYS: tuple[str, ...] = tuple(_BUILDERS)


def _key(name: object) -> str:
    if not isinstance(name, str) or name not in _BUILDERS:
        raise InputError(
            f"Unknown reference fluid {name!r}; choose one of {', '.join(_BUILDERS)}.",
            field="fluid")
    return name


@cache
def reference_fluid_model(name: str) -> FluidModel:
    """Return the :class:`FluidModel` of a reference fluid.

    Parameters
    ----------
    name : str
        One of :data:`REFERENCE_FLUID_KEYS`.

    Returns
    -------
    FluidModel
        Equation of state with its declared validity range.

    Raises
    ------
    InputError
        If the name is unknown.

    References
    ----------
    .. [1] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    .. [2] STEWART-1991-O2. Stewart et al., J. Phys. Chem. Ref. Data 20, 917
           (1991), https://srd.nist.gov/jpcrdreprint/1.555897.pdf
    .. [3] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    .. [4] SETZMANN-WAGNER-1991-CH4. Setzmann and Wagner, J. Phys. Chem. Ref.
           Data 20, 1061 (1991),
           https://srd.nist.gov/jpcrdreprint/1.555898.pdf
    .. [5] ORTIZVEGA-2023-HE. Ortiz Vega et al., NIST IR 8474 (2023),
           https://doi.org/10.6028/NIST.IR.8474
    """

    return _BUILDERS[_key(name)]()


def reference_fluid(name: str) -> RealFluid:
    """Return a :class:`RealFluid` calculator for a reference fluid.

    Parameters
    ----------
    name : str
        One of :data:`REFERENCE_FLUID_KEYS`.

    Returns
    -------
    RealFluid
        New calculator (with its own saturation cache).

    Raises
    ------
    InputError
        If the name is unknown.

    References
    ----------
    .. [1] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361
           (2000), https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
    """

    return RealFluid(reference_fluid_model(name))


def reference_fluid_info(name: str) -> ReferenceFluidInfo:
    """Return the acentric factor and critical density of a reference fluid.

    Parameters
    ----------
    name : str
        One of :data:`REFERENCE_FLUID_KEYS`.

    Returns
    -------
    ReferenceFluidInfo
        Constants used by the cubic fallbacks.

    Raises
    ------
    InputError
        If the name is unknown.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. NIST Chemistry WebBook SRD 69, fluid pages
           (acentric factor), https://webbook.nist.gov/chemistry/fluid/
    """

    return _INFO[_key(name)]
