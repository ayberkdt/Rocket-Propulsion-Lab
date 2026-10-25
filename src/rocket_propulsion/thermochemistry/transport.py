"""Gas viscosity, thermal conductivity and Prandtl number, NASA CEA formulation.

Species properties come from the least-squares fits of NASA's CEA transport
database (``data/nasa_cea_transport.json``, extracted by
``scripts/sync_nasa_transport.py`` from ``data/trans.inp`` at a pinned
commit [2]_), evaluated with RP-1311 Part I equation (5.1) [1]_. Mixture
viscosity and frozen thermal conductivity follow equations (5.3)--(5.7); the
reaction contributions to thermal conductivity and specific heat of a mixture
in local chemical equilibrium follow equations (5.8)--(5.18), and the Prandtl
numbers equations (5.19a,b). Two implementation constants that RP-1311 states
only in words, the collision-cross-section ratio ``A* = 1.1`` of equation
(5.13) and the rigid-sphere estimate of missing interaction viscosities, were
checked against CEA's own source at the same commit [3]_.

General-purpose fallbacks for gases outside the CEA database are Sutherland's
law in the form of the U.S. Standard Atmosphere 1976 [4]_ and the
Chapman--Enskog viscosity with the Neufeld et al. (1972) fit of the
Lennard-Jones collision integral as reproduced in [5]_.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis", chapter 5.
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
       ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
       https://github.com/nasa/cea
.. [3] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
       ``source/equilibrium.f90`` (``compute_transport_properties``) and
       ``source/transport.f90``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
       https://github.com/nasa/cea
.. [4] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", equation (51).
       https://ntrs.nasa.gov/citations/19770009539
.. [5] IDAES-CHAPMAN-ENSKOG. IDAES-PSE documentation, "Chapman Enskog Viscosity
       Model". https://idaes-pse.readthedocs.io/en/main/explanations/components/property_package/general/transport_properties/chapman_enskog.html
"""

from __future__ import annotations

import itertools
import json
from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache
from math import exp, isfinite, log, sqrt
from pathlib import Path

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import EquilibriumResult
from .species import get_species, species_catalog

DATA_FILE = Path(__file__).parent / "data" / "nasa_cea_transport.json"

#: Unit of the trans.inp viscosity fits (RP-1311 Part I, section 5.1): 1 uP = 1e-7 Pa s.
MICROPOISE_TO_PA_S = 1.0e-7
#: Unit of the trans.inp conductivity fits: 1 uW/(cm K) = 1e-4 W/(m K).
MICROWATT_PER_CM_K_TO_W_PER_M_K = 1.0e-4
#: Collision cross-section ratio substituted for A*_kl in RP-1311 eq. (5.13).
COLLISION_CROSS_SECTION_RATIO = 1.1
#: Default mole fraction below which a species is left out of the mixture rules.
DEFAULT_MINIMUM_MOLE_FRACTION = 1.0e-10
#: Default largest total mole fraction of species without transport data that may be omitted.
DEFAULT_MAXIMUM_OMITTED_MOLE_FRACTION = 1.0e-3


@dataclass(frozen=True, slots=True)
class TransportFitInterval:
    """One temperature interval of an RP-1311 equation (5.1) fit [1]_.

    Attributes
    ----------
    minimum_temperature_k, maximum_temperature_k : float
        Interval bounds in K.
    a, b, c, d : float
        Coefficients of ``ln(y) = a ln T + b/T + c/T**2 + d`` with ``y`` in
        micropoise (viscosity) or microwatt/(cm K) (conductivity).

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equation (5.1). https://ntrs.nasa.gov/citations/19950013764
    """

    minimum_temperature_k: float
    maximum_temperature_k: float
    a: float
    b: float
    c: float
    d: float

    def log_value(self, temperature_k: float) -> float:
        """Return ``ln(y)`` of equation (5.1) in the fit's CGS unit [1]_.

        Parameters
        ----------
        temperature_k : float
            Temperature in K.

        Returns
        -------
        float
            Natural logarithm of the property in micropoise or uW/(cm K).

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis",
               equation (5.1). https://ntrs.nasa.gov/citations/19950013764
        """
        t = temperature_k
        return self.a * log(t) + self.b / t + self.c / (t * t) + self.d


def _select_interval(
    intervals: tuple[TransportFitInterval, ...],
    temperature_k: float,
    *,
    owner: str,
    quantity: str,
    allow_extrapolation: bool,
) -> tuple[TransportFitInterval, bool]:
    """Return ``(interval, extrapolated)`` for ``temperature_k``.

    Inside the fitted range the interval containing ``temperature_k`` is used
    (the upper interval at a shared bound, as CEA does). Outside it the first
    or last interval is used only when ``allow_extrapolation`` is true, which
    reproduces CEA's own selection rule in ``source/transport.f90``
    (``st_calc_eta``); otherwise :class:`DomainError` is raised.
    """
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise InputError("Temperature must be finite and positive.", field="temperature_k")
    lower = intervals[0].minimum_temperature_k
    upper = intervals[-1].maximum_temperature_k
    extrapolated = temperature_k < lower or temperature_k > upper
    if extrapolated and not allow_extrapolation:
        raise DomainError(
            f"{quantity} fit of {owner} covers [{lower:g}, {upper:g}] K.",
            code="transport_temperature_out_of_range",
            field="temperature_k",
            details={"species": owner, "minimum_temperature_k": lower,
                     "maximum_temperature_k": upper},
        )
    selected = intervals[0]
    for previous, interval in itertools.pairwise(intervals):
        if temperature_k > previous.maximum_temperature_k:
            selected = interval
    return selected, extrapolated


@dataclass(frozen=True, slots=True)
class SpeciesTransportData:
    """CEA viscosity and conductivity fits of one gas [1]_ [2]_.

    Attributes
    ----------
    name : str
        trans.inp species name (``"H2O"``, ``"C2H2,acetylene"``).
    viscosity : tuple of TransportFitInterval
        Viscosity fit intervals in increasing temperature order.
    thermal_conductivity : tuple of TransportFitInterval
        Conductivity fit intervals; may be empty for a few species.
    data_source : str
        Per-record source comment of trans.inp (``"SVEHLA (1994)"``).

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           section 5.1. https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
           https://github.com/nasa/cea
    """

    name: str
    viscosity: tuple[TransportFitInterval, ...]
    thermal_conductivity: tuple[TransportFitInterval, ...]
    data_source: str

    def viscosity_pa_s(self, temperature_k: float, *, allow_extrapolation: bool = False) -> float:
        """Dynamic viscosity from equation (5.1) [1]_, in Pa s.

        Parameters
        ----------
        temperature_k : float
            Temperature in K.
        allow_extrapolation : bool, default False
            Evaluate the nearest end interval outside the fitted range, as CEA
            does [2]_, instead of raising.

        Returns
        -------
        float
            Viscosity in Pa s.

        Raises
        ------
        InputError
            If the temperature is not finite and positive.
        DomainError
            If the species has no viscosity fit (``transport_data_missing``) or
            the temperature is outside the fit (``transport_temperature_out_of_range``).

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis",
               equation (5.1). https://ntrs.nasa.gov/citations/19950013764
        .. [2] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
               ``source/transport.f90``, ``st_calc_eta``.
               https://github.com/nasa/cea
        """
        if not self.viscosity:
            raise DomainError(f"{self.name} has no viscosity fit.", code="transport_data_missing")
        interval, _ = _select_interval(
            self.viscosity, temperature_k, owner=self.name, quantity="Viscosity",
            allow_extrapolation=allow_extrapolation,
        )
        return exp(interval.log_value(temperature_k)) * MICROPOISE_TO_PA_S

    def thermal_conductivity_w_m_k(
        self, temperature_k: float, *, allow_extrapolation: bool = False
    ) -> float:
        """Thermal conductivity from equation (5.1) [1]_, in W/(m K).

        Parameters
        ----------
        temperature_k : float
            Temperature in K.
        allow_extrapolation : bool, default False
            Evaluate the nearest end interval outside the fitted range, as CEA
            does [2]_, instead of raising.

        Returns
        -------
        float
            Frozen (translational plus internal) conductivity in W/(m K).

        Raises
        ------
        InputError
            If the temperature is not finite and positive.
        DomainError
            If the species has no conductivity fit (``transport_data_missing``)
            or the temperature is outside the fit.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis",
               equation (5.1). https://ntrs.nasa.gov/citations/19950013764
        .. [2] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
               ``source/transport.f90``, ``st_calc_lambda``.
               https://github.com/nasa/cea
        """
        if not self.thermal_conductivity:
            raise DomainError(
                f"{self.name} has no thermal-conductivity fit.", code="transport_data_missing"
            )
        interval, _ = _select_interval(
            self.thermal_conductivity, temperature_k, owner=self.name,
            quantity="Thermal-conductivity", allow_extrapolation=allow_extrapolation,
        )
        return exp(interval.log_value(temperature_k)) * MICROWATT_PER_CM_K_TO_W_PER_M_K


@dataclass(frozen=True, slots=True)
class BinaryInteractionData:
    """CEA fit of the interaction viscosity ``eta_ij`` of one species pair [1]_ [2]_.

    Attributes
    ----------
    pair : tuple of str
        The two trans.inp species names.
    viscosity : tuple of TransportFitInterval
        Fit intervals of ``eta_ij`` (same form as equation (5.1)).
    data_source : str
        Per-record source comment of trans.inp.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           sections 5.1 and 5.2.1. https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
           https://github.com/nasa/cea
    """

    pair: tuple[str, str]
    viscosity: tuple[TransportFitInterval, ...]
    data_source: str

    def viscosity_pa_s(self, temperature_k: float, *, allow_extrapolation: bool = False) -> float:
        """Interaction viscosity ``eta_ij`` in Pa s [1]_.

        Parameters
        ----------
        temperature_k : float
            Temperature in K.
        allow_extrapolation : bool, default False
            Use the nearest end interval outside the fitted range.

        Returns
        -------
        float
            ``eta_ij`` in Pa s.

        Raises
        ------
        DomainError
            If the temperature is outside the fit and extrapolation is off.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis",
               equations (5.1) and (5.7). https://ntrs.nasa.gov/citations/19950013764
        """
        interval, _ = _select_interval(
            self.viscosity, temperature_k, owner="-".join(self.pair),
            quantity="Interaction-viscosity", allow_extrapolation=allow_extrapolation,
        )
        return exp(interval.log_value(temperature_k)) * MICROPOISE_TO_PA_S


def _intervals(records: list[dict]) -> tuple[TransportFitInterval, ...]:
    """Build fit intervals from bundled JSON records."""
    return tuple(
        TransportFitInterval(
            minimum_temperature_k=float(record["minimum_temperature_k"]),
            maximum_temperature_k=float(record["maximum_temperature_k"]),
            a=float(record["coefficients"][0]),
            b=float(record["coefficients"][1]),
            c=float(record["coefficients"][2]),
            d=float(record["coefficients"][3]),
        )
        for record in records
    )


@lru_cache(maxsize=1)
def _load() -> tuple[dict[str, SpeciesTransportData], dict[frozenset[str], BinaryInteractionData]]:
    """Read the bundled CEA transport data once; fail on duplicate species keys."""
    document = json.loads(DATA_FILE.read_text(encoding="utf-8"))
    species: dict[str, SpeciesTransportData] = {}
    for record in document["species"]:
        key = record["name"].upper()
        if key in species:
            raise DomainError(f"Duplicate transport record {key}.", code="transport_data_invalid")
        species[key] = SpeciesTransportData(
            name=record["name"],
            viscosity=_intervals(record["viscosity"]),
            thermal_conductivity=_intervals(record["thermal_conductivity"]),
            data_source=record["trans_inp_comment"],
        )
    binaries: dict[frozenset[str], BinaryInteractionData] = {}
    for record in document["binary_interactions"]:
        first, second = record["pair"]
        binaries[frozenset((first.upper(), second.upper()))] = BinaryInteractionData(
            pair=(first, second),
            viscosity=_intervals(record["viscosity"]),
            data_source=record["trans_inp_comment"],
        )
    return species, binaries


def species_transport_catalog() -> dict[str, SpeciesTransportData]:
    """Return every bundled species fit keyed by upper-case trans.inp name [1]_.

    Returns
    -------
    dict of str to SpeciesTransportData
        A new dictionary; keys are upper case to match the species catalog of
        :mod:`.species` (``"C2H2,ACETYLENE"``).

    References
    ----------
    .. [1] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
           https://github.com/nasa/cea
    """
    return dict(_load()[0])


def get_species_transport(name: str) -> SpeciesTransportData:
    """Return the transport fits of one species (case-insensitive name) [1]_.

    Parameters
    ----------
    name : str
        trans.inp or species-catalog name.

    Returns
    -------
    SpeciesTransportData
        Bundled fits.

    Raises
    ------
    DomainError
        If the CEA database has no record (code ``transport_data_missing``).

    References
    ----------
    .. [1] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
           https://github.com/nasa/cea
    """
    record = _load()[0].get(str(name).strip().upper())
    if record is None:
        raise DomainError(
            f"No CEA transport data for {name}.", code="transport_data_missing", field="species",
        )
    return record


def binary_interaction(first: str, second: str) -> BinaryInteractionData | None:
    """Return the CEA interaction-viscosity fit of a pair, or ``None`` [1]_.

    Parameters
    ----------
    first, second : str
        Species names in either order (case-insensitive).

    Returns
    -------
    BinaryInteractionData or None
        Bundled fit, or ``None`` when trans.inp has no record for the pair.

    References
    ----------
    .. [1] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
           https://github.com/nasa/cea
    """
    return _load()[1].get(frozenset((first.strip().upper(), second.strip().upper())))


def rigid_sphere_interaction_viscosity(
    viscosity_i_pa_s: float, viscosity_j_pa_s: float,
    molar_mass_i_kg_mol: float, molar_mass_j_kg_mol: float,
) -> float:
    """Interaction viscosity implied by RP-1311 equation (5.5) [1]_ [2]_.

    Parameters
    ----------
    viscosity_i_pa_s, viscosity_j_pa_s : float
        Pure-species viscosities in Pa s.
    molar_mass_i_kg_mol, molar_mass_j_kg_mol : float
        Molar masses in kg/mol.

    Returns
    -------
    float
        ``eta_ij = eta_i 2 M_j / ((M_i + M_j) phi_ij)`` with ``phi_ij`` of
        equation (5.5); the value is symmetric in ``i`` and ``j``.

    Notes
    -----
    Substituting equation (5.5) into equation (5.7) gives
    ``eta_ij = 4 sqrt(2) eta_i sqrt(M_j/(M_i+M_j)) / (1 + sqrt(eta_i/eta_j) (M_j/M_i)**0.25)**2``,
    which is the estimate CEA applies to pairs without a trans.inp record
    (``source/equilibrium.f90``, ``compute_transport_properties``, constant
    ``5.656854d0``) [2]_. Using it in equation (5.7) therefore reproduces
    equation (5.5) exactly.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (5.5) and (5.7). https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
           ``source/equilibrium.f90``, ``compute_transport_properties``.
           https://github.com/nasa/cea
    """
    phi = 0.25 * (
        1.0 + sqrt(viscosity_i_pa_s / viscosity_j_pa_s)
        * (molar_mass_j_kg_mol / molar_mass_i_kg_mol) ** 0.25
    ) ** 2 * sqrt(2.0 * molar_mass_j_kg_mol / (molar_mass_i_kg_mol + molar_mass_j_kg_mol))
    return viscosity_i_pa_s * 2.0 * molar_mass_j_kg_mol / (
        (molar_mass_i_kg_mol + molar_mass_j_kg_mol) * phi
    )


@dataclass(frozen=True, slots=True)
class MixtureTransportProperties:
    """Mixture transport properties of RP-1311 Part I chapter 5 [1]_.

    Attributes
    ----------
    temperature_k : float
        Temperature in K.
    mole_fractions : tuple of (str, float)
        Species used and their mole fractions renormalized over that set
        (``x_i`` of equations (5.3)--(5.18)), in decreasing order.
    viscosity_pa_s : float
        Mixture viscosity, equation (5.3).
    frozen_thermal_conductivity_w_m_k : float
        ``lambda_fr``, equation (5.4).
    reaction_thermal_conductivity_w_m_k : float
        ``lambda_re``, equations (5.8)--(5.13); zero for a frozen evaluation.
    equilibrium_thermal_conductivity_w_m_k : float
        ``lambda_eq = lambda_fr + lambda_re``, equation (5.2).
    frozen_cp_j_kg_k, reaction_cp_j_kg_k, equilibrium_cp_j_kg_k : float
        Gas-only specific heats of equations (5.15), (5.16)--(5.18) and (5.14).
    frozen_prandtl, equilibrium_prandtl : float
        Equations (5.19a) and (5.19b).
    reaction_count : int
        ``NR``, the number of independent reactions among the species used.
    omitted_species : tuple of (str, float)
        Species left out for lack of CEA transport data, with their input
        mole fractions.
    warnings : tuple of WarningMessage
        Omissions and extrapolations.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           chapter 5. https://ntrs.nasa.gov/citations/19950013764
    """

    temperature_k: float
    mole_fractions: tuple[tuple[str, float], ...]
    viscosity_pa_s: float
    frozen_thermal_conductivity_w_m_k: float
    reaction_thermal_conductivity_w_m_k: float
    equilibrium_thermal_conductivity_w_m_k: float
    frozen_cp_j_kg_k: float
    reaction_cp_j_kg_k: float
    equilibrium_cp_j_kg_k: float
    frozen_prandtl: float
    equilibrium_prandtl: float
    reaction_count: int
    omitted_species: tuple[tuple[str, float], ...]
    warnings: tuple[WarningMessage, ...]


def _solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting; raises on a singular matrix."""
    size = len(vector)
    rows = [row[:] + [value] for row, value in zip(matrix, vector, strict=True)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    for column in range(size):
        pivot = max(range(column, size), key=lambda index: abs(rows[index][column]))
        if abs(rows[pivot][column]) <= 1.0e-13 * max(scale, 1.0e-300):
            raise DomainError(
                "Singular reaction matrix in the transport calculation.",
                code="transport_singular_matrix",
            )
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for index in range(column + 1, size):
            factor = rows[index][column] / rows[column][column]
            if factor:
                for k in range(column, size + 1):
                    rows[index][k] -= factor * rows[column][k]
    solution = [0.0] * size
    for index in range(size - 1, -1, -1):
        total = rows[index][size] - sum(
            rows[index][k] * solution[k] for k in range(index + 1, size)
        )
        solution[index] = total / rows[index][index]
    return solution


def independent_reactions(
    formulas: list[Mapping[str, float]], *, tolerance: float = 1.0e-9
) -> list[list[float]]:
    """Return a basis of reactions ``sum_k alpha_ik A_k = 0`` among species [1]_.

    Parameters
    ----------
    formulas : list of mapping of str to float
        Element counts of each species, in the order of preference for the
        component (basis) species; put the most abundant species first.
    tolerance : float, default 1e-9
        Pivot tolerance of the row reduction.

    Returns
    -------
    list of list of float
        ``NR`` rows of stoichiometric coefficients ``alpha_ik`` (equation
        (5.10)); each row forms one non-component species (coefficient -1)
        from the component species, so the rows are independent and span all
        element-conserving reactions.

    Notes
    -----
    The rows are the null-space basis of the element matrix obtained by
    reduced row-echelon elimination. Equations (5.8)--(5.18) are invariant to
    the choice of basis: replacing the reactions by ``B alpha`` with ``B``
    non-singular changes ``Delta_r H`` to ``B Delta_r H`` and ``g`` to
    ``B g B^T``, leaving ``Delta_r H^T g^-1 Delta_r H`` unchanged.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (5.9) and (5.10). https://ntrs.nasa.gov/citations/19950013764
    """
    elements = sorted({element for formula in formulas for element in formula})
    columns = len(formulas)
    matrix = [[float(formula.get(element, 0.0)) for formula in formulas] for element in elements]
    pivots: list[int] = []
    row = 0
    for column in range(columns):
        if row >= len(matrix):
            break
        best = max(range(row, len(matrix)), key=lambda index: abs(matrix[index][column]))
        if abs(matrix[best][column]) <= tolerance:
            continue
        matrix[row], matrix[best] = matrix[best], matrix[row]
        pivot = matrix[row][column]
        matrix[row] = [value / pivot for value in matrix[row]]
        for index in range(len(matrix)):
            if index != row and matrix[index][column]:
                factor = matrix[index][column]
                matrix[index] = [a - factor * b for a, b in zip(matrix[index], matrix[row], strict=True)]
        pivots.append(column)
        row += 1
    reactions = []
    for free in (column for column in range(columns) if column not in pivots):
        coefficients = [0.0] * columns
        coefficients[free] = -1.0
        for position, pivot_column in enumerate(pivots):
            coefficients[pivot_column] = matrix[position][free]
        reactions.append(coefficients)
    return reactions


def mixture_transport_properties(
    mole_fractions: Mapping[str, float],
    temperature_k: float,
    *,
    reacting: bool = True,
    minimum_mole_fraction: float = DEFAULT_MINIMUM_MOLE_FRACTION,
    maximum_omitted_mole_fraction: float = DEFAULT_MAXIMUM_OMITTED_MOLE_FRACTION,
    allow_extrapolation: bool = False,
) -> MixtureTransportProperties:
    """Viscosity, conductivity, specific heat and Prandtl number of a gas mixture [1]_.

    Parameters
    ----------
    mole_fractions : mapping of str to float
        Non-negative mole fraction (or amount) per species-catalog key; they are
        normalized internally.
    temperature_k : float
        Temperature in K.
    reacting : bool, default True
        Add the reaction contributions of local chemical equilibrium,
        equations (5.8)--(5.18). ``False`` gives frozen values only.
    minimum_mole_fraction : float, default 1e-10
        Species below this normalized mole fraction are left out (they
        contribute in proportion to their mole fraction).
    maximum_omitted_mole_fraction : float, default 1e-3
        Largest total mole fraction of species without CEA transport data that
        may be omitted (with a warning); more raises :class:`DomainError`.
        Species absent from the bundled thermodynamic catalog of
        :mod:`.species` count as omitted too.
    allow_extrapolation : bool, default False
        Evaluate species and interaction fits outside their temperature range
        as CEA does, with a warning. When false, a species fit outside its
        range raises, and an interaction fit outside its range is replaced by
        the equation (5.5) estimate (RP-1311 uses (5.5) whenever ``eta_ij`` is
        not available) with a warning.

    Returns
    -------
    MixtureTransportProperties
        Mixture properties and the species actually used.

    Raises
    ------
    InputError
        If the temperature or a mole fraction is invalid or nothing is left.
    DomainError
        If too much of the mixture lacks transport data
        (``transport_data_missing``), a fit or thermodynamic polynomial does not
        cover the temperature, or the reaction matrix is singular.

    Notes
    -----
    Pure-species viscosity ``eta_i`` and conductivity ``lambda_i`` come from
    equation (5.1). The mixture viscosity and frozen conductivity are
    equations (5.3) and (5.4) with ``phi_ij`` from the interaction viscosity
    ``eta_ij`` of trans.inp through equation (5.7) when available, otherwise
    from equation (5.5), and ``psi_ij`` from equation (5.6). The reaction
    conductivity is equation (5.8) with ``lambda_r`` from equations (5.11)
    to (5.13) and ``A*_kl = 1.1``; CEA's source multiplies its reaction sum by
    ``0.6 = 3/5`` and divides by ``1.1`` exactly as equation (5.13) requires
    [2]_. The reaction specific heat is equations (5.16)--(5.18), the frozen
    one equation (5.15) with ``C_p`` of each species from the NASA Glenn
    polynomials of :mod:`.species` [3]_, and the Prandtl numbers equations
    (5.19a,b). Reactions are an independent basis built by
    :func:`independent_reactions`, components chosen in decreasing mole
    fraction (the result does not depend on that choice).

    Differences from CEA 3.3.4, which are documented and measured in
    ``docs/TRANSPORT_AND_DERIVATIVES.md``: CEA estimates properties of
    species missing from trans.inp from an undocumented correlation, whereas
    here they are omitted (fail-closed above ``maximum_omitted_mole_fraction``);
    and CEA selects its transport species by cumulative abundance, whereas
    here every species above ``minimum_mole_fraction`` is used.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           chapter 5, equations (5.1)--(5.19). https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
           ``source/equilibrium.f90``, ``compute_transport_properties``.
           https://github.com/nasa/cea
    .. [3] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise InputError("Temperature must be finite and positive.", field="temperature_k")
    if not isfinite(minimum_mole_fraction) or not 0.0 <= minimum_mole_fraction < 1.0:
        raise InputError("minimum_mole_fraction must be in [0, 1).", field="minimum_mole_fraction")
    if not isfinite(maximum_omitted_mole_fraction) or not 0.0 <= maximum_omitted_mole_fraction <= 1.0:
        raise InputError(
            "maximum_omitted_mole_fraction must be in [0, 1].",
            field="maximum_omitted_mole_fraction",
        )
    amounts: dict[str, float] = {}
    for key, value in mole_fractions.items():
        value = float(value)
        if not isfinite(value) or value < 0.0:
            raise InputError(f"Mole fraction of {key} must be finite and non-negative.",
                             field="mole_fractions")
        normalized_key = str(key).strip().upper()
        amounts[normalized_key] = amounts.get(normalized_key, 0.0) + value
    total = sum(amounts.values())
    if total <= 0.0:
        raise InputError("At least one mole fraction must be positive.", field="mole_fractions")
    catalog, _ = _load()
    thermo_catalog = species_catalog()
    warnings: list[WarningMessage] = []
    omitted = tuple(sorted(
        ((key, value / total) for key, value in amounts.items()
         if value / total >= minimum_mole_fraction
         and (key not in catalog or key not in thermo_catalog)),
        key=lambda item: -item[1],
    ))
    omitted_total = sum(value for _, value in omitted)
    if omitted_total > maximum_omitted_mole_fraction:
        raise DomainError(
            "Too much of the mixture has no CEA transport data.",
            code="transport_data_missing",
            field="mole_fractions",
            details={"omitted_species": dict(omitted), "omitted_mole_fraction": omitted_total},
        )
    if omitted:
        warnings.append(WarningMessage(
            code="transport_species_omitted",
            message=(
                f"Species without CEA transport or bundled thermodynamic data were omitted "
                f"(total mole fraction {omitted_total:.3g}): "
                + ", ".join(key for key, _ in omitted)
            ),
        ))
    used = sorted(
        ((key, value) for key, value in amounts.items()
         if value / total >= minimum_mole_fraction
         and key in catalog and key in thermo_catalog),
        key=lambda item: -item[1],
    )
    if not used:
        raise InputError("No species with transport data remain.", field="mole_fractions")
    used_total = sum(value for _, value in used)
    names = [key for key, _ in used]
    x = [value / used_total for _, value in used]
    count = len(names)
    thermo = [get_species(name) for name in names]
    molar_mass = [species.molecular_mass_kg_mol for species in thermo]
    properties = [species.properties(temperature_k) for species in thermo]

    extrapolated: list[str] = []

    def evaluate(intervals, owner, quantity, scale):
        interval, outside = _select_interval(
            intervals, temperature_k, owner=owner, quantity=quantity,
            allow_extrapolation=allow_extrapolation,
        )
        if outside:
            extrapolated.append(owner)
        return exp(interval.log_value(temperature_k)) * scale

    eta = [[0.0] * count for _ in range(count)]
    conductivity = [0.0] * count
    for index, name in enumerate(names):
        record = catalog[name]
        if not record.viscosity or not record.thermal_conductivity:
            raise DomainError(f"{record.name} lacks a viscosity or conductivity fit.",
                              code="transport_data_missing", field="mole_fractions")
        eta[index][index] = evaluate(record.viscosity, record.name, "Viscosity", MICROPOISE_TO_PA_S)
        conductivity[index] = evaluate(
            record.thermal_conductivity, record.name, "Thermal-conductivity",
            MICROWATT_PER_CM_K_TO_W_PER_M_K,
        )
    estimated_pairs: list[str] = []
    for i in range(count):
        for j in range(i + 1, count):
            pair = binary_interaction(names[i], names[j])
            usable = pair is not None and (
                allow_extrapolation
                or pair.viscosity[0].minimum_temperature_k <= temperature_k
                <= pair.viscosity[-1].maximum_temperature_k
            )
            if pair is not None and not usable:
                estimated_pairs.append("-".join(pair.pair))
            if usable:
                value = evaluate(pair.viscosity, "-".join(pair.pair), "Interaction-viscosity",
                                 MICROPOISE_TO_PA_S)
            else:
                value = rigid_sphere_interaction_viscosity(
                    eta[i][i], eta[j][j], molar_mass[i], molar_mass[j]
                )
            eta[i][j] = eta[j][i] = value
    if estimated_pairs:
        warnings.append(WarningMessage(
            code="transport_interaction_estimated",
            message="Interaction-viscosity fits do not cover the temperature; equation (5.5) "
            "was used instead for: " + ", ".join(estimated_pairs),
            field="temperature_k",
        ))
    if extrapolated:
        warnings.append(WarningMessage(
            code="transport_fit_extrapolated",
            message="Transport fits were extrapolated beyond their temperature range: "
            + ", ".join(sorted(set(extrapolated))),
            field="temperature_k",
        ))

    viscosity = 0.0
    frozen_conductivity = 0.0
    for i in range(count):
        sum_phi = x[i]
        sum_psi = x[i]
        for j in range(count):
            if j == i:
                continue
            phi = eta[i][i] * 2.0 * molar_mass[j] / (eta[i][j] * (molar_mass[i] + molar_mass[j]))
            psi = phi * (1.0 + 2.41 * (molar_mass[i] - molar_mass[j])
                         * (molar_mass[i] - 0.142 * molar_mass[j])
                         / (molar_mass[i] + molar_mass[j]) ** 2)
            sum_phi += x[j] * phi
            sum_psi += x[j] * psi
        viscosity += x[i] * eta[i][i] / sum_phi
        frozen_conductivity += x[i] * conductivity[i] / sum_psi

    mean_molar_mass = sum(xi * mi for xi, mi in zip(x, molar_mass, strict=True))
    frozen_cp = sum(
        xi * prop.cp_molar_j_mol_k for xi, prop in zip(x, properties, strict=True)
    ) / mean_molar_mass

    reaction_conductivity = 0.0
    reaction_cp = 0.0
    reactions = (
        independent_reactions([dict(species.formula) for species in thermo]) if reacting else []
    )
    if reactions:
        rt = MOLAR_GAS_CONSTANT_J_MOL_K * temperature_k
        enthalpy_over_rt = [
            sum(alpha[k] * properties[k].enthalpy_molar_j_mol for k in range(count)) / rt
            for alpha in reactions
        ]
        size = len(reactions)
        g = [[0.0] * size for _ in range(size)]
        d = [[0.0] * size for _ in range(size)]
        for k in range(count - 1):
            for m in range(k + 1, count):
                rt_over_pd = 5.0 * molar_mass[k] * molar_mass[m] / (
                    3.0 * COLLISION_CROSS_SECTION_RATIO * eta[k][m]
                    * (molar_mass[k] + molar_mass[m])
                )
                terms = [alpha[k] / x[k] - alpha[m] / x[m] for alpha in reactions]
                weight = x[k] * x[m]
                for i in range(size):
                    if not terms[i]:
                        continue
                    for j in range(size):
                        product = weight * terms[i] * terms[j]
                        d[i][j] += product
                        g[i][j] += rt_over_pd * product
        lambda_r = _solve(g, enthalpy_over_rt)
        capital_x = _solve(d, enthalpy_over_rt)
        reaction_conductivity = MOLAR_GAS_CONSTANT_J_MOL_K * sum(
            h * value for h, value in zip(enthalpy_over_rt, lambda_r, strict=True)
        )
        reaction_cp = MOLAR_GAS_CONSTANT_J_MOL_K * sum(
            h * value for h, value in zip(enthalpy_over_rt, capital_x, strict=True)
        ) / mean_molar_mass

    equilibrium_conductivity = frozen_conductivity + reaction_conductivity
    equilibrium_cp = frozen_cp + reaction_cp
    return MixtureTransportProperties(
        temperature_k=temperature_k,
        mole_fractions=tuple(zip(names, x, strict=True)),
        viscosity_pa_s=viscosity,
        frozen_thermal_conductivity_w_m_k=frozen_conductivity,
        reaction_thermal_conductivity_w_m_k=reaction_conductivity,
        equilibrium_thermal_conductivity_w_m_k=equilibrium_conductivity,
        frozen_cp_j_kg_k=frozen_cp,
        reaction_cp_j_kg_k=reaction_cp,
        equilibrium_cp_j_kg_k=equilibrium_cp,
        frozen_prandtl=prandtl_number(frozen_cp, viscosity, frozen_conductivity),
        equilibrium_prandtl=prandtl_number(equilibrium_cp, viscosity, equilibrium_conductivity),
        reaction_count=len(reactions),
        omitted_species=omitted,
        warnings=tuple(warnings),
    )


def equilibrium_transport_properties(
    result: EquilibriumResult, **options: object
) -> MixtureTransportProperties:
    """Transport properties of an equilibrium result's gas composition [1]_.

    Parameters
    ----------
    result : EquilibriumResult
        Result of an equilibrium provider whose species keys are in the bundled
        species catalog (for example :class:`.LocalEquilibriumProvider`).
    **options
        Keyword options of :func:`mixture_transport_properties`.

    Returns
    -------
    MixtureTransportProperties
        Properties at ``result.temperature_k`` and the result's mole fractions.

    Raises
    ------
    DomainError, InputError
        As :func:`mixture_transport_properties`.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           chapter 5. https://ntrs.nasa.gov/citations/19950013764
    """
    return mixture_transport_properties(
        {amount.species: amount.mole_fraction for amount in result.species},
        result.temperature_k,
        **options,
    )


def prandtl_number(cp_j_kg_k: float, viscosity_pa_s: float, conductivity_w_m_k: float) -> float:
    """Prandtl number ``cp mu / lambda``, RP-1311 equations (5.19a,b) [1]_.

    Parameters
    ----------
    cp_j_kg_k : float
        Specific heat in J/(kg K).
    viscosity_pa_s : float
        Dynamic viscosity in Pa s.
    conductivity_w_m_k : float
        Thermal conductivity in W/(m K).

    Returns
    -------
    float
        Dimensionless Prandtl number.

    Raises
    ------
    InputError
        If any argument is not finite and positive.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (5.19a) and (5.19b). https://ntrs.nasa.gov/citations/19950013764
    """
    for name, value in (("cp_j_kg_k", cp_j_kg_k), ("viscosity_pa_s", viscosity_pa_s),
                        ("conductivity_w_m_k", conductivity_w_m_k)):
        if not isfinite(value) or value <= 0.0:
            raise InputError(f"{name} must be finite and positive.", field=name)
    return cp_j_kg_k * viscosity_pa_s / conductivity_w_m_k


#: U.S. Standard Atmosphere 1976 eq. (51) constant beta, kg/(s m K^0.5).
USSA1976_SUTHERLAND_BETA = 1.458e-6
#: U.S. Standard Atmosphere 1976 Sutherland constant S, K (eq. (51), table 2B).
USSA1976_SUTHERLAND_CONSTANT_K = 110.4


def sutherland_viscosity(
    temperature_k: float,
    *,
    beta: float = USSA1976_SUTHERLAND_BETA,
    sutherland_constant_k: float = USSA1976_SUTHERLAND_CONSTANT_K,
) -> float:
    """Sutherland's law ``mu = beta T**1.5 / (T + S)`` [1]_.

    Parameters
    ----------
    temperature_k : float
        Temperature in K.
    beta : float, default 1.458e-6
        Constant in kg/(s m K^0.5); the default is the U.S. Standard Atmosphere
        1976 value for air.
    sutherland_constant_k : float, default 110.4
        Sutherland's constant ``S`` in K; default for air from the same source.

    Returns
    -------
    float
        Dynamic viscosity in Pa s.

    Raises
    ------
    InputError
        If an argument is not finite and positive.

    Notes
    -----
    Equation (51) of the U.S. Standard Atmosphere 1976, whose text states that
    the equation fails at very high and very low temperatures; the standard
    applies it to air from 86 km down to sea level. The printed table 2B
    misprints ``S``; the NASA errata appended to the NTRS copy confirms
    110.4 K, the value also stated in the text of section 1.3.11 [1]_.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", section 1.3.11,
           equation (51). https://ntrs.nasa.gov/citations/19770009539
    """
    for name, value in (("temperature_k", temperature_k), ("beta", beta),
                        ("sutherland_constant_k", sutherland_constant_k)):
        if not isfinite(value) or value <= 0.0:
            raise InputError(f"{name} must be finite and positive.", field=name)
    return beta * temperature_k ** 1.5 / (temperature_k + sutherland_constant_k)


#: Neufeld et al. (1972) Omega(2,2)* fit constants A-F, as reproduced by IDAES.
NEUFELD_VISCOSITY_COEFFICIENTS = (1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787)
#: Chapman-Enskog prefactor 26.69 uP A^2 / sqrt(g/mol K) (Poling et al. eq. 9-3.9 via IDAES).
CHAPMAN_ENSKOG_VISCOSITY_PREFACTOR = 26.69


def neufeld_collision_integral(reduced_temperature: float) -> float:
    """Lennard-Jones viscosity collision integral ``Omega(2,2)*`` [1]_.

    Parameters
    ----------
    reduced_temperature : float
        ``T* = k T / epsilon``.

    Returns
    -------
    float
        ``A T*^-B + C exp(-D T*) + E exp(-F T*)`` with the Neufeld, Janzen and
        Aziz (1972) constants.

    Raises
    ------
    InputError
        If ``T*`` is not finite and positive.

    Notes
    -----
    The constants were read from the IDAES documentation and source, which
    cite Poling, Prausnitz and O'Connell, *The Properties of Gases and
    Liquids*, eq. 9-4.3; the original J. Chem. Phys. 57, 1100 (1972) paper
    was not opened, and no validity range of ``T*`` was verified [1]_.

    References
    ----------
    .. [1] IDAES-CHAPMAN-ENSKOG. IDAES-PSE documentation, "Chapman Enskog
           Viscosity Model". https://idaes-pse.readthedocs.io/en/main/explanations/components/property_package/general/transport_properties/chapman_enskog.html
    """
    if not isfinite(reduced_temperature) or reduced_temperature <= 0.0:
        raise InputError("Reduced temperature must be finite and positive.",
                         field="reduced_temperature")
    a, b, c, d, e, f = NEUFELD_VISCOSITY_COEFFICIENTS
    t = reduced_temperature
    return a * t ** (-b) + c * exp(-d * t) + e * exp(-f * t)


def chapman_enskog_viscosity(
    temperature_k: float,
    *,
    molar_mass_kg_mol: float,
    collision_diameter_m: float,
    well_depth_k: float,
) -> float:
    """Dilute-gas Chapman--Enskog viscosity with a Lennard-Jones potential [1]_.

    Parameters
    ----------
    temperature_k : float
        Temperature in K.
    molar_mass_kg_mol : float
        Molar mass in kg/mol.
    collision_diameter_m : float
        Lennard-Jones ``sigma`` in m.
    well_depth_k : float
        Lennard-Jones ``epsilon/k`` in K.

    Returns
    -------
    float
        ``26.69 sqrt(M T) / (sigma**2 Omega_v)`` micropoise (``M`` in g/mol,
        ``sigma`` in angstrom), converted to Pa s.

    Raises
    ------
    InputError
        If an argument is not finite and positive.

    Notes
    -----
    ``sigma`` and ``epsilon/k`` must come from the same source; they are
    under-determined by viscosity data alone, as the cited page notes [1]_.

    References
    ----------
    .. [1] IDAES-CHAPMAN-ENSKOG. IDAES-PSE documentation, "Chapman Enskog
           Viscosity Model". https://idaes-pse.readthedocs.io/en/main/explanations/components/property_package/general/transport_properties/chapman_enskog.html
    """
    for name, value in (("temperature_k", temperature_k), ("molar_mass_kg_mol", molar_mass_kg_mol),
                        ("collision_diameter_m", collision_diameter_m),
                        ("well_depth_k", well_depth_k)):
        if not isfinite(value) or value <= 0.0:
            raise InputError(f"{name} must be finite and positive.", field=name)
    omega = neufeld_collision_integral(temperature_k / well_depth_k)
    molar_mass_g_mol = molar_mass_kg_mol * 1.0e3
    sigma_angstrom = collision_diameter_m * 1.0e10
    micropoise = CHAPMAN_ENSKOG_VISCOSITY_PREFACTOR * sqrt(molar_mass_g_mol * temperature_k) / (
        sigma_angstrom ** 2 * omega
    )
    return micropoise * MICROPOISE_TO_PA_S


def modified_eucken_conductivity(
    viscosity_pa_s: float, cp_molar_j_mol_k: float, molar_mass_kg_mol: float
) -> float:
    """Conductivity estimate CEA applies to gases without trans.inp data [1]_.

    Parameters
    ----------
    viscosity_pa_s : float
        Species viscosity in Pa s.
    cp_molar_j_mol_k : float
        Ideal-gas molar heat capacity in J/(mol K).
    molar_mass_kg_mol : float
        Molar mass in kg/mol.

    Returns
    -------
    float
        ``mu (R/M) (15/4 + 1.32 (cp/R - 5/2))`` in W/(m K).

    Raises
    ------
    InputError
        If an argument is not finite and positive.

    Notes
    -----
    CEA's ``compute_transport_properties`` sets the missing conductivity of a
    species to ``eta R (0.00375 + 0.00132 (cp/R - 2.5)) / M`` with ``eta`` in
    micropoise, ``R`` in J/(kmol K) and the result in uW/(cm K); in SI units
    this is the expression above (a modified Eucken form). RP-1311 Part I does
    not document it; it is offered only as a labelled estimate [1]_.

    References
    ----------
    .. [1] NASA-CEA-SRC-TRANSPORT. NASA Glenn Research Center, CEA source
           ``source/equilibrium.f90``, ``compute_transport_properties``.
           https://github.com/nasa/cea
    """
    for name, value in (("viscosity_pa_s", viscosity_pa_s), ("cp_molar_j_mol_k", cp_molar_j_mol_k),
                        ("molar_mass_kg_mol", molar_mass_kg_mol)):
        if not isfinite(value) or value <= 0.0:
            raise InputError(f"{name} must be finite and positive.", field=name)
    cp_over_r = cp_molar_j_mol_k / MOLAR_GAS_CONSTANT_J_MOL_K
    return viscosity_pa_s * MOLAR_GAS_CONSTANT_J_MOL_K / molar_mass_kg_mol * (
        3.75 + 1.32 * (cp_over_r - 2.5)
    )
