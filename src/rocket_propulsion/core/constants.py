"""Physical constants, each with the primary source that defines its value.

Every physical constant used by the library lives here, so each value can be
traced to one source. ``tests/test_constants_and_references.py`` fails if a
registered value is written as a bare literal anywhere else in the package;
other modules import the module-level names (``STANDARD_GRAVITY_M_S2`` and
the like) instead.

Where a model is *defined* with a particular constant, the model's own value
is used even if a newer value exists. The U.S. Standard Atmosphere 1976 [2]_
is defined with its own gas constant, ``USSA1976_GAS_CONSTANT``, rather than
the 2019 SI value ``MOLAR_GAS_CONSTANT_SI_2019`` [1]_, and its tabulated layer
pressures are reproduced only with the former.

References
----------
.. [1] BIPM, "The International System of Units (SI Brochure)", 9th ed., 2019,
       https://www.bipm.org/en/publications/si-brochure
.. [2] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
       NASA-TM-X-74335, https://ntrs.nasa.gov/citations/19770009539
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PhysicalConstant:
    """One registered physical constant with its unit and defining source.

    Attributes
    ----------
    name : str
        Descriptive name, unique within ``REGISTERED_CONSTANTS``.
    value : float
        Numerical value in ``unit``.
    unit : str
        SI unit of ``value``.
    source : str
        Document that defines the value.
    url : str
        HTTPS location of ``source``; it is also listed in
        ``docs/REFERENCES.md``.
    exact : bool
        True when the value is exact by definition in its source, False when
        it is a rounded derived value.
    note : str
        Optional usage note.

    References
    ----------
    .. [1] BIPM, "The International System of Units (SI Brochure)", 9th ed.,
           2019, https://www.bipm.org/en/publications/si-brochure
    .. [2] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335, https://ntrs.nasa.gov/citations/19770009539
    """

    name: str
    value: float
    unit: str
    source: str
    url: str
    exact: bool
    note: str = ""


STANDARD_GRAVITY = PhysicalConstant(
    "standard acceleration of gravity g_n",
    9.80665,
    "m/s^2",
    "3rd CGPM (1901), restated in the BIPM SI Brochure, 9th ed. (2019), and used "
    "as g0 by the U.S. Standard Atmosphere 1976",
    "https://www.bipm.org/en/publications/si-brochure",
    exact=True,
    note="Conventional value; Isp = F / (mdot g_n).",
)

MOLAR_GAS_CONSTANT_SI_2019 = PhysicalConstant(
    "molar gas constant R = N_A k",
    8.31446261815324,
    "J/(mol K)",
    "BIPM SI Brochure, 9th ed. (2019): N_A and k are exact, so R = N_A k is exact",
    "https://www.bipm.org/en/publications/si-brochure",
    exact=True,
    note="Product 6.02214076e23 * 1.380649e-23, rounded to double precision.",
)

USSA1976_GAS_CONSTANT = PhysicalConstant(
    "U.S. Standard Atmosphere 1976 gas constant R*",
    8.31432,
    "J/(mol K)",
    "U.S. Standard Atmosphere, 1976, Part 1 (R* = 8.31432e3 N m / (kmol K))",
    "https://ntrs.nasa.gov/citations/19770009539",
    exact=True,
    note="Defining constant of the 1976 model; differs from the 2019 SI value by 1.7e-5.",
)

USSA1976_MOLAR_MASS = PhysicalConstant(
    "U.S. Standard Atmosphere 1976 sea-level molar mass of air M0",
    0.0289644,
    "kg/mol",
    "U.S. Standard Atmosphere, 1976, Part 1 (M0 = 28.9644 kg/kmol)",
    "https://ntrs.nasa.gov/citations/19770009539",
    exact=True,
)

AIR_SPECIFIC_GAS_CONSTANT = PhysicalConstant(
    "specific gas constant of dry air R*/M0 (rounded)",
    287.05,
    "J/(kg K)",
    "U.S. Standard Atmosphere, 1976: R*/M0 = 287.0531 J/(kg K), rounded to 287.05",
    "https://ntrs.nasa.gov/citations/19770009539",
    exact=False,
    note="Default for ideal-gas air examples only; pass the gas's own value otherwise.",
)

USSA1976_SEA_LEVEL_PRESSURE = PhysicalConstant(
    "U.S. Standard Atmosphere 1976 sea-level pressure P0",
    101325.0,
    "Pa",
    "U.S. Standard Atmosphere, 1976, Part 1 (P0 = 1.013250e5 N/m^2)",
    "https://ntrs.nasa.gov/citations/19770009539",
    exact=True,
)

USSA1976_SEA_LEVEL_TEMPERATURE = PhysicalConstant(
    "U.S. Standard Atmosphere 1976 sea-level temperature T0",
    288.15,
    "K",
    "U.S. Standard Atmosphere, 1976, Part 1 (T0 = 288.15 K)",
    "https://ntrs.nasa.gov/citations/19770009539",
    exact=True,
)

CELSIUS_OFFSET = PhysicalConstant(
    "Celsius temperature offset",
    273.15,
    "K",
    "BIPM SI Brochure, 9th ed. (2019): t/degC = T/K - 273.15",
    "https://www.bipm.org/en/publications/si-brochure",
    exact=True,
)

REGISTERED_CONSTANTS: tuple[PhysicalConstant, ...] = (
    STANDARD_GRAVITY,
    MOLAR_GAS_CONSTANT_SI_2019,
    USSA1976_GAS_CONSTANT,
    USSA1976_MOLAR_MASS,
    AIR_SPECIFIC_GAS_CONSTANT,
    USSA1976_SEA_LEVEL_PRESSURE,
    USSA1976_SEA_LEVEL_TEMPERATURE,
    CELSIUS_OFFSET,
)

STANDARD_GRAVITY_M_S2 = STANDARD_GRAVITY.value
MOLAR_GAS_CONSTANT_J_MOL_K = MOLAR_GAS_CONSTANT_SI_2019.value
USSA1976_GAS_CONSTANT_J_MOL_K = USSA1976_GAS_CONSTANT.value
USSA1976_MOLAR_MASS_KG_MOL = USSA1976_MOLAR_MASS.value
AIR_SPECIFIC_GAS_CONSTANT_J_KG_K = AIR_SPECIFIC_GAS_CONSTANT.value
CELSIUS_OFFSET_K = CELSIUS_OFFSET.value
USSA1976_SEA_LEVEL_PRESSURE_PA = USSA1976_SEA_LEVEL_PRESSURE.value
USSA1976_SEA_LEVEL_TEMPERATURE_K = USSA1976_SEA_LEVEL_TEMPERATURE.value
