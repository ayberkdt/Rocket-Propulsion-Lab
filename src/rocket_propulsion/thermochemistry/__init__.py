"""Temperature-dependent species and ideal-mixture thermochemistry.

The package evaluates NASA Glenn polynomial species properties [3]_, ideal-gas
mixtures, an experimental local Gibbs-energy equilibrium solver [1]_, an
adapter for a user-installed NASA CEA executable [2]_, and frozen versus
equilibrium nozzle-expansion comparisons. Each public callable documents its
own assumptions and sources.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part II".
       https://ntrs.nasa.gov/citations/19960044559
.. [3] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

# Higher-level submodules import the names above, so they load last.
from . import (
    combustion_efficiency,
    condensed_equilibrium,
    finite_area_combustor,
    kinetics,
    nozzle_kinetics,
)
from .cea import (
    CEA_MANUAL_URL,
    CeaEquilibriumProvider,
    CeaPlotTable,
    CeaReactant,
    CeaRocketProblem,
    CeaRunArtifact,
    CeaSubprocessAdapter,
    cea_output_warnings,
    parse_cea_plot,
    parse_fortran_number,
    render_cea_equilibrium_input,
    render_cea_rocket_input,
)
from .cea_library import CeaLibraryInfo, cea_library_info, solve_cea_rocket
from .equilibrium import (
    EquilibriumConvergence,
    EquilibriumProblem,
    EquilibriumProvider,
    EquilibriumResiduals,
    EquilibriumResult,
    EquilibriumSpeciesAmount,
    LocalEquilibriumProvider,
    ReactantAmount,
)
from .equilibrium_derivatives import EquilibriumDerivatives, equilibrium_derivatives
from .mixtures import Mixture, MixtureComponent, MixtureProperties
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K, NasaPolynomialInterval
from .rocket_equilibrium import (
    PROPELLANT_PAIRS,
    ChemistryNozzleComparison,
    EquilibriumExpansionResult,
    OxidizerFuelSweep,
    OxidizerFuelSweepPoint,
    PropellantPair,
    compare_nozzle_chemistry,
    generate_oxidizer_fuel_sweep,
    get_propellant_pair,
    reactants_for_mass_ratio,
)
from .rocket_performance import (
    RocketExit,
    RocketPerformance,
    RocketStation,
    solve_rocket_performance,
)
from .species import (
    CondensedReactant,
    Species,
    SpeciesProperties,
    condensed_reactant_catalog,
    get_reactant,
    get_species,
    species_catalog,
)
from .transport import (
    MixtureTransportProperties,
    equilibrium_transport_properties,
    mixture_transport_properties,
    species_transport_catalog,
)

__all__ = [
    "CEA_MANUAL_URL",
    "PROPELLANT_PAIRS",
    "UNIVERSAL_GAS_CONSTANT_J_MOL_K",
    "CeaEquilibriumProvider",
    "CeaLibraryInfo",
    "CeaPlotTable",
    "CeaReactant",
    "CeaRocketProblem",
    "CeaRunArtifact",
    "CeaSubprocessAdapter",
    "ChemistryNozzleComparison",
    "CondensedReactant",
    "EquilibriumConvergence",
    "EquilibriumDerivatives",
    "EquilibriumExpansionResult",
    "EquilibriumProblem",
    "EquilibriumProvider",
    "EquilibriumResiduals",
    "EquilibriumResult",
    "EquilibriumSpeciesAmount",
    "LocalEquilibriumProvider",
    "Mixture",
    "MixtureComponent",
    "MixtureProperties",
    "MixtureTransportProperties",
    "NasaPolynomialInterval",
    "OxidizerFuelSweep",
    "OxidizerFuelSweepPoint",
    "PropellantPair",
    "ReactantAmount",
    "RocketExit",
    "RocketPerformance",
    "RocketStation",
    "Species",
    "SpeciesProperties",
    "cea_library_info",
    "cea_output_warnings",
    "combustion_efficiency",
    "compare_nozzle_chemistry",
    "condensed_equilibrium",
    "condensed_reactant_catalog",
    "equilibrium_derivatives",
    "equilibrium_transport_properties",
    "finite_area_combustor",
    "generate_oxidizer_fuel_sweep",
    "get_propellant_pair",
    "get_reactant",
    "get_species",
    "kinetics",
    "mixture_transport_properties",
    "nozzle_kinetics",
    "parse_cea_plot",
    "parse_fortran_number",
    "reactants_for_mass_ratio",
    "render_cea_equilibrium_input",
    "render_cea_rocket_input",
    "solve_cea_rocket",
    "solve_rocket_performance",
    "species_catalog",
    "species_transport_catalog",
]
