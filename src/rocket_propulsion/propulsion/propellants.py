"""Curated screening records for widely used chemical rocket propellant systems.

Notes
-----
Vacuum Isp entries are deliberately broad preliminary-design ranges for each
propellant family [2]_; categories follow the liquid/solid/hybrid
classification of [3]_. Where a NASA briefing prints a single representative
value it is stored separately with its page [1]_. Records that the bundled
thermochemistry can compute name their ``thermochemistry_pair``; ideal
performance for those comes from the CEA-gated equilibrium solver, not from
the screening range. The remaining records state why they cannot be
computed. Mixture-ratio notes describe screening regions, not operating
instructions.

References
----------
.. [1] KSC-2010-233. L. G. Delgado, "Expendable Launch Vehicles Briefing and
       Basic Rocketry Physics", NASA Kennedy Space Center, 2010, NTRS
       20100036816, PDF p.46 "Examples of Specific Impulse".
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [3] NASA, "2020 NASA Technology Taxonomy", TX01 Propulsion Systems (liquid,
       solid, hybrid, cryogenic, storable and gelled categories).
"""

from dataclasses import dataclass

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class PropellantRecord:
    """High-level screening data for one chemical propulsion-system family.

    Attributes
    ----------
    key : str
        Stable identifier, for example ``'lox-lh2'``.
    name : str
        Display name.
    category : {'liquid', 'solid', 'hybrid'}
        Propulsion-system category.
    fuel, oxidizer : str
        Fuel and oxidizer descriptions.
    storage : str
        Storage condition of fuel and oxidizer.
    ignition : str
        Ignition mode.
    vacuum_isp_min_s, vacuum_isp_max_s : float
        Bounds of the representative vacuum specific impulse range, s.
    mixture_ratio_note : str
        Descriptive oxidizer/fuel screening region.
    strengths, limitations : tuple of str
        Qualitative screening notes.
    thermochemistry_pair : str or None
        Key of the matching
        :data:`~rocket_propulsion.thermochemistry.PROPELLANT_PAIRS` preset
        when the bundled NASA data can compute this system, else ``None``.
    computability_note : str
        Why the system can or cannot be computed by the bundled solver.
    published_isp_s : float or None
        Single representative Isp printed by [1]_, when the briefing lists the
        family; the briefing does not state vacuum or sea-level conditions.
    published_isp_source : str
        Locator of ``published_isp_s``.

    Notes
    -----
    The screening ranges are broad preliminary-design bands [2]_ and are not
    substitutes for equilibrium chemistry, engine test data or a validated
    cycle/nozzle model; categories follow [3]_.

    References
    ----------
    .. [1] KSC-2010-233. L. G. Delgado, "Expendable Launch Vehicles Briefing and
           Basic Rocketry Physics", NASA Kennedy Space Center, 2010, NTRS
           20100036816, PDF p.46.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [3] NASA, "2020 NASA Technology Taxonomy", TX01 Propulsion Systems (liquid,
           solid, hybrid, cryogenic, storable and gelled categories).
    """

    key: str
    name: str
    category: str
    fuel: str
    oxidizer: str
    storage: str
    ignition: str
    vacuum_isp_min_s: float
    vacuum_isp_max_s: float
    mixture_ratio_note: str
    strengths: tuple[str, ...]
    limitations: tuple[str, ...]
    thermochemistry_pair: str | None = None
    computability_note: str = ""
    published_isp_s: float | None = None
    published_isp_source: str = ""


PROPELLANTS: tuple[PropellantRecord, ...] = (
    PropellantRecord(
        key="lox-lh2",
        name="LOX / LH₂",
        category="liquid",
        fuel="Liquid hydrogen",
        oxidizer="Liquid oxygen",
        storage="Cryogenic / cryogenic",
        ignition="External ignition",
        vacuum_isp_min_s=430.0,
        vacuum_isp_max_s=465.0,
        mixture_ratio_note="Typical O/F study region: about 5–6.5",
        strengths=("Highest common chemical Isp", "Clean exhaust", "Strong upper-stage heritage"),
        limitations=("Very low fuel density", "Deep cryogenic storage", "Large insulated tanks"),
        thermochemistry_pair="LOX_LH2",
        computability_note="Liquid reactants at NASA CEA storage-state enthalpy; local equilibrium solver gated against NASA CEA 3.3.4.",
        published_isp_s=455.0,
        published_isp_source="KSC-2010-233 (NTRS 20100036816), PDF p.46",
    ),
    PropellantRecord(
        key="lox-rp1",
        name="LOX / RP-1",
        category="liquid",
        fuel="Rocket-grade kerosene (RP-1)",
        oxidizer="Liquid oxygen",
        storage="Ambient fuel / cryogenic oxidizer",
        ignition="External ignition",
        vacuum_isp_min_s=320.0,
        vacuum_isp_max_s=360.0,
        mixture_ratio_note="Typical O/F study region: about 2.3–2.8",
        strengths=("High bulk density", "Strong first-stage performance", "Mature supply chain"),
        limitations=("Coking potential", "Lower Isp than hydrogen", "LOX handling required"),
        thermochemistry_pair="LOX_RP1",
        computability_note="Liquid reactants at NASA CEA storage-state enthalpy; local equilibrium solver gated against NASA CEA 3.3.4.",
        published_isp_s=358.0,
        published_isp_source="KSC-2010-233 (NTRS 20100036816), PDF p.46",
    ),
    PropellantRecord(
        key="lox-lch4",
        name="LOX / LCH₄",
        category="liquid",
        fuel="Liquid methane",
        oxidizer="Liquid oxygen",
        storage="Cryogenic / cryogenic",
        ignition="External ignition",
        vacuum_isp_min_s=340.0,
        vacuum_isp_max_s=380.0,
        mixture_ratio_note="Typical O/F study region: about 3.2–3.8",
        strengths=("Clean-burning hydrocarbon", "Good reusable-engine fit", "Better density than LH₂"),
        limitations=("Cryogenic storage", "Lower Isp than hydrogen", "Boil-off management"),
        thermochemistry_pair="LOX_CH4",
        computability_note="Liquid reactants at NASA CEA storage-state enthalpy; local equilibrium solver gated against NASA CEA 3.3.4.",
    ),
    PropellantRecord(
        key="nto-mmh",
        name="NTO / MMH",
        category="liquid",
        fuel="Monomethylhydrazine",
        oxidizer="Nitrogen tetroxide",
        storage="Storable / storable",
        ignition="Hypergolic",
        vacuum_isp_min_s=300.0,
        vacuum_isp_max_s=340.0,
        mixture_ratio_note="Typical O/F study region: about 1.6–2.1",
        strengths=("Reliable restart", "Long-duration storage", "No separate igniter"),
        limitations=("Highly toxic", "Corrosive oxidizer", "Strict ground-safety burden"),
        thermochemistry_pair="NTO_MMH",
        computability_note="Liquid reactants at NASA CEA storage-state enthalpy; local equilibrium solver gated against NASA CEA 3.3.4.",
        published_isp_s=305.0,
        published_isp_source="KSC-2010-233 (NTRS 20100036816), PDF p.46, listed as hydrazine/N2O4",
    ),
    PropellantRecord(
        key="ap-htpb-al",
        name="AP / HTPB / Al composite",
        category="solid",
        fuel="HTPB binder with aluminium",
        oxidizer="Ammonium perchlorate",
        storage="Cast solid grain",
        ignition="Pyrotechnic or hot-gas igniter",
        vacuum_isp_min_s=240.0,
        vacuum_isp_max_s=290.0,
        mixture_ratio_note="Premixed composite; formulation-specific",
        strengths=("High readiness", "Simple feed system", "High thrust density"),
        limitations=("Limited throttling/restart", "Grain integrity constraints", "Particle-rich exhaust"),
        computability_note=(
            "Not computed: aluminized exhaust has condensed Al2O3, the local solver is "
            "gas-only, and HTPB is not in NASA CEA thermo.inp; formulation data must be supplied."
        ),
        published_isp_s=250.0,
        published_isp_source="KSC-2010-233 (NTRS 20100036816), PDF p.46, listed as solid motor",
    ),
    PropellantRecord(
        key="n2o-htpb",
        name="N₂O / HTPB",
        category="hybrid",
        fuel="Solid HTPB",
        oxidizer="Nitrous oxide",
        storage="Solid fuel / liquefied oxidizer",
        ignition="External ignition",
        vacuum_isp_min_s=220.0,
        vacuum_isp_max_s=280.0,
        mixture_ratio_note="Operating O/F changes during regression",
        strengths=("Separated fuel and oxidizer", "Potential shutdown", "Simpler than bipropellant feed"),
        limitations=("Regression-rate limits", "Mixture-ratio shift", "Oxidizer decomposition hazards"),
        computability_note=(
            "Not computed: HTPB is not in NASA CEA thermo.inp and the hybrid O/F shifts "
            "with grain regression."
        ),
    ),
    PropellantRecord(
        key="lox-paraffin",
        name="LOX / paraffin",
        category="hybrid",
        fuel="Solid paraffin-based fuel",
        oxidizer="Liquid oxygen",
        storage="Solid fuel / cryogenic oxidizer",
        ignition="External ignition",
        vacuum_isp_min_s=260.0,
        vacuum_isp_max_s=330.0,
        mixture_ratio_note="Operating O/F and regression are geometry-dependent",
        strengths=("High paraffin regression rate", "Potential throttling", "Relatively dense fuel"),
        limitations=("Cryogenic oxidizer", "Combustion efficiency sensitivity", "Fuel mechanical design"),
        computability_note=(
            "Not computed: paraffin wax is not in NASA CEA thermo.inp and the hybrid O/F "
            "shifts with grain regression."
        ),
    ),
)


def list_propellants(category: str = "all") -> tuple[PropellantRecord, ...]:
    """Return the immutable propellant records, optionally filtered by category.

    Parameters
    ----------
    category : {'all', 'liquid', 'solid', 'hybrid'}, default 'all'
        Category filter; case and surrounding whitespace are ignored.

    Returns
    -------
    tuple of PropellantRecord
        Matching records in catalogue order.

    Raises
    ------
    DomainError
        If the category is not recognized.

    References
    ----------
    .. [1] NASA, "2020 NASA Technology Taxonomy", TX01 Propulsion Systems (liquid,
           solid, hybrid, cryogenic, storable and gelled categories).
    """

    normalized = category.strip().lower()
    if normalized == "all":
        return PROPELLANTS
    if normalized not in {"liquid", "solid", "hybrid"}:
        raise DomainError("Propellant category must be all, liquid, solid, or hybrid.")
    return tuple(record for record in PROPELLANTS if record.category == normalized)

