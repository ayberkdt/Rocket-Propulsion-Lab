"""Curated screening records for widely used chemical rocket propellant systems.

Notes
-----
Vacuum Isp entries are deliberately broad preliminary-design ranges for each
propellant family [1]_ [2]_; categories follow the liquid/solid/hybrid
classification of [3]_. They are not substitutes for equilibrium chemistry, engine
test data or a validated cycle/nozzle model, and mixture-ratio notes describe
screening regions, not operating instructions.

References
----------
.. [1] NASA KSC Engineering, "Examples of Specific Impulse", NTRS 20100036816
       (representative LOX/LH2, LOX/RP-1, storable and solid systems).
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

    Notes
    -----
    Vacuum Isp entries are deliberately broad preliminary-design ranges for each
    propellant family [1]_ [2]_; categories follow the liquid/solid/hybrid
    classification of [3]_. They are not substitutes for equilibrium chemistry, engine
    test data or a validated cycle/nozzle model, and mixture-ratio notes describe
    screening regions, not operating instructions.

    References
    ----------
    .. [1] NASA KSC Engineering, "Examples of Specific Impulse", NTRS 20100036816
           (representative LOX/LH2, LOX/RP-1, storable and solid systems).
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

