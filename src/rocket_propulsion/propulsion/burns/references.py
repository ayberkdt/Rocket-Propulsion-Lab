"""Primary-source registry for propulsion-burn equations and model claims.

The registry keeps citation identity stable across Python docstrings, JSON
artifacts, reports, and the future Sidera adapter.  URLs point to the issuing
organization rather than a secondary summary.
"""

from __future__ import annotations

from dataclasses import dataclass

from rocket_propulsion.core.errors import InputError


@dataclass(frozen=True, slots=True)
class EngineeringReference:
    """One immutable primary engineering reference.

    References
    ----------
    NASA Scientific and Technical Information program, NTRS metadata and
    public technical-report records: https://ntrs.nasa.gov/
    """

    reference_id: str
    title: str
    organization: str
    year: int
    url: str
    scope: str


NASA_SP_125 = EngineeringReference(
    "NASA-SP-125",
    "Design of Liquid Propellant Rocket Engines, Second Edition",
    "NASA",
    1967,
    "https://ntrs.nasa.gov/citations/19710019929",
    "Liquid-engine system design, thrust chamber, feed system, and performance foundations.",
)

NASA_SP_8112 = EngineeringReference(
    "NASA-SP-8112",
    "Pressurization Systems for Liquid Rockets",
    "NASA",
    1975,
    "https://ntrs.nasa.gov/citations/19760015212",
    "Blowdown-system architecture, ullage expansion, pressure decay, and design limits.",
)

NASA_SP_8080 = EngineeringReference(
    "NASA-SP-8080",
    "Liquid Rocket Pressure Regulators, Relief Valves, Check Valves, Burst Disks, and Explosive Valves",
    "NASA",
    1973,
    "https://ntrs.nasa.gov/citations/19740002611",
    "Pressure-regulator and valve design criteria for liquid rocket pressurization systems.",
)

NASA_NODAL_FEED_PRESSURIZATION = EngineeringReference(
    "NASA-20240003493",
    "Nodal Modeling of Liquid Propellant Feed and Pressurization System",
    "NASA",
    2024,
    "https://ntrs.nasa.gov/citations/20240003493",
    "Mass, energy, equation-of-state, and branch-flow modeling of coupled feed and pressurization networks.",
)

NASA_LIQUID_FEEDLINE_DYNAMICS = EngineeringReference(
    "NASA-19740028545",
    "Liquid Rocket Propellant Feedline Dynamics",
    "NASA",
    1973,
    "https://ntrs.nasa.gov/citations/19740028545",
    "Analytical and experimental pressure/flow dynamics of compliant liquid-rocket feed lines.",
)

NASA_NESC_TRANSIENT_PRESSURE = EngineeringReference(
    "NASA-20220006583",
    "NESC Technical Bulletin 22-03: Treatment of Transient Pressure Events in Space Flight Pressurized Systems",
    "NASA",
    2022,
    "https://ntrs.nasa.gov/citations/20220006583",
    "Valve-actuation, priming, discharge, vibration, and flow-disturbance pressure transients.",
)

NASA_GRC_MASS_FLOW_CHOKING = EngineeringReference(
    "NASA-GRC-MASS-FLOW-CHOKING",
    "Mass Flow Choking",
    "NASA Glenn Research Center",
    2021,
    "https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html",
    "Compressible mass-flow function and sonic choking limit for an ideal-gas restriction.",
)

NASA_CR_131400 = EngineeringReference(
    "NASA-CR-131400",
    "Reliability Model of a Monopropellant Auxiliary Propulsion System",
    "NASA",
    1971,
    "https://ntrs.nasa.gov/citations/19730012094",
    "Time-history modeling of an unregulated blowdown hydrazine propulsion system.",
)

NASA_CR_122347 = EngineeringReference(
    "NASA-CR-122347",
    "Monopropellant Hydrazine Resisto Jet",
    "NASA",
    1971,
    "https://ntrs.nasa.gov/citations/19720011120",
    "Supply-pressure effects on thrust and Isp, pulse operation, rise, and tailoff.",
)

NASA_TM_107318 = EngineeringReference(
    "NASA-TM-107318",
    "RL10A-3-3A Rocket Engine Modeling Project",
    "NASA",
    1997,
    "https://ntrs.nasa.gov/citations/19970010379",
    "Component-model assembly and comparison of start/steady/shutdown predictions with data.",
)

NASA_TRANSIENT_2004 = EngineeringReference(
    "NASA-20040000363",
    "Transient Mathematical Modeling for Liquid Rocket Engine Systems",
    "NASA",
    2003,
    "https://ntrs.nasa.gov/citations/20040000363",
    "Governing-equation, numerical, test-correlation, and validity guidance for transients.",
)

NASA_CR_140800 = EngineeringReference(
    "NASA-CR-140800",
    "Investigation of the Starting Transients of High Performance Solid-Propellant Motors",
    "NASA",
    1974,
    "https://ntrs.nasa.gov/citations/19750003988",
    "Solid-motor ignition transients and the need to preserve measured transient histories.",
)

NASA_THRUST_EQUATION = EngineeringReference(
    "NASA-GRC-THRUST-EQUATION",
    "Rocket Thrust Equation",
    "NASA Glenn Research Center",
    2024,
    "https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/",
    "Momentum and pressure-thrust identity for liquid and solid rockets.",
)

CCSDS_502_0_B_3 = EngineeringReference(
    "CCSDS-502.0-B-3",
    "Orbit Data Messages, Recommended Standard",
    "CCSDS",
    2023,
    "https://ccsds.org/Pubs/502x0b3e1.pdf",
    "Versioned, unit-explicit interchange and maneuver metadata context.",
)

NASA_STD_7009B = EngineeringReference(
    "NASA-STD-7009B",
    "Standard for Models and Simulations",
    "NASA",
    2024,
    "https://standards.nasa.gov/standard/NASA/NASA-STD-7009",
    "Model credibility, verification, validation, uncertainty qualification, and use assessment.",
)

NASA_HDBK_7009 = EngineeringReference(
    "NASA-HDBK-7009",
    "NASA Handbook for Models and Simulations: An Implementation Guide for NASA-STD-7009",
    "NASA",
    2013,
    "https://ntrs.nasa.gov/citations/20140002378",
    "Validation evidence, model discrepancy, uncertainty characterization, and use assessment.",
)

NASA_SP_2011_3421 = EngineeringReference(
    "NASA-SP-2011-3421",
    "Probabilistic Risk Assessment Procedures Guide for NASA Managers and Practitioners",
    "NASA",
    2011,
    "https://ntrs.nasa.gov/citations/20120001369",
    "Uncertainty characterization, Monte Carlo/Latin-hypercube propagation, and sensitivity.",
)

NASA_SP_2009_569 = EngineeringReference(
    "NASA-SP-2009-569",
    "Bayesian Inference for NASA Probabilistic Risk and Reliability Analysis",
    "NASA",
    2009,
    "https://ntrs.nasa.gov/citations/20090023159",
    "Separation and treatment of data, parameter uncertainty, and model knowledge.",
)

NASA_MULTIPLE_VALIDATION_CALIBRATION = EngineeringReference(
    "NASA-20150006032",
    "Calibration of Predictor Models Using Multiple Validation Experiments",
    "NASA",
    2015,
    "https://ntrs.nasa.gov/citations/20150006032",
    "Calibration from multiple validation experiments and separation of measurement and model-form uncertainty.",
)

NASA_AEROSPIKE_PARAMETRIC_MODEL = EngineeringReference(
    "NASA-20000031654",
    "Parametric Model of an Aerospike Rocket Engine",
    "NASA",
    2000,
    "https://ntrs.nasa.gov/citations/20000031654",
    "Trajectory-coupled engine performance tables over mixture ratio, power level, vectoring, and altitude.",
)

NASA_LIQUID_ENGINE_RELIABILITY = EngineeringReference(
    "NASA-20050207429",
    "Key Reliability Drivers of Liquid Propulsion Engines and a Reliability Model for Sensitivity Analysis",
    "NASA",
    2005,
    "https://ntrs.nasa.gov/citations/20050207429",
    "Engine-out, start/cutoff transient, duration, restart, and propulsion reliability drivers.",
)

NASA_COMMON_CAUSE_FAILURE = EngineeringReference(
    "NASA-20160007073",
    "Common Cause Failure Modeling in Space Launch Vehicles",
    "NASA",
    2016,
    "https://ntrs.nasa.gov/citations/20160007073",
    "Dependent failures that can defeat redundant launch-vehicle systems and their data limits.",
)

NIST_TN_1297 = EngineeringReference(
    "NIST-TN-1297",
    "Guidelines for Evaluating and Expressing the Uncertainty of NIST Measurement Results",
    "NIST",
    1994,
    "https://www.nist.gov/pml/nist-technical-note-1297",
    "Type A/Type B measurement uncertainty, combination, and reporting guidance.",
)

NIST_VARIANCE_COMPONENTS = EngineeringReference(
    "NIST-SEMATECH-7.4.4",
    "NIST/SEMATECH e-Handbook of Statistical Methods: Variance Components",
    "NIST",
    2013,
    "https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm",
    "One-way random-effects models and within-group/between-group variance components.",
)

NIST_CONSENSUS_MEAN = EngineeringReference(
    "NIST-DATAPLOT-CONSENSUS-MEAN",
    "NIST Dataplot: Consensus Mean",
    "NIST",
    2025,
    "https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm",
    "Heteroscedastic one-way random-effects models and Mandel-Paule estimation.",
)

JCGM_100_2008 = EngineeringReference(
    "JCGM-100-2008",
    "Evaluation of Measurement Data — Guide to the Expression of Uncertainty in Measurement",
    "JCGM",
    2008,
    "https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf",
    "Propagation of correlated measurement uncertainty and covariance reporting.",
)

NASA_NEXT_IPS_SUMMARY = EngineeringReference(
    "NASA-20090004685",
    "NASA's Evolutionary Xenon Thruster (NEXT) Ion Propulsion System Information Summary",
    "NASA",
    2008,
    "https://ntrs.nasa.gov/citations/20090004685",
    "Electric-thruster throttle tables: input power, thrust, specific impulse and efficiency.",
)

OREKIT_PROPULSION_MODEL = EngineeringReference(
    "OREKIT-13.1.5-PROPULSION-MODEL",
    "Orekit 13.1.5 PropulsionModel API",
    "Orekit Project",
    2026,
    "https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html",
    "Propagation-time acceleration, mass derivative, parameter-driver, and event-detector contract.",
)


BURN_REFERENCES: tuple[EngineeringReference, ...] = (
    CCSDS_502_0_B_3,
    JCGM_100_2008,
    NASA_CR_122347,
    NASA_CR_131400,
    NASA_CR_140800,
    NASA_COMMON_CAUSE_FAILURE,
    NASA_AEROSPIKE_PARAMETRIC_MODEL,
    NASA_LIQUID_ENGINE_RELIABILITY,
    NASA_LIQUID_FEEDLINE_DYNAMICS,
    NASA_HDBK_7009,
    NASA_GRC_MASS_FLOW_CHOKING,
    NASA_MULTIPLE_VALIDATION_CALIBRATION,
    NASA_NODAL_FEED_PRESSURIZATION,
    NASA_NESC_TRANSIENT_PRESSURE,
    NASA_NEXT_IPS_SUMMARY,
    NASA_SP_125,
    NASA_SP_2009_569,
    NASA_SP_2011_3421,
    NASA_SP_8112,
    NASA_SP_8080,
    NASA_STD_7009B,
    NASA_THRUST_EQUATION,
    NASA_TM_107318,
    NASA_TRANSIENT_2004,
    NIST_CONSENSUS_MEAN,
    NIST_TN_1297,
    NIST_VARIANCE_COMPONENTS,
    OREKIT_PROPULSION_MODEL,
)


def burn_reference(reference_id: str) -> EngineeringReference:
    """Return one registered source by stable identifier.

    Raises
    ------
    InputError
        If the identifier is not in the reviewed burn reference set.

    References
    ----------
    NASA STI/NTRS catalog: https://ntrs.nasa.gov/
    CCSDS publication catalog: https://ccsds.org/publications/
    """

    for reference in BURN_REFERENCES:
        if reference.reference_id == reference_id:
            return reference
    raise InputError(f"Unknown burn reference identifier: {reference_id}.")
