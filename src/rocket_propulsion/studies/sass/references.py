"""Evidence registry for the provisional UZM 508E propulsion study tools.

References
----------
NASA-SE-APP: requirement/source traceability, Appendix D.
ITU-UZM508: official course identity; it does not publish an assignment brief.
"""

from __future__ import annotations

from dataclasses import dataclass

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class Reference:
    """A source and the location actually used, rather than a generic bibliography.

    References
    ----------
    NASA-SE-APP: Appendix D, identify the definitive requirement source.
    """

    key: str
    title: str
    publisher: str
    url: str
    locator: str
    accessed: str = "2026-10-18"


REFERENCES = (
    Reference(
        "NASA-TRL", "Technology Readiness Levels", "NASA ESTO",
        "https://esto.nasa.gov/trl/",
        "TRL 5 hardware definition and exit criterion: relevant environment and analytical agreement.",
    ),
    Reference(
        "NASA-7009B", "Standard for Models and Simulations, NASA-STD-7009B", "NASA",
        "https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf",
        "Approved 2024-03-05; 4.2.3-4.2.7, verification, validation and uncertainty; 4.3.8 reporting.",
    ),
    Reference(
        "SAMSON-2019", "Fabrication and Testing of the Cold Gas Propulsion System Flight Unit "
        "for the Adelis-SAMSON Nano-Satellites", "Zaberchik, Lev, Edlerman and Kaidar; Aerospace",
        "https://www.mdpi.com/2226-4310/6/8/91",
        "2019, 6(8), 91; DOI 10.3390/aerospace6080091; component/system tests, not our hardware.",
    ),
    Reference(
        "NASA-GFSSP6", "Generalized Fluid System Simulation Program, Version 6.0", "NASA",
        "https://ntrs.nasa.gov/citations/20160005401",
        "NASA/TP-2016-218218; analytical and experimental examples need case-specific pedigree review.",
    ),
    Reference(
        "ITU-UZM508", "UZM 508E - Spacecraft Architecture and Subsystems Design",
        "Istanbul Technical University",
        "https://ninova.itu.edu.tr/en/courses/faculty-of-aeronautics-and-astronautics/37923/uzm-508e/",
        "Course title and coordinator; public description contains only 'Sassd'.",
    ),
    Reference(
        "ITU-UZM-GROUP", "Graduate course group: Zorunlu Secmeli Ders III", "ITU",
        "https://obs.itu.edu.tr/public/DersPlan/_DersGrupSearch?grupId=13265",
        "UZM 508E row: 3 credits, 7.5 ECTS, 3 theoretical hours.",
    ),
    Reference(
        "ITU-UZB419", "Spacecraft Systems Design I course catalog", "ITU",
        "https://obs.itu.edu.tr/public/DersKatalog/DersKatalogBilgiBransDersKodu?bransKodu=UZB&dersNo=419",
        "Related undergraduate course ONLY: objectives, course plan weeks 1-5 and final project.",
    ),
    Reference(
        "NASA-SE-APP", "Systems Engineering Handbook: Appendix", "NASA",
        "https://www.nasa.gov/reference/system-engineering-handbook-appendix/",
        "Appendix B: Trade Study Report; Appendix D: Requirements Verification Matrix.",
    ),
    Reference(
        "NASA-SE-6", "Systems Engineering Handbook: Crosscutting Technical Management", "NASA",
        "https://www.nasa.gov/reference/6-0-crosscutting-technical-management/",
        "Sections 6.3 (interfaces), 6.6 (technical data), 6.8 (decision analysis).",
    ),
    Reference(
        "NASA-ROCKET", "Ideal Rocket Equation", "NASA Glenn Research Center",
        "https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/",
        "Delta u = Isp g0 ln(MR), inverse mass ratio; ideal external-force-free assumptions.",
    ),
    Reference(
        "NASA-ISP", "Specific Impulse", "NASA Glenn Research Center",
        "https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/",
        "Total Impulse and Specific Impulse sections: I = integral F dt; Isp = F/(mdot g0).",
    ),
    Reference(
        "NASA-GAS", "Equation of State (Ideal Gas)", "NASA Glenn Research Center",
        "https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/",
        "Deriving an Equation: pV=nRT and p=rho R T; absolute temperature and ideal-gas limit.",
    ),
    Reference(
        "NASA-SP8112", "Pressurization Systems for Liquid Rockets", "NASA",
        "https://ntrs.nasa.gov/citations/19760015212",
        "NASA SP-8112 (1975): preliminary blowdown and regulated pressurization design.",
    ),
    Reference(
        "MIT-522-L5", "16.522 Session 5: Review of Classical Astrodynamics", "MIT OCW",
        "https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf",
        "Equation (11), specific orbital energy; equation (19), thrust acceleration F/m.",
    ),
    Reference(
        "MIT-522-L6", "16.522 Session 6: Analytical Approximations for Low Thrust Maneuvers",
        "MIT OCW",
        "https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/7f725e54b9be201164d56ebbd5e08023_MIT16_522S15_Lecture6.pdf",
        "Pages 1-2, equations (4)-(8): near-circular tangential spiral and Hohmann transfer.",
    ),
    Reference(
        "MIT-512-L35", "16.512 Rocket Propulsion, Lectures 35-36", "MIT OCW",
        "https://ocw.mit.edu/courses/16-512-rocket-propulsion-fall-2005/b0b082753985a3c02eb408af4d5a46a4_lecture_35_36.pdf",
        "Page 5, equation (2): transfer coast time is half the ellipse period, pi*a^(3/2)/sqrt(mu).",
    ),
    Reference(
        "NASA-SOA26", "State-of-the-Art of Small Spacecraft Technology: In-Space Propulsion",
        "NASA",
        "https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/",
        "2026 web edition: Table 4-1 technology envelopes; ECAPS 1 N HPGP device row.",
    ),
    Reference(
        "MOOG-COLD", "Cold Gas Thrusters: Performance Characteristics", "Moog",
        "https://www.moog.com/products/propulsion/space-propulsion/spacecraft-propulsion/thrusters/cold-gas-thruster.html",
        "Model 058E151 column: nominal thrust/Isp/inlet pressure; thruster mass and valve power.",
    ),
    Reference(
        "ARIANE-1N", "1 N Monopropellant Hydrazine Thruster", "ArianeGroup",
        "https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html",
        "1N Monopropellant Thruster Characteristics: nominal thrust, Isp, reported flow and mass.",
    ),
    Reference(
        "ARIANE-10N", "10 N Bipropellant Thruster", "ArianeGroup",
        "https://www.space-propulsion.com/spacecraft-propulsion/bipropellant-thrusters/10-bipropellant-thrusters.html",
        "10N characteristics: nominal thrust/Isp, nominal mixture-ratio interval, valve-specific mass.",
    ),
    Reference(
        "BUSEK-BIT3", "BIT-3 RF Ion Thruster product page", "Busek",
        "https://www.busek.com/bit3",
        "Performance panel: 56-75 W, up to 1.1 mN, up to 2150 s, system dry mass 1.40 kg with gimbal.",
    ),
    Reference(
        "BUSEK-RF", "RF Ion Thrusters product family page", "Busek",
        "https://www.busek.com/rf-ion-thrusters",
        "BIT-3 panel differs from /bit3: 56-80 W, up to 1.25 mN and 2300 s. No silent merge.",
    ),
    Reference(
        "BIPM-SI", "SI Brochure, 9th edition", "BIPM",
        "https://www.bipm.org/en/publications/si-brochure",
        "Conventional standard gravity and defining molar-gas constant; see core.constants.",
    ),
)


def reference(key: str) -> Reference:
    """Resolve a reviewed source; unknown identifiers are errors.

    References
    ----------
    NASA-SE-APP: Appendix D, definitive source identification.
    """
    for item in REFERENCES:
        if item.key == key:
            return item
    raise DomainError(f"Unknown SASS reference: {key!r}.")
