# Engineering references

- **SETZMANN-WAGNER-1991-CH4-NIST-MIRROR**: Official NIST-hosted copy of
  Setzmann and Wagner's 1991 methane reference EOS paper; the model's original
  NIST DOI/PDF reference remains registered separately. Accessed 2026-10-29.
  <https://www.nist.gov/system/files/documents/srd/jpcrd425.pdf>

- **PYTHON-3.12-LRU-CACHE**: Bounded least-recently-used function memoization.
  <https://docs.python.org/3.12/library/functools.html#functools.lru_cache>

- **PYTHON-3.12-INTEGER-LIMIT**: Integer/string conversion limits and default policy.
  <https://docs.python.org/3.12/library/stdtypes.html#integer-string-conversion-length-limitation>
- **PYTHON-3.12-DECIMAL-ADJUSTED**: Decimal adjusted exponent.
  <https://docs.python.org/3.12/library/decimal.html#decimal.Decimal.adjusted>

- **HALL-MCCABE-RM3440**: McCabe, planar Hall coefficients and nozzle-design context, ARC R&M3440 (report 1964, issued 1967). <https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023>
- **KLIEGEL-QUAN-1966**: Full potential equations, Cartesian Hall coefficients and comparison of transonic approximations, 02874-6002-R000. Theory context only. <https://ntrs.nasa.gov/api/citations/19670009614/downloads/19670009614.pdf>
- **SAUER-NACA-TM1147**: Historical near-critical-flow and curved-sonic-line theory, 1947. <https://ntrs.nasa.gov/api/citations/20050019415/downloads/20050019415.pdf>

- **NASA-GLENN-MASS-FLOW**: Isentropic mass-flow conservation and density/velocity context. Used as theory context, not a validation dataset or hardware evidence. Accessed 2026-10-24. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/>

- **PYTHON-3.12-JSON**: Python JSON encoder/decoder hooks, duplicate-member and nonstandard numeric token behavior. Accessed 2026-10-19. <https://docs.python.org/3.12/library/json.html>

- **SCIPY-1.14.1-REGULAR-GRID**: SciPy, RegularGridInterpolator, rectilinear grids and linear interpolation. Accessed 2026-10-18. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.interpolate.RegularGridInterpolator.html>

- **SCIPY-1.14.1-PCHIP**: SciPy, PchipInterpolator algorithm and implementation. Accessed 2026-10-18. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.interpolate.PchipInterpolator.html>, <https://raw.githubusercontent.com/scipy/scipy/v1.14.1/scipy/interpolate/_cubic.py>.

- **PYTHON-3.12-HYPOT**: Python 3.12, `math.hypot`, Euclidean norm and floating-point accuracy. Accessed 2026-10-18. <https://docs.python.org/3.12/library/math.html#math.hypot>

- Python 3.12 faulthandler and CPython 3.12.1 implementation (test diagnostics):
  <https://docs.python.org/3.12/library/faulthandler.html>,
  <https://raw.githubusercontent.com/python/cpython/v3.12.1/Modules/faulthandler.c>.
- CPython native traceback race reports (diagnostic context, not propulsion validation):
  <https://github.com/python/cpython/issues/116008>,
  <https://github.com/python/cpython/issues/131580>.

- Python 3.12 threading: Thread, Lock, Event and daemon shutdown semantics,
  <https://docs.python.org/3.12/library/threading.html> (HTTP capability probe
  implementation only; no propulsion or physical validation claim).

- NASA Glenn, **Specific Heats**, ideal-gas cp=cv+R and gamma=cp/cv. <https://www.grc.nasa.gov/www/k-12/airplane/specheat.html>

- **NASA-TRL**: NASA ESTO, hardware TRL definitions and exit criteria. <https://esto.nasa.gov/trl/>
- **NASA-7009B**: NASA-STD-7009B, approved 2024-03-05; verification, empirical validation, uncertainty and reporting. <https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf>
- **SAMSON-2019**: Zaberchik et al., Aerospace 2019, 6(8), 91; measured cold-gas performance and component/system testing of their hardware. <https://www.mdpi.com/2226-4310/6/8/91>
- **NASA-GFSSP6**: NASA/TP-2016-218218, fluid network examples and comparisons; not an executed comparison against this library. <https://ntrs.nasa.gov/citations/20160005401>
- NASA/JSC, 20170000809, one-page cold-gas system abstract; not raw data. <https://ntrs.nasa.gov/citations/20170000809>
- DLR, DX2025 benchmark, section 4.3 distinguishes simulated trajectories from experimental measurements. <https://elib.dlr.de/219953/1/DX2025benchmark.pdf>

The equations in this repository are standard perfect-gas and ideal rocket
relations. These primary NASA resources anchor terminology and screening values:

- NASA Glenn, **NASA Glenn Coefficients for Calculating Thermodynamic Properties
  of Individual Species**, NASA/TP-2002-211556 — seven-term heat-capacity
  polynomials, enthalpy/entropy integration constants, temperature intervals,
  and the source species database:
  <https://ntrs.nasa.gov/citations/20020085330>

  Equations (1)--(3), printed p.2; Appendix A, printed p.7, defines the fit's
  R=8.314510 J/(mol K); Appendix C gives fixed-column data layout. The official
  NASA/CEA subset comparison is pinned to commit
  `985fd33368b2334c4162a586c499e0ab41d12bd3`; its `data/thermo.inp` has the same
  SHA-256 as the NASA-CEA-THERMO pin below:
  <https://github.com/nasa/cea/blob/985fd33368b2334c4162a586c499e0ab41d12bd3/data/thermo.inp>

- NASA Glenn, **Ideal Rocket Equation** — Tsiolkovsky equation, mass ratio,
  effective exhaust velocity, and specific impulse:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>
- NASA Glenn, **Rocket Thrust Equation** — momentum thrust, pressure thrust,
  choked throat, and exit-area effects:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/>
- NASA Glenn, **Specific Impulse** — equivalent velocity, total impulse, and Isp:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/>
- NASA Technology Taxonomy, **TX01 Propulsion Systems** — liquid, solid, hybrid,
  cryogenic, storable, and gelled propulsion categories:
  <https://www.nasa.gov/wp-content/uploads/2015/03/2020_nasa_technology_taxonomy_lowres.pdf>
- NASA KSC Engineering, **Examples of Specific Impulse** — representative
  comparison of LOX/LH2, LOX/RP-1, storable, and solid systems:
  <https://ntrs.nasa.gov/api/citations/20100036816/downloads/20100036816.pdf>
- NASA RP-1311, **Computer Program for Calculation of Complex Chemical
  Equilibrium Compositions and Applications, Part II** — CEA input/output,
  examples, and rocket problem definitions:
  <https://ntrs.nasa.gov/citations/19960044559>
- NASA RP-1311, **Part I: Analysis** — equilibrium formulation and
  thermodynamic methods: <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-SP8120-SEPARATION**: NASA SP-8120, **Liquid Rocket Engine Nozzles**, 1976 —
  printed p.17, §2.1.2.1.3: the early rule that separation is a danger when the
  exit-to-ambient pressure ratio equals 0.4, and the later wall-pressure fit;
  printed p.18, Figure 7: divergence efficiency against area ratio and half-angle.
  PDF pp.31-32 inspected 2026-10-25. <https://ntrs.nasa.gov/citations/19770009165>
- **NASA-TM-77396-SCHMUCKER**: R. H. Schmucker, **Flow Processes in Overexpanded
  Chemical Rocket Nozzles. Part 1: Flow Separation**, NASA TM-77396
  (TUM-LRT-TB-7 translation), 1984 — printed p.37, Eq. (7): Summerfield's
  `p_i/p_a = 0.4`, a conservative guideline for `p_c/p_a` 15 to 20; printed
  p.38, Eq. (11) and Fig. 25: `p_i/p_a = (1.88 M_i - 1)^(-0.64)` fitted to
  averaged separation data of chemical engines; printed p.39: not to be used
  above `M_i = 5`. PDF pp.43-45 inspected 2026-10-25.
  <https://ntrs.nasa.gov/citations/19840011402>
- Correction (2026-10-25): NASA SP-8041 is **Captive-Fired Testing of Solid Rocket
  Motors** (1971), <https://ntrs.nasa.gov/citations/19710021390>. Earlier versions
  of this library cited it, under a wrong title, for the nozzle separation rule;
  that rule is now cited to **NASA-SP8120-SEPARATION**.
- NASA Glenn, **Inlets and Nozzles Design and Analysis Software** — method of
  characteristics design-tool context:
  <https://www.nasa.gov/glenn/research/inlets-and-nozzles/design-analysis-software/>
- NASA SP-125, **Design of Liquid Propellant Rocket Engines** — Bartz
  heat-transfer use in preliminary liquid-engine design:
  <https://ntrs.nasa.gov/citations/19710019929>

## Physical constants

Registered with value, unit and source in `rocket_propulsion.core.constants`;
`tests/test_constants_and_references.py` refuses bare copies elsewhere.

- BIPM, **The International System of Units (SI Brochure)**, 9th ed., 2019 —
  exact \(N_A\) and \(k\) (so \(R = N_A k = 8.314\,462\,618\ldots\) J/(mol K)),
  the Celsius offset 273.15 K and the conventional
  \(g_n = 9.80665\) m/s²: <https://www.bipm.org/en/publications/si-brochure>
- NOAA/NASA/USAF, **U.S. Standard Atmosphere, 1976**, NOAA-S/T 76-1562 /
  NASA-TM-X-74335 — defining constants \(R^\ast = 8.31432\) J/(mol K),
  \(M_0 = 28.9644\) kg/kmol, \(P_0 = 101\,325\) Pa, \(T_0 = 288.15\) K and the
  layer table below 84.852 km: <https://ntrs.nasa.gov/citations/19770009539>

## Uncertainty and calibration references

- JCGM 100:2008, **Evaluation of Measurement Data — Guide to the Expression of
  Uncertainty in Measurement** — propagation of correlated uncertainty and
  covariance reporting:
  <https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf>
- NIST TN 1297, **Guidelines for Evaluating and Expressing the Uncertainty of
  NIST Measurement Results**, 1994 — Type A/Type B uncertainty:
  <https://www.nist.gov/pml/nist-technical-note-1297>
- NIST/SEMATECH e-Handbook §7.4.4, **Variance Components** — one-way
  random-effects models:
  <https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm>
- NIST Dataplot, **Consensus Mean** — heteroscedastic random effects and
  Mandel–Paule estimation:
  <https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm>
- NASA CR-184099 — transient engine simulation covering steady balance,
  throttle transients, starts and shutdowns (cited by the Sidera plan):
  <https://ntrs.nasa.gov/citations/19910011919>

## Numerical methods

- NIST DLMF §3.5(v), **Gauss quadrature** (Gauss–Legendre nodes and weights):
  <https://dlmf.nist.gov/3.5#v>
- F. N. Fritsch and R. E. Carlson, **Monotone Piecewise Cubic Interpolation**,
  SIAM J. Numer. Anal. 17(2), 1980; F. N. Fritsch and J. Butland, SIAM J. Sci.
  Stat. Comput. 5(2), 1984 — PCHIP slopes (bibliographic citations).
- C. L. Lawson and R. J. Hanson, **Solving Least Squares Problems**, 1974;
  P. B. Stark and R. L. Parker, **Bounded-Variable Least-Squares**, Comput.
  Stat. 10, 1995 — bounded thrust allocation (bibliographic citations).
- IETF RFC 8259, **The JavaScript Object Notation (JSON) Data Interchange
  Format** — the document grammar the readers enforce (no NaN/Infinity):
  <https://www.rfc-editor.org/rfc/rfc8259>
- D. Goldberg, **What Every Computer Scientist Should Know About
  Floating-Point Arithmetic**, ACM Computing Surveys 23(1), 1991 — the
  cancellation-free quadratic root (bibliographic citation).

## Finite-burn, transient, and blowdown references

- NASA, **Simulation Analysis of Computer-Controlled Pressurization for
  Mixture Ratio Control** — independent oxidizer/fuel flow measurement,
  pressure control, and mixture-ratio response:
  <https://ntrs.nasa.gov/citations/20050206420>
- NASA, **Space Simulation Experiments on Reaction Control System Thruster
  Plumes** — measured short-pulse bipropellant flow-ratio transients:
  <https://ntrs.nasa.gov/citations/19730002016>
- NASA SP-8059, **Spacecraft Attitude Control During Thrusting Maneuvers** —
  thrust-vector pointing, misalignment, centre-of-mass moment arms, and
  attitude-control design criteria:
  <https://ntrs.nasa.gov/citations/19710016722>
- NASA/TM-20240012525, **MSFC-TVC TB-03: Derivation of Thrust Vector Control
  Actuator-Force / Gimbal-Torque Transformation Matrix** — order-sensitive
  three-dimensional two-axis gimbal geometry:
  <https://ntrs.nasa.gov/citations/20240012525>
- NASA SP-8112, **Pressurization Systems for Liquid Rockets** — pressure-fed
  architectures, ullage expansion, blowdown pressure decay and sizing context:
  <https://ntrs.nasa.gov/citations/19760015212>
- NASA CR-131400, **Reliability Model of a Monopropellant Auxiliary Propulsion
  System** — performance time histories for an unregulated hydrazine blowdown
  system: <https://ntrs.nasa.gov/citations/19730012094>
- NASA CR-122347, **Monopropellant Hydrazine Resisto Jet** — measured
  supply-pressure effects on thrust/Isp plus pulse rise and tailoff behavior:
  <https://ntrs.nasa.gov/citations/19720011120>
- NASA TM-107318, **RL10A-3-3A Rocket Engine Modeling Project** — component
  transient model assembly and comparison with start/steady/shutdown data:
  <https://ntrs.nasa.gov/citations/19970010379>
- NASA, **Transient Mathematical Modeling for Liquid Rocket Engine Systems** —
  governing-equation, numerical and test-correlation guidance for transient
  simulations: <https://ntrs.nasa.gov/citations/20040000363>
- NASA CR-140800, **Investigation of the Starting Transients of High Performance
  Solid-Propellant Motors** — ignition transient and measured-history context:
  <https://ntrs.nasa.gov/citations/19750003988>
- CCSDS 502.0-B-3, **Orbit Data Messages** — versioned, unit-explicit exchange
  and maneuver-metadata context: <https://ccsds.org/Pubs/502x0b3e1.pdf>
- NASA-STD-7009B, **Standard for Models and Simulations** — model credibility,
  validation, sensitivity and uncertainty-qualification requirements:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
- NASA-HDBK-7009, **NASA Handbook for Models and Simulations** — implementation
  guidance for validation evidence, uncertainty characterization and intended
  use assessment: <https://ntrs.nasa.gov/citations/20140002378>
- NASA, **Calibration of Predictor Models Using Multiple Validation
  Experiments** — multiple validation cases, referent uncertainty and
  model-form uncertainty in predictor calibration:
  <https://ntrs.nasa.gov/citations/20150006032>
- NASA/SP-2011-3421, **Probabilistic Risk Assessment Procedures Guide for NASA
  Managers and Practitioners** — uncertainty characterization, Monte Carlo/LHS
  propagation and sensitivity context:
  <https://ntrs.nasa.gov/citations/20120001369>
- NASA/SP-2009-569, **Bayesian Inference for NASA Probabilistic Risk and
  Reliability Analysis** — data, parameter knowledge and reducible uncertainty
  treatment: <https://ntrs.nasa.gov/citations/20090023159>
- NASA, **Key Reliability Drivers of Liquid Propulsion Engines and a Reliability
  Model for Sensitivity Analysis** — engine count, engine-out design,
  start/cutoff transients, duration, restart and health-management reliability
  drivers: <https://ntrs.nasa.gov/citations/20050207429>
- NASA, **Common Cause Failure Modeling in Space Launch Vehicles** — dependent
  failures that defeat redundancy and limitations in applicable aerospace
  failure data: <https://ntrs.nasa.gov/citations/20160007073>
- Orekit 13.1.5, **PropulsionModel API** — the public maneuver-propulsion
  boundary for acceleration, mass derivatives, parameter drivers and model
  event detectors:
  <https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html>
- NASA, **Parametric Model of an Aerospike Rocket Engine** — trajectory-facing
  engine performance tables parameterized by mixture ratio, power level,
  thrust-vectoring level and altitude:
  <https://ntrs.nasa.gov/citations/20000031654>
- NASA SP-8080, **Liquid Rocket Pressure Regulators and Valves** — regulator,
  relief/check valve and pressure-system design criteria:
  <https://ntrs.nasa.gov/citations/19740002611>
- NASA, **Nodal Modeling of Liquid Propellant Feed and Pressurization System**
  — mass, energy, equation-of-state and branch-flow formulation for coupled
  feed networks: <https://ntrs.nasa.gov/citations/20240003493>
- NASA Glenn, **Mass Flow Choking** — ideal-gas compressible restriction flow
  and sonic maximum-flow relation:
  <https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html>
- NASA, **Liquid Rocket Propellant Feedline Dynamics** — liquid-feed system
  dynamic modes and trajectory-relevant feedline modeling context:
  <https://ntrs.nasa.gov/citations/19740028545>
- NASA NESC Technical Bulletin 22-03, **Treatment of Transient Pressure
  Events** — screening and design treatment for transient pressure and water
  hammer: <https://ntrs.nasa.gov/citations/20220006583>

## Propagation-interface references

- NASA, **NASA's Evolutionary Xenon Thruster (NEXT) Ion Propulsion System
  Information Summary** — electric-thruster throttle tables (input power,
  thrust, specific impulse, efficiency):
  <https://ntrs.nasa.gov/citations/20090004685>
- J. M. Parker and M. J. Wilson, **The Minimum Impulse Thruster**, 2005 —
  minimum impulse bit set by valve response and electrical pulse width
  (record located through the NTRS search index; the full text could not be
  opened from the development environment):
  <https://ntrs.nasa.gov/citations/20070032005>
- NASA/TP-2018-219822, **Navigation Filter Best Practices** — prior/consider
  covariance, Gauss-Markov process noise and maneuver execution errors:
  <https://ntrs.nasa.gov/citations/20180003657>
- Orekit 13.1.5 **PropulsionModel / EventDetector** API — the propagation-time
  force, mass-derivative, parameter-driver and event contract the bridge
  follows: <https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html>
- C. R. Gates, **A Simplified Model of Midcourse Maneuver Execution Errors**,
  JPL Technical Report 32-504, 1963 — original four-parameter execution-error
  model (bibliographic citation; no verified public URL recorded).
- B. Wie, **Space Vehicle Dynamics and Control**, 2nd ed., AIAA, 2008 —
  thruster torque and rigid-body mass properties (bibliographic citation).

Stable identifiers and scopes are defined in
`rocket_propulsion.propulsion.burns.references`. The transient public API is
tested for `References` sections; tabulated artifacts store reference IDs in
their hashed content.

## Interpretation rules

- Propellant Isp entries are deliberately broad preliminary-design ranges.
- They are not substitutes for equilibrium chemistry, engine test data, or a
  validated cycle/nozzle model.
- Mixture-ratio notes are descriptive screening regions, not operating instructions.
- All ideal models explicitly document omitted losses in their Python docstrings.
- Pressure-performance exponents are hardware-specific calibration inputs; the
  implementation does not transfer coefficients from a cited test article to
  unrelated engines and refuses pressure extrapolation.


## UZM 508E / SASS study sources (reviewed 2026-10-18)

Sources below support the provisional study package; ITU course identity is not an assignment brief. Product records retain source-specific variants and missing-data semantics.

- **ITU-UZM508**: Istanbul Technical University, *UZM 508E - Spacecraft Architecture and Subsystems Design*. Course title and coordinator; public description contains only 'Sassd'. Accessed 2026-10-18. <https://ninova.itu.edu.tr/en/courses/faculty-of-aeronautics-and-astronautics/37923/uzm-508e/>

- **ITU-UZM-GROUP**: ITU, *Graduate course group: Zorunlu Secmeli Ders III*. UZM 508E row: 3 credits, 7.5 ECTS, 3 theoretical hours. Accessed 2026-10-18. <https://obs.itu.edu.tr/public/DersPlan/_DersGrupSearch?grupId=13265>

- **ITU-UZB419**: ITU, *Spacecraft Systems Design I course catalog*. Related undergraduate course ONLY: objectives, course plan weeks 1-5 and final project. Accessed 2026-10-18. <https://obs.itu.edu.tr/public/DersKatalog/DersKatalogBilgiBransDersKodu?bransKodu=UZB&dersNo=419>

- **NASA-SE-APP**: NASA, *Systems Engineering Handbook: Appendix*. Appendix B: Trade Study Report; Appendix D: Requirements Verification Matrix. Accessed 2026-10-18. <https://www.nasa.gov/reference/system-engineering-handbook-appendix/>

- **NASA-SE-6**: NASA, *Systems Engineering Handbook: Crosscutting Technical Management*. Sections 6.3 (interfaces), 6.6 (technical data), 6.8 (decision analysis). Accessed 2026-10-18. <https://www.nasa.gov/reference/6-0-crosscutting-technical-management/>

- **NASA-ROCKET**: NASA Glenn Research Center, *Ideal Rocket Equation*. Delta u = Isp g0 ln(MR), inverse mass ratio; ideal external-force-free assumptions. Accessed 2026-10-18. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>

- **NASA-ISP**: NASA Glenn Research Center, *Specific Impulse*. Total Impulse and Specific Impulse sections: I = integral F dt; Isp = F/(mdot g0). Accessed 2026-10-18. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/>

- **NASA-GAS**: NASA Glenn Research Center, *Equation of State (Ideal Gas)*. Deriving an Equation: pV=nRT and p=rho R T; absolute temperature and ideal-gas limit. Accessed 2026-10-18. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/>

- **NASA-SP8112**: NASA, *Pressurization Systems for Liquid Rockets*. NASA SP-8112 (1975): preliminary blowdown and regulated pressurization design. Accessed 2026-10-18. <https://ntrs.nasa.gov/citations/19760015212>

- **MIT-522-L5**: MIT OCW, *16.522 Session 5: Review of Classical Astrodynamics*. Equation (11), specific orbital energy; equation (19), thrust acceleration F/m. Accessed 2026-10-18. <https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf>

- **MIT-522-L6**: MIT OCW, *16.522 Session 6: Analytical Approximations for Low Thrust Maneuvers*. Pages 1-2, equations (4)-(8): near-circular tangential spiral and Hohmann transfer. Accessed 2026-10-18. <https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/7f725e54b9be201164d56ebbd5e08023_MIT16_522S15_Lecture6.pdf>

- **MIT-512-L35**: MIT OCW, *16.512 Rocket Propulsion, Lectures 35-36*. Page 5, equation (2): transfer coast time is half the ellipse period, pi*a^(3/2)/sqrt(mu). Accessed 2026-10-18. <https://ocw.mit.edu/courses/16-512-rocket-propulsion-fall-2005/b0b082753985a3c02eb408af4d5a46a4_lecture_35_36.pdf>

- **NASA-SOA26**: NASA, *State-of-the-Art of Small Spacecraft Technology: In-Space Propulsion*. 2026 web edition: Table 4-1 technology envelopes; ECAPS 1 N HPGP device row. Accessed 2026-10-18. <https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/>

- **MOOG-COLD**: Moog, *Cold Gas Thrusters: Performance Characteristics*. Model 058E151 column: nominal thrust/Isp/inlet pressure; thruster mass and valve power. Accessed 2026-10-18. <https://www.moog.com/products/propulsion/space-propulsion/spacecraft-propulsion/thrusters/cold-gas-thruster.html>

- **ARIANE-1N**: ArianeGroup, *1 N Monopropellant Hydrazine Thruster*. 1N Monopropellant Thruster Characteristics: nominal thrust, Isp, reported flow and mass. Accessed 2026-10-18. <https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html>

- **ARIANE-10N**: ArianeGroup, *10 N Bipropellant Thruster*. 10N characteristics: nominal thrust/Isp, nominal mixture-ratio interval, valve-specific mass. Accessed 2026-10-18. <https://www.space-propulsion.com/spacecraft-propulsion/bipropellant-thrusters/10-bipropellant-thrusters.html>

- **BUSEK-BIT3**: Busek, *BIT-3 RF Ion Thruster product page*. Performance panel: 56-75 W, up to 1.1 mN, up to 2150 s, system dry mass 1.40 kg with gimbal. Accessed 2026-10-18. <https://www.busek.com/bit3>

- **BUSEK-RF**: Busek, *RF Ion Thrusters product family page*. BIT-3 panel differs from /bit3: 56-80 W, up to 1.25 mN and 2300 s. No silent merge. Accessed 2026-10-18. <https://www.busek.com/rf-ion-thrusters>

- **BIPM-SI**: BIPM, *SI Brochure, 9th edition*. Conventional standard gravity and defining molar-gas constant; see core.constants. Accessed 2026-10-18. <https://www.bipm.org/en/publications/si-brochure>

- **SAND2001-0417**: J. C. Helton and F. J. Davis, **Latin Hypercube Sampling and the Propagation of Uncertainty in Analyses of Complex Systems**, SAND2001-0417, Sandia National Laboratories, 2002 — Sect. 6: rank correlation (RCC), partial rank correlation (PRCC) and standardized rank regression coefficients (SRRC), rank-regression R^2, and RCC = SRRC for independent inputs (PDF pp. 43-44 inspected 2026-10-26). <https://doi.org/10.2172/806696>. Journal version (Crossref record checked): Reliability Engineering & System Safety 81(1), 23-69, 2003, doi:10.1016/S0951-8320(03)00058-9.
- **SALTELLI-2010-CPC**: A. Saltelli, P. Annoni, I. Azzini, F. Campolongo, M. Ratto and S. Tarantola, **Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index**, Computer Physics Communications 181(2), 259-270, 2010, doi:10.1016/j.cpc.2009.09.018 (bibliographic record via Crossref; publisher page not accessible). Estimators and output centring as implemented in SALib `analyze/sobol.py` (`first_order`, `total_order`, `Y = (Y - Y.mean()) / Y.std()`), read 2026-10-26: <https://github.com/SALib/SALib/blob/main/src/SALib/analyze/sobol.py>
- **SCIPY-1.14.1-FMINBOUND**: SciPy, *fminbound*: bounded scalar minimization by Brent's method (golden-section steps safeguarding parabolic interpolation), default xtol 1e-5; the stopping rule transcribed in `rocket_performance._maximize_bounded` follows `scipy/optimize/_optimize.py` (`_minimize_scalar_bounded`). Accessed 2026-10-26. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.optimize.fminbound.html>
- **SCIPY-1.14.1-IVP**: SciPy, *solve_ivp*. Event sign-change detection, terminal and direction attributes; multiple roots within one step may be missed. Accessed 2026-10-18. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.solve_ivp.html>

- **PYTHON-3.12-FRACTION**: Python, *fractions*. Float inputs retain their exact binary value as a rational. Accessed 2026-10-18. <https://docs.python.org/3.12/library/fractions.html>
- **PYTHON-3.12-NEXTAFTER**: Python, *math.nextafter*. Adjacent representable floating-point values for the direct inventory helper's safe time horizon. Accessed 2026-10-19. <https://docs.python.org/3.12/library/math.html#math.nextafter>
- **PYTHON-3.12-DECIMAL**: Python, *decimal*. Precision contexts and exact float conversion for independent numerical checks. Accessed 2026-10-18. <https://docs.python.org/3.12/library/decimal.html>

- **PYTHON-3.12-EXCEPTIONS**: Python, *Built-in Exceptions*. OverflowError ("Raised when the result of an arithmetic operation is too large to be represented"; most floating-point operations are not checked) and ZeroDivisionError (second argument of a division is zero). Used for `core.validation.numeric_range_guard`. Accessed 2026-10-25. <https://docs.python.org/3.12/library/exceptions.html>
- **PYTHON-3.12-NUMBERS**: Python, *numbers — Numeric abstract base classes*. Real, Rational and Integral protocols; application-level physical exclusions remain local contracts. Accessed 2026-10-20. <https://docs.python.org/3.12/library/numbers.html>
- **NASA-CR-2004-213099**: NASA, *Development and Characterization of High-Efficiency, High-Specific Impulse Xenon Hall Thrusters*. Printed page 14, equations (2-9)-(2-10): jet power and system efficiency referenced to spacecraft bus input power. Used for equations, not as validation data for an unselected motor. Accessed 2026-10-20. <https://ntrs.nasa.gov/api/citations/20040084644/downloads/20040084644.pdf>

- **NIST-WEIGHTED-MEAN**: NIST Dataplot, *Weighted Mean*, equation 2-20, normalized weighted mean. Accessed 2026-10-21. <https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weigmean.pdf>
- **MIT-TOTAL-VARIANCE**: MIT RES.6-012 Spring 2018, lecture 13.7 derivation transcript, conditional and between-group variance identity. Accessed 2026-10-21. <https://ocw.mit.edu/courses/res-6-012-introduction-to-probability-spring-2018/c025b4baa48ff7801b96c3915ee62e5c_mHonq7Gjjqg.pdf>
- **PYTHON-3.12-MATH**: Python 3.12 math, isqrt and ldexp primitives for local binary64 statistical exports. Accessed 2026-10-21. <https://docs.python.org/3.12/library/math.html>
- **NIST-DLMF-7.6**: NIST DLMF, positive error-function power series, equation 7.6.2. Accessed 2026-10-21. <https://dlmf.nist.gov/7.6#E2>
- **NIST-DLMF-7.8**: NIST DLMF, error-function and Mills ratio bounds. Accessed 2026-10-21. <https://dlmf.nist.gov/7.8>
- **NIST-TRIPPF**: NIST Dataplot, triangular percent point function. Accessed 2026-10-21. <https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/trippf.pdf>
- **SCIPY-1.14.1-LOGUNIFORM**: SciPy, log-uniform CDF and quantile formulas. Accessed 2026-10-21. <https://docs.scipy.org/doc/scipy-1.14.1/tutorial/stats/continuous_loguniform.html>
- **NIST-MEAN-INTERVAL**: NIST, confidence limits for the mean and long-run interpretation. Accessed 2026-10-21. <https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm>
- **NIST-KNOWN-SIGMA-INTERVAL**: NIST, known and estimated population deviation in confidence intervals. Accessed 2026-10-21. <https://www.itl.nist.gov/div898/handbook/prc/section1/prc14.htm>
- **NUMPY-2.2-LINEAR-QUANTILE**: NumPy, linear/Hyndman-Fan type-7 sample quantile definition. Accessed 2026-10-21. <https://numpy.org/doc/2.2/reference/generated/numpy.quantile.html>
- **NIST-PEARSON**: NIST Dataplot, centered-product correlation definition. Accessed 2026-10-21. <https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/correlat.htm>
- **SCIPY-1.14.1-PEARSON**: SciPy, Pearson coefficient and undefined constant-input behavior. Accessed 2026-10-21. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.stats.pearsonr.html>

- **NASA-STEADY-ROCKET-THRUST**: NASA Glenn, momentum and pressure thrust equation. Accessed 2026-10-21. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/>

- **NASA-TMX64749-CYCLE**: R. H. Schmucker, 1973, equivalent cycle power balances and whole pumped mass flows. Theory source, not selected hardware validation. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19730016059>
- **NASA-ENGINESIM-THEORY**: NASA Glenn, ideal gas turbine temperature/pressure and shaft-work relations. Accessed 2026-10-24. <https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf>
- **NASA-SP8120-DISCHARGE**: NASA SP-8120, Liquid Rocket Engine Nozzles, 1976, discharge coefficient and flow effects, printed pp.10,19-20,93. Accessed 2026-10-24. <https://ntrs.nasa.gov/api/citations/19770009165/downloads/19770009165.pdf>
- **ROCKETISP-CSTAR-CD-CONVENTION**: RocketIsp 0.1.13 project documentation, its declared intrinsic cstar/discharge relation; not independent NASA qualification. Accessed 2026-10-24. <https://rocketisp.readthedocs.io/en/latest/efficiencies.html>
- **NIST-TN1297-COVARIANCE**: NIST TN1297 Appendix A, first-order measurement uncertainty with sensitivity and covariance. Accessed 2026-10-24. <https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-law-propagation-uncertainty>
- **DREWNIAK-TWO-PHASE-2022**: Drewniak, Boyer and Yetter, 2022, Eq.(3), separate gas/particle momentum and pressure. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/20220006969>
- **KIBBEY-TWO-PHASE-2012**: Kibbey, 2012, coupled gas/particle flow and momentum debit. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/20120015678>
- **WHITESIDES-COLD-FLOW-1995**: Whitesides et al., NASA-CR-204231/AIAA95-2724, Table I, air/water pressure-response experiment; calibration is not independent validation. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19970015335>
- **NASA-SP125-BARTZ-PAGES**: SP-125, printed pp.100-101, equations 4-10, 4-13 and 4-14; original PDF pages 109-110 visually inspected. Accessed 2026-10-24. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- **NASA-HOPPER-2012-THERMAL**: Linne et al., AIAA 2012-0806, printed p.8, equations 1/2/4 and Table VIII. Original experimental comparison is not validation of this implementation. Accessed 2026-10-24. <https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf>
- **NASA-GFSSP6-NONCIRCULAR**: GFSSP v6 manual, printed p.40/PDF p.80, equation 28 and cross-section-dependent Poiseuille convention. Accessed 2026-10-24. <https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf>
- **EPANET22-FRICTION-SOURCE**: EPANET tag v2.2, hydcoeffs.c, frictionFactor and base-ten conversion; primary software convention, not cooling-channel validation. Accessed 2026-10-24. <https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c>
- **SWAMEE-JAIN-1976-IDENTITY**: Swamee and Jain, Explicit Equations for Pipe-Flow Problems, J. Hydraulics Division 102(5), 657-664; publisher abstract consulted, full text not recovered. Accessed 2026-10-24. <https://doi.org/10.1061/JYCEAJ.0004542>
- **AUTODESK-LWPOLYLINE**: Autodesk, lightweight-polyline subclass, vertex/count and closure group codes. Accessed 2026-10-24. <https://help.autodesk.com/cloudhelp/2018/ENU/AutoCAD-DXF/files/GUID-748FC305-F3F2-4F74-825A-61F04D757A50.htm>
- **AUTODESK-INSUNITS**: Autodesk, drawing insertion unit convention (6 = metres). Accessed 2026-10-24. <https://help.autodesk.com/cloudhelp/2025/ENU/AutoCAD-Core/files/GUID-A58A87BB-482B-4042-A00A-EEF55A2B4FD8.htm>
- **LOC-STL-FAMILY**: Library of Congress, oriented triangular surface and format limitations; original specification not recovered by its compilers. Accessed 2026-10-24. <https://www.loc.gov/preservation/digital/formats/fdd/fdd000504.shtml>
- **LOC-STL-ASCII**: Library of Congress, ASCII STL records and metadata limitations. Accessed 2026-10-24. <https://www.loc.gov/preservation/digital/formats/fdd/fdd000506.shtml>
- **W3C-XML10**: W3C, XML 1.0 fifth edition, permitted characters and literal character data. Accessed 2026-10-24. <https://www.w3.org/TR/xml/>
- **NASA-MOC-CARTESIAN-1993**: Aldo and Argrow, Cartesian characteristic and compatibility equations (1a,b), PDF pp.3-4; printed sigma convention inconsistency retained and checked against cylindrical continuity. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19940019183>
- **NASA-TRR33-AXISYMMETRIC**: Guentert and Neumann, 1959, radius source and averaged characteristic coefficients, printed pp.3-4/PDF pp.5-6; full variable-composition results not reproduced here. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19980232191>
- **NASA-TMX1502-PLANAR**: Vanco and Goldman, 1968, planar sharp-throat nozzle, finite-grid sensitivity and printed Table II. Screened future full-nozzle reference; not validation of the current mapped field. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19680005278>
- **NASA-TND4960-AXIS-CONTINUITY**: Anderson, 1969, Appendix C equation C1 and streamline definitions, printed p.30/PDF p.33; cylindrical continuity source only. Accessed 2026-10-24. <https://ntrs.nasa.gov/citations/19690005234>

## Foundational textbooks and reports (bibliographic details verified 2026-10-25)

- **NACA-TR-1135**: Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*, NACA Report 1135, 1953 — one-dimensional isentropic flow, normal and oblique shocks, Prandtl–Meyer expansion. <https://ntrs.nasa.gov/citations/19930091059>
- **NACA-TN-1336**: B. L. Hicks, D. J. Montgomery and R. H. Wasserman, *On the One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with Friction and Heat Addition*, NACA TN 1336, 1947 — constant-area duct flow with friction (Fanno) and heat addition (Rayleigh). <https://ntrs.nasa.gov/citations/19930081963>
- **RAO-1958**: G. V. R. Rao, *Exhaust Nozzle Contour for Optimum Thrust*, Jet Propulsion 28(6), 377–382, 1958. doi:10.2514/8.7324
- **SUTTON-BIBLARZ-2017**: G. P. Sutton and O. Biblarz, *Rocket Propulsion Elements*, 9th ed., Wiley, 2017, ISBN 978-1-118-75365-1 — ideal rocket, nozzle, thrust-coefficient, staging and feed-system relations (bibliographic citation).
- **MCKAY-1979**: M. D. McKay, R. J. Beckman and W. J. Conover, *A Comparison of Three Methods for Selecting Values of Input Variables in the Analysis of Output from a Computer Code*, Technometrics 21(2), 239–245, 1979 — Latin hypercube sampling. doi:10.2307/1268522
- **PAULE-MANDEL-1982**: R. C. Paule and J. Mandel, *Consensus Values and Weighting Factors*, J. Res. Natl. Bur. Stand. 87(5), 377–385, 1982 — Mandel–Paule between-group variance. doi:10.6028/jres.087.022

Not cited because the bibliographic record could not be verified in this
repository's review: D. R. Bartz (1957) original article (use **NASA-SP125-BARTZ-PAGES**)
and I. M. Hall (1962) journal article (use **HALL-MCCABE-RM3440** and **KLIEGEL-QUAN-1966**).

## NASA CEA data and program (verified 2026-10-25)

- **NASA-CEA-THERMO**: NASA Glenn Research Center, CEA thermodynamic database `data/thermo.inp`, repository commit 2f79a647e85737742fb08a127da00d1ea1d069da (2026-09-25), SHA-256 fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363, Apache-2.0. Source of the 26 added gas species and the six condensed reactants (H2(L), O2(L), CH4(L), RP-1, N2O4(L), CH6N2(L)) in `thermochemistry/data/nasa_cea_thermo.json`; fixed-column layout per RP-1311 Part II Appendix A. The parser reproduces the ten bundled TP-2002-211556 species exactly except the H-atom b1 constant (25473.70801 here vs 25473.6599). <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp>
- **NASA-CEA-PY**: NASA Glenn Research Center, *CEA: Chemical Equilibrium with Applications*, Python package 3.3.4 (released 2026-08-28), Apache-2.0, packaged `thermo.lib` SHA-256 8e5df1cca92d4a48663d1ee5a1372e6508c59cddc2247ceeca32f041a03ec52a. Independent implementation used for the nine-point release gate (`tests/reference_data/cea_gate_points.json`). <https://pypi.org/project/cea/>
- **NASA-CEA-EXAMPLE8-OUT**: NASA CEA regression output for RP-1311 example 8, `test/main_interface/reference_output/example8.out`, same commit; confirms that the recorded RP-1311 exit values belong to the area-ratio-75 column. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/test/main_interface/reference_output/example8.out>
- **KSC-2010-233**: L. G. Delgado, *Expendable Launch Vehicles Briefing and Basic Rocketry Physics*, NASA Kennedy Space Center, KSC-2010-233, 2010, NTRS 20100036816, PDF p.46 "Examples of Specific Impulse": LH2/LOX 455 s, kerosene (RP-1)/LOX 358 s, hydrazine/N2O4 305 s, solid motor 250 s; vacuum or sea-level condition not stated. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20100036816>

## Thermal hardware and turbopump sources (verified 2026-10-25)

Each entry was opened and read on 2026-10-25; the locator names the table,
page or query that states every value used by
`thermodynamics/materials.py`, `thermodynamics/wall_conduction.py` and
`propulsion/turbopump.py`. Catalog values are reference-condition data, not
validation of this library.

- **NIST-WEBBOOK-FLUID**: E. W. Lemmon, M. O. McLinden and D. G. Friend, *Thermophysical Properties of Fluid Systems*, NIST Chemistry WebBook SRD 69. `fluid.cgi` isobaric queries at 1.0000 MPa (`Digits=5`): parahydrogen ID B5000001 row 20.000 K, methane ID C74828 row 111.00 K, water ID C7732185 row 300.00 K; critical points from the final row of each full `Type=SatT` saturation table. Accessed 2026-10-25. <https://webbook.nist.gov/chemistry/fluid/>
- **NISTIR-6646-RP1**: J. W. Magee et al., *Thermophysical Properties Measurements and Models for Rocket Propellant RP-1: Phase I*, NISTIR 6646, 2007. Table 8 (printed p. 35, PDF p. 41) density row 23.283 °C/83.72 kPa; Table 10 (printed p. 48, PDF p. 54) point 2056 radiation-corrected conductivity; Table 11 (printed p. 61, PDF p. 67) kinematic viscosity at 298.15 K. No heat-capacity table is printed (section 5). Accessed 2026-10-25. <https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf>
- **AURUBIS-C10200**: Aurubis, *Technical Datasheet C10200 (Cu-OF)*, revision 18 08 US, p. 1 physical-properties table (melting point, density, thermal conductivity at 20 °C). Manufacturer typical values. Accessed 2026-10-25. <https://www.aurubis.be/en/dam/jcr:90029b00-b968-486e-bde0-ce56f16022cd/c10200-cu-of-us.pdf>
- **NASA-TM-2005-213566**: D. L. Ellis, *GRCop-84: A High-Temperature Copper Alloy for High-Heat-Flux Applications*, NASA/TM-2005-213566, printed p. 8 (PDF p. 12): thermal conductivity between 305 and 320 W/m·K over the MCC-liner operating range. Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf>
- **NASA-CR-2000-210055**: D. L. Ellis and D. J. Keller, *Thermophysical Properties of GRCop-84*, NASA/CR-2000-210055, printed p. 13 (PDF p. 17): with a 100 MPa yield-strength criterion GRCop-84 can be used up to 973 K (NARloy-Z 773 K). Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf>
- **SMC-045-INCONEL-718**: Special Metals Corporation, *INCONEL alloy 718*, Publication SMC-045, 2007: p. 1 use range −423 to 1300 °F and Table 2 (density, melting range); Table 5 p. 2 (thermal conductivity, annealed and aged). Typical values. Accessed 2026-10-25. <https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf>
- **AKSTEEL-304-PDB**: AK Steel, *Type 304/304L Stainless Steel* product data bulletin, 2018, PDF p. 4: density, thermal conductivity at 100 °C, melting range 1399–1454 °C, continuous scaling limit 899 °C. Accessed 2026-10-25. <https://Protolabs.com/media/1020110/304-304l-stainless.pdf>
- **NIST-SP811-B8**: NIST SP 811, Appendix B.8, factors listed alphabetically: °C and °F to K, g/cm³ to kg/m³. Accessed 2026-10-25. <https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8>
- **NIST-SP811-B9**: NIST SP 811, Appendix B.9, factors by kind of quantity: lb/in³ (2.767 990 E+04 kg/m³), Btu_IT·in/(h·ft²·°F) (1.442 279 E−01 W/(m·K)), U.S. gpm (6.309 020 E−05 m³/s), ft (3.048 E−01 m). Accessed 2026-10-25. <https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9>
- **NASA-SP125-REGEN-WALL**: D. K. Huzel and D. H. Huang, *Design of Liquid Propellant Rocket Engines*, NASA SP-125, printed p. 104 (PDF p. 113), equations (4-19)–(4-22): `q = (k/t)(Twg − Twc) = hc(Twc − Tco)`; original page image inspected. Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- **NASA-SP8109-PUMP**: NASA, *Liquid Rocket Engine Centrifugal Flow Turbopumps*, NASA SP-8109, December 1973: equations (1)–(2) printed p. 6 (PDF p. 18), specific-speed experience 450–2100 printed p. 7, equation (8) and 2 %/10 % head-loss basis printed p. 11 (PDF p. 23), H-Q slope recommendation printed p. 67 (PDF p. 79), affinity caution printed p. 83, NPSH glossary printed p. 88 (PDF p. 100); page images inspected. Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf>
- **NASA-GFSSP6-AREA-CHANGE**: NASA GFSSP v6 user manual, section 3.1.7.5 equation (33) `Kf = K1/(2 g_c rho_u A^2)` with the upstream area, printed p. 42 (PDF p. 82); section 3.1.7.7 square reduction, equations (38)-(39), printed p. 44 (PDF p. 84); section 3.1.7.8 square expansion, equations (40)-(41), printed p. 45 (PDF p. 85); the manual's reference 18 is the source of these K1 forms. Accessed 2026-10-26. <https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf>
- **NASA-GLENN-CD-NOZZLE**: NASA Glenn Research Center, *Nozzle Design - Converging/Diverging (CD) Nozzle*: momentum `r V dV = -dp` with continuity `mdot = r V A = constant` in the area-velocity derivation. Accessed 2026-10-26. <https://www.grc.nasa.gov/www/k-12/airplane/nozzled.html>
- **NASA-GFSSP6-PUMP**: NASA GFSSP v6 user manual, section 3.1.7.14 equation (53) `Δp = A0 + B0 ṁ + C0 ṁ²`, printed p. 51 (PDF p. 91); section 5.4.4 turbopump option, Head/Speed² and Torque/(Density·Speed²) versus flow rate/speed, printed p. 110 (PDF p. 150). Accessed 2026-10-25. <https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf>
- **DOE-PUMP-SOURCEBOOK**: U.S. DOE and Hydraulic Institute, *Improving Pumping System Performance: A Sourcebook for Industry*, 2nd ed., 2006: NPSHA/NPSHR and 25 % margin design rule printed p. 38 (PDF p. 42); specific speed in rpm, gpm and ft printed p. 39; impeller-trimming relations printed p. 49 and 70 % practice printed p. 50 (PDF pp. 53–54). Accessed 2026-10-25. <https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf>

Not catalogued because no primary value was verified in this review: RP-1
isobaric heat capacity and directly published dynamic viscosity; GRCop-84
density and solidus; NARloy-Z conductivity, density and melting range; Type
347 stainless steel properties.

## Calibration field sources (verified 2026-10-25)

- **JCGM-100-2008-SEC5**: JCGM 100:2008, *Evaluation of measurement data — Guide to the expression of uncertainty in measurement*, printed pp.19–22 (PDF pp.31–34) inspected 2026-10-25: §5.1.2 equation (10), law of propagation for uncorrelated inputs; §5.1.6 equation (12) and NOTE 2, relative combined variance for products and quotients; §5.2.2 equations (13)–(16), correlated inputs with correlation coefficient −1 ≤ r ≤ +1; NOTE 1, for r = +1 the combined standard uncertainty is the linear sum of the c_i u(x_i); NOTE 4, correlation terms for the power-law form. Used by `calibration_field` for the linear-envelope consensus uncertainty and the thrust/Isp flow-scale bound. <https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf>
- **NIST-TN1297-COVARIANCE** (re-opened 2026-10-25): NIST TN 1297 Appendix A, equation (A-3) with the covariance term u(x_i, x_j); consistent with JCGM 100:2008 equation (13). <https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-law-propagation-uncertainty>

## Feed dynamics sources (verified 2026-10-25)

Read during the 2026-10-25 implementation of `water_hammer.py`,
`regulator_hysteresis.py` and `ullage_stratification.py`; locators are the
pages actually inspected.

- **NASA-TM-20210022275-TRANSIENTS**: K. Imtiaz, V. Goyal, P. Babuska, E. Barbour and J. Smith, *Treatment of Transient Pressure Events in Space Flight Pressurized Systems*, NASA/TM-20210022275 (NESC-NPP-21-01679), October 2021. Printed p. 13 (§6.1, priming described qualitatively); p. 19, eq. 1 (Joukowsky `dP = -rho a dV`) and eq. 2 (wave speed with bulk modulus, pipe modulus, D, t and boundary factor C; thick-walled forms for D/t <= 25); p. 20, Figures 4-5 (boundary factor C, cases I-IV); p. 22, Case Study 8-1 (0.088 lb/s, 0.62 in², a = 39,370 in/s, rho = 0.036 lb/in³, Joukowsky 14.5 psi recovered by the method of characteristics); pp. 24-25 (gradual closure lower than instantaneous); p. 54, eqs. A-1 to A-6; pp. 55-57, Appendix B, eqs. B-1 to B-17 (continuity and momentum equations, characteristic compatibility equations, B = rho a / A, R = rho f dx / (2 A² D), interior-point solution). <https://ntrs.nasa.gov/citations/20210022275>
- **NASA-20190030443-PRIMING**: A. LeClair, A. Boehm and A. Majumdar (NASA/MSFC), *Fluid Transient Analysis of Propellant Feedlines during a Priming Event*, TFAWS 2019 presentation. Slides 2-3 (priming of an evacuated line; dead-end surge "as high as" rho a dV), slide 17 (linear and power-law valve-opening profiles), slides 5-21 (comparisons with Lee 2005, Prickett et al. 1992 and Moore et al. 2019 data). No closed-form priming model is given; none is implemented. <https://ntrs.nasa.gov/citations/20190030443>
- **NASA-SP8080-FORCE-BALANCE**: L. E. Tomlinson (author), NASA SP-8080, *Liquid Rocket Pressure Regulators, Relief Valves, Check Valves, Burst Disks, and Explosive Valves*, March 1973. Printed p. 9 (§2.1, lockup pressure set a minimum above the maximum regulated pressure); p. 11 (§2.1.1.1, supply-pressure and flow-demand error limited to about 25 % of the regulated-pressure variation; mechanical hysteresis and friction as error sources); pp. 72-73 (§3.1.1.1, eq. 2 direct-acting force balance, eq. 3 regulated-pressure error); p. 76 (§3.1.2.1, orifice discharge coefficient 0.6 unless better known; §3.1.2.3, hysteresis and friction error at most 25 % of tolerance); p. 104 (Glossary: lockup, Coulomb damping); p. 105 (response: 63 % of signalled position). PDF pp. 21, 23, 84-88, 116-117 inspected. <https://ntrs.nasa.gov/citations/19740002611>
- **NASA-SP8112-STRATIFICATION**: NASA SP-8112, *Pressurization Systems for Liquid Rockets*, 1975. Printed p. 52, eq. 13 and ullage fractions; p. 53, eq. 14 (`M_u = V_t (rho_f)_u`, density at the lesser of minimum propellant or tank-surface temperature and maximum ullage pressure); p. 90, §2.3.1.2 (stratified ullage from inlet temperature at the top to propellant temperature at the interface; stratification reduces heat transfer to the propellant and the gas mass). PDF pp. 64-65 and 102 inspected. <https://ntrs.nasa.gov/citations/19760015212>
- **NASA-19980236009-PRESSURIZATION**: A. Majumdar and T. Steadman (Sverdrup Technology), *Numerical Modeling of Pressurization of a Propellant Tank*, draft for the 37th AIAA Aerospace Sciences Meeting (AIAA 99-0879). Nomenclature (g_c = 32.174 lbm ft/(lbf s²)); p. 3, eq. 2 (branch friction `K_f mdot abs(mdot)`); p. 4, eqs. 9-13 (`h_c = K_H C (k/L_s) (Gr Pr)**n`); p. 5 (C = 0.27, n = 0.25, K_H = 1 ullage-to-propellant; C = 0.54, n = 0.25 ullage-to-wall, both "according to Ring"; tank diameter as length scale; eqs. 15-17 wall). The exponent on the length in the scanned eq. 12 is not legible. <https://ntrs.nasa.gov/citations/19980236009>
- **NASA-20190030427-GRASHOF**: A. Bandyopadhyay, A. K. Majumdar, A. C. LeClair and J. G. Valenzuela, *Multi-Node Modeling of Cryogenic Tank Pressurization System using Generalized Fluid System Simulation Program*. PDF p. 1 nomenclature (`Gr = g beta dT L**3 rho**2 / mu**2`, `Ra = Gr Pr`); p. 6 (`Nu = h L / k`; the inline Grashof there is printed with `L2`, inconsistent with p. 1 and dimensionally wrong); p. 7 (K_H set to 0.5 to match a measured boil-off rate). <https://ntrs.nasa.gov/citations/20190030427>
- **NASA-20240006787-MULTINODE-ULLAGE**: A. LeClair and A. Majumdar (NASA/MSFC), *Nodal Modeling of Tank Pressurization and Draining using a Multi-Node-Ullage Approach*, extended abstract, AIAA SciTech 2025. Pp. 2-3 (ullage in six layers of equal height, wall solid nodes, re-discretized each step); p. 4, eqs. 1-2 (`Nu = 0.54 Ra**0.25` ullage-wall, `Nu = 0.27 Ra**0.25` ullage-propellant); p. 6 (natural-convection correlation under-predicts measured helium consumption). <https://ntrs.nasa.gov/citations/20240006787>
- **NASA-TND2797-STRATIFIED-OUTFLOW**: D. A. Mandell and W. H. Roudebush, *Parametric Investigation of Liquid Hydrogen Tank Pressurization During Outflow*, NASA TN D-2797, April 1965. Printed pp. 2-3, assumptions (1)-(13) and eqs. 1-3 (one-dimensional gas and wall energy equations, no axial conduction, wall at the liquid temperature at the interface); p. 8, eqs. 23-25 (pressurant mass, ideal mass `m_i = pi r² (L_f - L_0) rho_g`, mass ratio). <https://ntrs.nasa.gov/citations/19650012768>
- **NASA-GRC-MASS-FLOW-EQ10**: NASA Glenn Research Center, *Mass Flow Choking*: eq. 10 (compressible mass flow in terms of total pressure, total temperature, area and Mach number) and the choked (M = 1) form. Accessed 2026-10-25. <https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html>
- **NASA-GRC-ISENTROPIC**: NASA Glenn Research Center, *Isentropic Flow*: eqs. 6-7 (`p/pt` and `T/Tt` against Mach number). Accessed 2026-10-25. <https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html>

Read but not used for equations: NESC Technical Bulletin 22-03 (names the
Joukowsky equation without printing it) and Holt, Majumdar, Steadman and
Hedayat, AIAA-2000-3719 (<https://ntrs.nasa.gov/citations/20000090555>, a
±3 psi pressure-switch control band in a GFSSP model). Not read and not cited:
Wylie and Streeter, *Fluid Transients* (the MOC here follows the NASA/TM
Appendix B, which presents the same characteristic method); E. Ring, *Rocket
Propellant and Pressurization Systems* (1964), known here only through the
NASA MSFC papers that quote its constants; the GFSSP 6 user manual
pressurization chapter (download exceeded the fetch limit).

## Real-fluid equations of state (verified 2026-10-25)

- **LEACHMAN-2009-H2**: J. W. Leachman, R. T Jacobsen, S. G. Penoncello and E. W. Lemmon, *Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen*, J. Phys. Chem. Ref. Data 38(3), 721-748, 2009, doi:10.1063/1.3160306 — R on p. 724; Table 3 p. 728; Table 4 p. 729; Tables 5-6 and Eqs. (31)-(32) p. 730; Table 7 p. 731; Table 8 p. 732; Eq. (33) p. 733; validity (13.8033/13.957-1000 K, 2000 MPa) in the abstract. Read from the NIST-hosted publisher reprint. Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.3160306.pdf>
- **STEWART-1991-O2**: R. B. Stewart, R. T Jacobsen and W. Wagner, *Thermodynamic Properties of Oxygen from the Triple Point to 300 K with Pressures to 80 MPa*, J. Phys. Chem. Ref. Data 20(5), 917-1021, 1991, doi:10.1063/1.555897 — R, M, fixed points and ideal-gas reference state p. 919; Eq. (7) form p. 924; Eq. (7) and (8) coefficients p. 925; Eq. (10) and Eq. (11) p. 926; Eqs. (13)-(14) p. 927; Table 4 and Eq. (15) p. 928 (the Schmidt-Wagner equation as printed here). Scanned reprint; values compared between OCR text and rendered page images. Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.555897.pdf>
- **SCHMIDT-WAGNER-1985-O2**: R. Schmidt and W. Wagner, *A New Form of the Equation of State for Pure Substances and its Application to Oxygen*, Fluid Phase Equilib. 19(3), 175-200, 1985 — original source of the oxygen equation; not opened (publisher paywall); coefficients taken from STEWART-1991-O2 Table 4, which states that it presents this formulation. Citation as displayed on the NIST WebBook oxygen page. <https://doi.org/10.1016/0378-3812(85)87016-3>
- **SPAN-2000-N2**: R. Span, E. W. Lemmon, R. T Jacobsen, W. Wagner and A. Yokozeki, *A Reference Equation of State for the Thermodynamic Properties of Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to 2200 MPa*, J. Phys. Chem. Ref. Data 29(6), 1361-1433, 2000, doi:10.1063/1.1349047 — physical constants p. 1364; Maxwell criterion p. 1365; ancillary Eqs. (4)-(6) p. 1368; Eq. (53) p. 1382; Eq. (55) p. 1384; Tables 17-18 and property relations Eqs. (56)-(65) p. 1385; Eqs. (68)-(86) p. 1386. Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.1349047.pdf>
- **SETZMANN-WAGNER-1991-CH4**: U. Setzmann and W. Wagner, *A New Equation of State and Tables of Thermodynamic Properties for Methane Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa*, J. Phys. Chem. Ref. Data 20(6), 1061-1155, 1991, doi:10.1063/1.555898 — physical constants p. 1064 (prints p_c = 4.5922 MPa, a misprint of the 4.5992 MPa used in Eq. (2.3), Table 26 and p. 1091); Eqs. (3.2), (3.4) p. 1082; Table 26 p. 1086; Eqs. (5.1)-(5.3), Table 33 and validity p. 1091; Tables 34-35 p. 1092. Scanned reprint; values compared between OCR text and rendered page images. Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.555898.pdf>
- **ORTIZVEGA-2023-HE**: D. O. Ortiz Vega, K. R. Hall, J. C. Holste, A. H. Harvey and E. W. Lemmon, *An Equation of State for the Thermodynamic Properties of Helium*, NIST IR 8474, National Institute of Standards and Technology, August 2023 — Table 1 p. 2; Eqs. (1)-(4) pp. 2-3; Table 2 p. 4; range of validity p. 15; check values Tables 3-4 p. 16. Accessed 2026-10-25. <https://doi.org/10.6028/NIST.IR.8474>
- **BELL-JAGER-2016-CUBIC**: I. H. Bell and A. Jäger, *Helmholtz Energy Transformations of Common Cubic Equations of State for Use with Pure Fluids and Mixtures*, J. Res. NIST 121, 238-263, 2016 — Eq. (9) p. 241; Table 1 and Eq. (13) p. 242; Sec. 2.2.1 p. 243; Eqs. (25)-(29) p. 244; reference list (Soave 1972, Peng-Robinson 1976, Robinson-Peng 1978) pp. 262-263. Accessed 2026-10-25. <https://doi.org/10.6028/jres.121.011>
- **PENG-ROBINSON-1976**: D.-Y. Peng and D. B. Robinson, *A New Two-Constant Equation of State*, Ind. Eng. Chem. Fundam. 15(1), 59-64, 1976 — not opened (publisher paywall); constants used as printed in BELL-JAGER-2016-CUBIC Table 1; citation as listed there (reference 17) and on the NIST TRC TDE help page (accessed 2026-10-25, https://trc.nist.gov/TDE/TDE_Help_4_0/eos-Peng-Robinson.htm). <https://doi.org/10.1021/i160057a011>
- **SOAVE-1972**: G. Soave, *Equilibrium Constants from a Modified Redlich-Kwong Equation of State*, Chem. Eng. Sci. 27(6), 1197-1203, 1972 — not opened (publisher paywall); constants used as printed in BELL-JAGER-2016-CUBIC Table 1; citation as listed there (reference 15). <https://doi.org/10.1016/0009-2509(72)80096-4>
- **ROBINSON-PENG-1978**: D. B. Robinson and D.-Y. Peng, *The Characterization of the Heptanes and Heavier Fractions for the GPA Peng-Robinson Programs*, Gas Processors Association Research Report RR-28, 1978 — not opened; the omega > 0.491 form of m(omega) is used as printed in BELL-JAGER-2016-CUBIC Table 1 (reference 23). <https://doi.org/10.6028/jres.121.011>
- **COOLPROP-FLUIDS-C00300E**: CoolProp developers, fluid definition files `dev/fluids/{ParaHydrogen,Hydrogen,Oxygen,Nitrogen,Methane,Helium}.json` and `CoolPropBibTeXLibrary.bib`, commit c00300e97ad4cdb30032635d7379209d5e45a121 — used only as an independent machine-readable transcription to cross-check every residual coefficient (all matched the papers); the files cite the papers above (BibTeX keys Leachman-JPCRD-2009, Schmidt-FPE-1985, Stewart-JPCRD-1991, Span-JPCRD-2000, Setzmann-JPCRD-1991, OrtizVega-NISTIR-2023). Accessed 2026-10-25. <https://github.com/CoolProp/CoolProp/tree/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids>
- **NIST-WEBBOOK-FLUID-RF-2026**: E. W. Lemmon, M. O. McLinden and D. G. Friend, *Thermophysical Properties of Fluid Systems*, NIST Chemistry WebBook SRD 69 — `fluid.cgi` Action=Data, Digits=12 tables for parahydrogen (B5000001), normal hydrogen (C1333740), oxygen (C7782447), nitrogen (C7727379), methane (C74828) and helium (C7440597): 106 single-phase rows (Type=IsoBar) and 36 saturation pairs (Type=IsoTherm), each URL stored in `tests/reference_data/nist_real_fluids.json`; equation-of-state citations, critical constants and acentric factors from the fluid pages. Accessed 2026-10-25. <https://webbook.nist.gov/chemistry/fluid/> <https://webbook.nist.gov/cgi/fluid.cgi>
- **NIST-WEBBOOK-H2-SPECIES**: NIST Chemistry WebBook SRD 69, *Hydrogen* species record — molecular weight 2.01588 g/mol (not printed in LEACHMAN-2009-H2). Accessed 2026-10-25. <https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI>

## Transport properties and equilibrium derivatives (verified 2026-10-25)

- **NASA-CEA-TRANS**: NASA Glenn Research Center, CEA transport-property database `data/trans.inp`, repository commit 2f79a647e85737742fb08a127da00d1ea1d069da, Apache-2.0, file SHA-256 379c514a7f7638371d8d9254a1d653cd078084a9371061aa0b45c5a5c3334d41; 66 species records, 41 binary interaction records, 379 coefficient lines. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/trans.inp>
- **NASA-CEA-SRC-TRANSPORT**: NASA Glenn Research Center, CEA source `source/equilibrium.f90`, subroutine `compute_transport_properties` (lines 5097-5763 at commit 2f79a647e85737742fb08a127da00d1ea1d069da): mixture rules, rigid-sphere estimate `5.656854d0` (lines 5576, 5583), `psi` factor 2.41/0.142 (5614-5617), `rtpd` with 1.1 (5670), reaction cp (5718), reaction conductivity factor 0.6 (5733), missing-conductivity estimate 0.00375/0.00132 (5537), transport species selection (5206 onward). Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/equilibrium.f90>
- **NASA-CEA-SRC-TRANSPORT-FITS**: NASA Glenn Research Center, CEA source `source/transport.f90` (interval selection `st_calc_eta`, `st_calc_lambda`, `bt_calc_eta`, lines 83-150) and `source/fits.f90` line 175 (`A*log(T) + B/T + C/(T*T) + D`), commit 2f79a647e85737742fb08a127da00d1ea1d069da. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/transport.f90>
- **NASA-CEA-SRC-TRANS-READER**: NASA Glenn Research Center, CEA source `source/database_compile.f90`, subroutine `utran`, formats 99001 `(2a16, 2x, a1, i1, a1, i1)` and 99002 `(1x, a1, 2f9.2, 4e15.8)` (lines 329-384), commit 2f79a647e85737742fb08a127da00d1ea1d069da. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/database_compile.f90>
- **NASA-CEA-REPO**: NASA Glenn Research Center, CEA repository (Apache-2.0). Accessed 2026-10-25. <https://github.com/nasa/cea>
- **NASA-CEA-PY-TRANSPORT**: NASA Glenn Research Center, `cea` Python package 3.3.4, `RocketSolver`/`EqSolver` with `transport=True`; packaged `trans.lib` SHA-256 c203d9330746405e1b57008cd697354745ab991694b4ca11ed556b18cb834111, `thermo.lib` SHA-256 8e5df1cca92d4a48663d1ee5a1372e6508c59cddc2247ceeca32f041a03ec52a, `cea.R` = 8314.51 J/(kmol K); outputs in millipoise, mW/(cm K), kJ/(kg K). Used only to generate `tests/reference_data/cea_transport_points.json`. Accessed 2026-10-25. <https://pypi.org/project/cea/>
- **NASA-RP1311-P1-TRANSPORT**: Gordon and McBride, NASA RP-1311 Part I (1994), chapter 5 "Thermal Transport Property Data", printed pp. 21-23 (PDF pp. 27-29), equations (5.1)-(5.19b); A* = 1.1 substitution stated under equation (5.13); fit units micropoise and microwatt/(cm K) in section 5.1. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-RP1311-P1-DERIVATIVES**: Gordon and McBride, NASA RP-1311 Part I (1994), sections 2.5 and 2.6, printed pp. 10-12 (PDF pp. 16-18), equations (2.49a)-(2.75). Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19950013764/downloads/19950013764.pdf>
- **USSA1976-SUTHERLAND**: NOAA/NASA/USAF, U.S. Standard Atmosphere 1976, section 1.3.11, printed p. 19, equation (51) with beta = 1.458e-6 kg/(s m K^0.5) and S = 110.4 K; Table III printed p. 101, Z = 0 m: mu = 1.7894e-5 N s/m^2; NTRS errata sheet appended to the PDF (S = 110.4 K; Table III conductivity column off by 1000). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19770009539>
- **IDAES-CHAPMAN-ENSKOG**: IDAES-PSE documentation, "Chapman Enskog Viscosity Model": mu = 26.69 sqrt(M T)/(sigma^2 Omega) micropoise (M g/mol, sigma angstrom), Neufeld form Omega = A T*^-B + C exp(-D T*) + E exp(-F T*), A-F = 1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787; cites Poling et al., The Properties of Gases and Liquids, 5th ed., eqs. 9-3.9 and 9-4.3, and Neufeld, Janzen and Aziz, J. Chem. Phys. 57, 1100 (1972) (the primary paper was not opened). Accessed 2026-10-25. <https://idaes-pse.readthedocs.io/en/main/explanations/components/property_package/general/transport_properties/chapman_enskog.html>
- **IDAES-CHAPMAN-ENSKOG-SRC**: IDAES-PSE source `idaes/models/properties/modular_properties/pure/ChapmanEnskog.py`, `collision_integral_neufeld_callback` constants and the 26.69 uP A^2/sqrt(g/mol K) prefactor; no T* validity range stated. Accessed 2026-10-25. <https://idaes-pse.readthedocs.io/en/main/_modules/idaes/models/properties/modular_properties/pure/ChapmanEnskog.html>
- **SVEHLA-TR-R132**: R. A. Svehla, "Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures", NASA TR R-132, 1962: Appendix B equation (1), printed p. 8; Table I(a) printed p. 23, N2 sigma = 3.798 A, eps/k = 71.4 K; Table III printed p. 90, N2 eta x 1e6 = 72.4, 177.7, 397.1, 624.0 poise at 100, 300, 1000, 2000 K. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19630012982>
- **NASA-TP-2002-211556-SPECIES**: NASA Glenn, NASA/TP-2002-211556, "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species" (not re-opened in this work package; title as recorded in the existing species module), used here for species cp and enthalpy in equations (5.15)-(5.16) and (2.59). <https://ntrs.nasa.gov/citations/20020085330>

## Advanced compressible flow (verified 2026-10-25)

- **ANDERSON-FOA-2017**: J. D. Anderson Jr., *Fundamentals of Aerodynamics*, 6th ed., McGraw-Hill Education, 2017, ISBN 978-1-259-12991-9 — Ch. 9 (shock-expansion theory, Examples 9.6, 9.12, 9.13), Ch. 11 (Prandtl-Glauert, Karman-Tsien, Laitone, critical Mach number, Examples 11.1-11.3), Ch. 12 (linearized supersonic flow, Examples 12.1-12.3), Sec. 13.6 (Taylor-Maccoll conical flow), Ch. 14 (Newtonian theory, hypersonic shock relations, Mach-number independence, Example 14.1) (bibliographic citation).
- **NASA-SP-3004**: J. L. Sims, *Tables for Supersonic Flow Around Right Circular Cones at Zero Angle of Attack*, NASA SP-3004, 1964 — Tables 7, 9 and 11 (surface Mach number, surface pressure coefficient, shock-wave angle), gamma = 1.4. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19640009035>
- **NACA-TR-1194**: M. D. Van Dyke, *A Study of Hypersonic Small-Disturbance Theory*, NACA Report 1194, 1954 — eq. (1) validity criterion, eq. (19) and p. 6 unified similarity rule, eqs. (33)-(35) thin-wedge solution (Linnell). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930092209>

## Incompressible and lifting aerodynamics (verified 2026-10-25)

- **NACA-TR-824**: I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr., "Summary of Airfoil Data", NACA Report 824, 1945 — pp. 4-5 (Eqs. 1-4, numbering systems, mean-line scaling), p. 14 (Eq. 16 and theoretical examples), p. 30 (zero-lift and lift-slope trends), pp. 100-101 (ordinate tables), p. 136 (NACA 2412 data). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930090976>
- **NACA-TR-460**: E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The Characteristics of 78 Related Airfoil Sections from Tests in the Variable-Density Wind Tunnel", NACA Report 460, 1933 — pp. 4-5 (four-digit thickness form, leading-edge radius, mean lines, surface construction). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091108>
- **NACA-TR-537**: E. N. Jacobs, R. M. Pinkerton, "Tests in the Variable-Density Wind Tunnel of Related Airfoils Having the Maximum Camber Unusually Far Forward", NACA Report 537, 1936 — pp. 521-522 (five-digit mean-line equations, Table I). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091610>
- **NACA-TR-610**: E. N. Jacobs, R. M. Pinkerton, H. Greenberg, "Tests of Related Forward-Camber Airfoils in the Variable-Density Wind Tunnel", NACA Report 610, 1937 — p. 2 (five-digit designation; 230 mean line camber 1.8 percent at 0.15c). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091685>
- **NACA-TR-411**: T. Theodorsen, "Theory of Wing Sections of Arbitrary Shape", NACA Report 411, 1931 — pp. 229-233 (Joukowski transformation Eq. I, coordinates Eq. II, circle flow Eqs. V-VII, velocity transformation p. 232, exact surface velocity Eq. XII). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091485>
- **NACA-TR-572**: R. F. Anderson, "Determination of the Characteristics of Tapered Wings", NACA Report 572, 1937 — pp. 3-5 (Glauert series for the section lift coefficient, Eq. 5 with general section slope m0). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091647>
- **DLMF-3.5**: NIST Digital Library of Mathematical Functions, Sec. 3.5(v) Gauss quadrature, Gauss-Legendre formula and Table 3.5.1. Accessed 2026-10-25. <https://dlmf.nist.gov/3.5>

## Viscous flow and aerothermodynamics (verified 2026-10-25)

- **SUTTON-GRAVES-1971-TRR376**: K. Sutton and R. A. Graves Jr., *A General Stagnation-Point Convective-Heating Equation for Arbitrary Gas Mixtures*, NASA TR R-376, November 1971 — Eqs. (33), (41)-(43), Tables I-II. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19720003329>
- **TAUBER-1989-TP2914**: M. E. Tauber, *A Review of High-Speed, Convective, Heat-Transfer Computation Methods*, NASA TP-2914, 1989 — Eqs. (39)-(41), (45)-(49), (52)-(59). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19890017745>
- **HOPKINS-1972-TND6945**: E. J. Hopkins, *Charts for Predicting Turbulent Skin Friction from the Van Driest Method (II)*, NASA TN D-6945, 1972 — Eqs. (1)-(18), Table 2. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19730001588>
- **ALLEN-EGGERS-NACA-R1381**: H. J. Allen and A. J. Eggers Jr., *A Study of the Motion and Aerodynamic Heating of Ballistic Missiles Entering the Earth's Atmosphere at High Supersonic Speeds*, NACA Report 1381 (supersedes NACA TN 4047, 1957) — Eqs. (1)-(51). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091020>
- **CHAPMAN-1958-TN4276**: D. R. Chapman, *An Approximate Analytical Method for Studying Entry into Planetary Atmospheres*, NACA TN 4276, May 1958 — polar-coordinate equations of motion, Eqs. (3)-(7). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930085059>
- **USSA-1976**: NOAA/NASA/USAF, *U.S. Standard Atmosphere, 1976*, NOAA-S/T 76-1562 — Table 2 (r0, g0), Table 4, Eqs. (17)-(19), (23), (33a)-(33b), (42), Table I (geometric-altitude density). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19770009539>

## Launch-vehicle aerodynamics (verified 2026-10-25)

- **NACA-TR-1048**: H. J. Allen and E. W. Perkins, *A Study of Effects of Viscosity on Flow over Slender Inclined Bodies of Revolution*, NACA Report 1048, 1951 (supersedes NACA TN 2044). Printed p. 8: cross force per unit length `f = q dS/dx sin2a cos(a/2) + 2 R c_d q sin^2 a`; p. 9: equations 14 (C_L, C_Dp, C_M) and 15 (small-angle form), symbols A_p, S_b, Q, L, A, X; figure 6 body dimensions (generating radius 25.000 in, nose 4.314 in, diameter 0.750 in, fineness 21.1 and 13.1); figure 7 body-I characteristics at M = 1.98; appendix p. 13: c_dc = 1.2 at low cross-flow Mach number, finite-length ratios 0.755 and 0.692, x(alpha=90) at the planform centroid. Copy read: pdas.com scan. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930090962> <https://www.pdas.com/refs/rep1048.pdf>
- **BARROWMAN-1967**: J. S. Barrowman, *The Practical Calculation of the Aerodynamic Characteristics of Slender Finned Vehicles*, M.S. thesis, Catholic University of America, March 1967 (NTRS record 20010047838). Section 2.40 assumptions; equations 3-1 to 3-12 (subsonic fin C_N_alpha, Diederich correlation), 3-24 to 3-30 (fin CP), 3-57 to 3-65 and 3-77 to 3-89 (slender-body C_N_alpha and CP), 3-90 (tangent-ogive volume, normalisation discrepancy noted), 3-93 to 3-98 (interference; 3-96 printed without the +pi/4 term), 3-106/3-107 (total C_N_alpha and CP); chapter 4: equations 4-2 to 4-16 (skin friction, compressibility, body friction), 4-40 to 4-46 (subsonic forebody pressure drag), 4-48 to 4-56 (base drag; 0.29 in 4-50 vs 0.029 in 4-48/4-49/4-52 and a missing square in 4-53 noted). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20010047838>
- **NACA-TR-1307**: W. C. Pitts, J. N. Nielsen and G. E. Kaattari, *Lift and Center of Pressure of Wing-Body-Tail Combinations at Subsonic, Transonic, and Supersonic Speeds*, NACA Report 1307, 1959. Printed p. 4, equation 14 (slender-body K_W(B), limits 1 and 2) and figure 2; printed p. 6, equation 21 (K_B(W), limits 0 and 2). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930091008>
- **NASA-SP-8035**: NASA, *Wind Loads During Ascent*, NASA Space Vehicle Design Criteria (Structures) SP-8035, June 1970. Section 1 (printed p. 1: maximum wind loads at 6-15 km where high dynamic pressure and severe winds combine; load sources); section 2.2.1.3 (printed p. 9: gust/wind-shear combination, unity correlation, MSFC 0.85 factor); section 4.2.2 (printed p. 17: wind vector combined with vehicle velocity for wind-induced angle of attack; dispersion loads root-sum-squared with buffeting and gust loads and superimposed on quasi-steady loads); Appendix A.2 (printed p. 24: one-minus-cosine gust equations 1-5, 7.65 m/s, thickness 30-275 m, a1 = 0.9215 m, a2 = 0.9137 m; equation 1 printed as 0.765 dH, inconsistent with the stated tangency). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19700027637>
- **SHEA-2020-SLS-ASCENT-DB**: P. R. Shea, J. T. Pinier, H. P. Houlden, A. L. Favaregh, M. J. Hemsch, D. J. Dalle, S. E. Rogers, J. G. Meeroff and H. C. Lee, *Ascent Aerodynamic Force and Moment Database Development for the Space Launch System* (SLS Block 1B; Ames UPWT wind-tunnel data and FUN3D CFD). Nomenclature (six body-axis coefficients); section V.B (natural-neighbour response surfaces per Mach number, interrogated at alpha/beta breakpoints for look-up tables); section VII (range analysis, 99.7% coverage bounds, separate Mach 0.5-1.3 and 1.6-2.5 uncertainties, equations 2-3 net modelling/validation terms, equation 4 final bound with UF = 1.05 and U_WTcorr = U_bal). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20200002382>
- **FAVAREGH-2016-SLS-UQ**: A. L. Favaregh, H. P. Houlden and J. T. Pinier, *Quantification of the Uncertainties for the Space Launch System Liftoff/Transition and Ascent Databases*, AIAA 2016-0795. Equation 1 (C_Q -+ U_C,Q bounds), 2 (RSS of U_EXP, U_DBM, U_bal), 3 (sigma_rep = R/d2), 4-5 (U_EXP = 3 sqrt(n) sigma_rep), 6 (modelling residuals), 7-8 and 11 (net modelling uncertainty), 9, 10, 12 (UF times RSS); ascent uncertainties a function of Mach only. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20160007692>
- **LEE-2024-ARTEMIS-UQ**: M. W. Lee, P. R. Shea, J. T. Pinier and D. T. Chan, *Uncertainty Quantification of Artemis I Space Launch System Integrated Aerodynamics Databases*, Journal of Spacecraft and Rockets (NTRS 20240000054, author manuscript "2023jsr_slsUqSummary_rev1_02"). Equations 1-3 (3-sigma calibration, range method with d2 = 1.128, 1.693, 2.059, 2.326 for n = 2-5, 3-sigma repeatability), 4 (database returns nominal and centred margins), 5-11 (transition, liftoff, full-stack, centerbody, booster-separation and panel-jettison models; uncertainty factors inside or outside the RSS); section II.B (natural-neighbour response surfaces, kriging for booster separation, table-interpolation uncertainty). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20240000054>
- **HEMSCH-2016-SLS-CFD-WT**: M. J. Hemsch, *Global Comparison of CFD and Wind-Tunnel-Derived Force and Moment Databases for the Space Launch System*, AIAA 2016-0794. Abstract (739 RANS solutions, line loads integrated to forces and moments; 0.5 <= M <= 5, -6 <= alpha, beta <= 6 deg); equations 1-4 (linear slender-body-theory models C_N = C_N0 + C_N_alpha alpha, C_Y = C_Y0 + C_Y_beta beta, C_m = C_m0 + (x_CN - x_BMC) C_N, C_n = C_n0 + (x_CY - x_BMC) C_Y, least squares per Mach number). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20160007695>
- **BRAUN-2017-SRP**: R. D. Braun, B. Sforzo and C. H. Campbell, *Advancing Supersonic Retropropulsion Using Mars-Relevant Flight Data: An Overview* (NTRS 20170008535). States that SpaceX performed the first supersonic retropropulsion manoeuvre on the Falcon 9 first stage on 29 September 2013 and that NASA and SpaceX formed a three-year (FY2014-) public-private partnership on SRP data analysis (flight reconstruction, CFD, imagery, Mars EDL design). Gives no implementable relations; cited only to document scope. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/20170008535>
- **NASA-GRC-DYNAMIC-PRESSURE**: NASA Glenn Research Center, *Dynamic Pressure*, Beginner's Guide to Aeronautics: `q = rho u^2 / 2`. Accessed 2026-10-25. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/dynamic-pressure/>
- **NASA-GRC-SPEED-OF-SOUND**: NASA Glenn Research Center, *Speed of Sound*, Beginner's Guide to Aeronautics: `a = sqrt(gamma R T)`, absolute temperature; Mach number as the ratio of speed to the speed of sound. Accessed 2026-10-25. <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/speed-of-sound-interactive/>
- **NASA-NPR-7123.1D**: NASA Procedural Requirements NPR 7123.1D, *NASA Systems Engineering Processes and Requirements*, Appendix E "Technology Readiness Levels" (TRL 4–6 definitions, hardware and software descriptions, success criteria), effective 2023-07-05. Accessed 2026-10-25. <https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE>

## Conduction heat transfer (verified 2026-10-25)

- **AHTT-V6**: J. H. Lienhard IV and J. H. Lienhard V, *A Heat Transfer Textbook*, 6th ed., version 6.00 (19 April 2024), Phlogiston Press; Chs. 1, 2, 4 and 5 (Eqs. (1.19)-(1.35), (2.10)-(2.36), (4.27)-(4.61), (5.1)-(5.62); Tables 2.1, 5.1, 5.2, 5.3; Examples 1.4-1.6, 2.1-2.10, 4.8-4.10, 5.2-5.7; Problems 2.1, 2.2, 2.15, 2.20, 2.21). Printed page = PDF page - 14. Accessed 2026-10-25. <https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf> (landing page <https://ahtt.mit.edu/>)
- **ENCIT-2012-KIRCHHOFF**: W. F. da Silva, E. D. Correa, R. M. Saldanha da Gama and M. L. Martins-Costa, "An a priori upper bound estimate for steady-state conduction heat transfer problems with a linear relationship between temperature and conductivity", Proc. ENCIT 2012, paper ENCIT2012-352, Sec. 3, Eqs. (5)-(6) (Kirchhoff transform, citing Arpaci 1966). Accessed 2026-10-25. <https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf>
- **DLMF-10-BESSEL**: NIST Digital Library of Mathematical Functions, Ch. 10 Bessel functions: Eqs. 10.2.2 (J series), 10.6.3 and 10.6.6 (derivatives), 10.9.1-10.9.2 (Bessel's integral), 10.17.1 (a_k), 10.25.1-10.25.2 (I series), 10.28.2 (Wronskian), 10.29.3 (I_0' = I_1, K_0' = -K_1), 10.31.1-10.31.2 (K_n series), 10.32.3 and 10.32.9 (integrals), 10.40.1-10.40.2 (large-argument expansions and the real-argument remainder bound). Accessed 2026-10-25. <https://dlmf.nist.gov/10>
- **DLMF-7-ERFC**: NIST DLMF Sec. 7.12(i), Eq. 7.12.1 (asymptotic erfc and its remainder bound for real z), and Sec. 7.18, Eqs. 7.18.1, 7.18.2, 7.18.7 (repeated integrals of erfc). Accessed 2026-10-25. <https://dlmf.nist.gov/7.18>, <https://dlmf.nist.gov/7.12>
- **DLMF-5-GAMMA**: NIST DLMF Eq. 5.2.3 (Euler's constant 0.57721 56649 01532 86060...) and Eq. 5.4.14 (psi(n+1) = H_n - gamma). Accessed 2026-10-25. <https://dlmf.nist.gov/5.2>, <https://dlmf.nist.gov/5.4>
- **DLMF-3.5-TRAPEZOID**: NIST DLMF Sec. 3.5(i), Eqs. 3.5.4-3.5.5 (exponential convergence of the trapezoidal rule for periodic analytic integrands and for analytic integrands over the real line). Accessed 2026-10-25. <https://dlmf.nist.gov/3.5>
- **YOVANOVICH-1988-ANNULAR-FIN**: M. M. Yovanovich, J. R. Culham and T. F. Lemczyk, "Simplified Solutions to Circular Annular Fins with Contact Resistance and End Cooling", J. Thermophysics and Heat Transfer 2(2), 152-157, April 1988, Eqs. (26)-(41) (one-dimensional annular fin, edge condition, heat rate, ideal rate and efficiency). Accessed 2026-10-25. <https://uwaterloo.ca/microelectronics-heat-transfer-laboratory/sites/default/files/2024-03/mhtl88-1.pdf>
- **NACA-TR-158**: D. R. Harper 3d and W. B. Brown, "Mathematical Equations for Heat Conduction in the Fins of Air-Cooled Engines", NACA Report 158, 1923; printed p. 26 assumption (3) (edge loss by adding half the fin thickness to the width) and printed p. 21 Eq. (5). Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19930091223/downloads/19930091223.pdf>
- **NIST-CODATA-2022-STEFAN-BOLTZMANN**: NIST Physical Measurement Laboratory, CODATA 2022 value of the Stefan-Boltzmann constant, 5.670 374 419... x 10^-8 W m^-2 K^-4 (exact). Accessed 2026-10-25. <https://physics.nist.gov/cgi-bin/cuu/Value?sigma>
- **PYTHON-3.12-MATH-ERF**: Python 3.12 standard library, `math.erf` and `math.erfc`. Accessed 2026-10-25. <https://docs.python.org/3.12/library/math.html#math.erfc>
- **AANDS-UBC-PDF**: M. Abramowitz and I. A. Stegun (eds.), *Handbook of Mathematical Functions*, NBS AMS 55, PDF copy hosted at UBC. Consulted for Bessel tables; this copy omits the numerical tables of Ch. 9, so no value from it is used. Accessed 2026-10-25. <https://personal.math.ubc.ca/~cbm/aands/abramowitz_and_stegun.pdf>

## Convection heat transfer (verified 2026-10-25)

- **LIENHARD-AHTT-6.00**: J. H. Lienhard V and J. H. Lienhard IV, *A Heat Transfer Textbook*, 6th ed., version 6.00 (2024), Phlogiston Press, free PDF (same work as the `heat_radiation` draft entry; merge the locator lists). Used by `heat_transfer/convection.py`: Eqs. (6.58), (6.68) and Example 6.5 (printed pp. 308-311); Eqs. (6.102), (6.111)-(6.113) (pp. 322, 325-327); Eqs. (7.22)-(7.30), Examples 7.1-7.2, Fig. 7.4 (pp. 359-364); Eq. (7.31) and Table 7.2 (pp. 365-366); Eqs. (7.32)-(7.47), Example 7.3 (pp. 366-372); Eqs. (7.48a)-(7.50), Table 7.3, Example 7.4 (pp. 372-374); Eqs. (7.52)-(7.55), Example 7.5 (pp. 377-379); Example 7.6 (pp. 381-382); Eqs. (7.60)-(7.63), Example 7.7, Table 7.4 (pp. 382-387); Eqs. (7.65)-(7.68), Example 7.8 (pp. 390-392); Eqs. (8.13a)-(8.13b), Example 8.1 (pp. 419-421); Eqs. (8.56)-(8.62b), Example 8.6 (pp. 445-448); Eq. (9.4), Example 9.2 (pp. 482-484); Eq. (9.6a) and lambda_d1 = 2.72 cm (pp. 486, 488); Eqs. (9.11)-(9.12), Examples 9.5-9.6 (pp. 493-496). In the PDF opened, printed page = PDF page - 14 in Chapters 6-9. Accessed 2026-10-25. <https://ahtt.mit.edu/>, <https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf>
- **NASA-TN-D-4332-TAYLOR-1968**: M. F. Taylor, *Correlation of Local Heat-Transfer Coefficients for Single-Phase Turbulent Flow of Hydrogen in Tubes with Temperature Ratios to 23*, NASA TN D-4332, January 1968. Used: Summary (p. 1), Table I range (p. 2), Symbols (p. 3: Re_b, Re_f, T*), Eqs. (9), (13), (14) (Hess-Kunz as reported) (pp. 6-7), Fig. 2, Eqs. (19)-(20) (p. 12), Table III (p. 13), region-2 accuracy (p. 14). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19680004771>
- **NASA-CR-678-AEROJET-1967**: Aerojet-General Corporation, *Heat Transfer to Cryogenic Hydrogen Flowing Turbulently in Straight and Curved Tubes at High Heat Fluxes*, NASA CR-678, February 1967. Used: Summary (p. 1), Eq. (1) with T_f = (T_w + T_b)/2 (p. 2), Eq. (3) and the C_L table (p. 3), Ito's Eq. (4) (p. 4), curved-tube findings (pp. 53-54). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19670008827>
- **NASA-TM-1998-206313-WADEL**: M. F. Wadel, *Comparison of High Aspect Ratio Cooling Channel Designs for a Rocket Combustion Chamber With Development of an Optimized Design*, NASA/TM-1998-206313, 1998. Used: nomenclature (pp. vii-viii) and the RTE coolant-side model, Eqs. (6)-(13) (pp. 15-17); the coolant pressure drop as viscous (Darcy, Eq. (14), with station-averaged density, velocity and diameter; velocity Eq. (15)) plus momentum (Eq. (16)) terms, and RTE's stated omission of sudden area changes (printed pp. 17-18, PDF pp. 25-26, scanned text partly garbled; read 2026-10-26); channel designs varying along the chamber, Appendix A Tables A-I to A-XV. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19980017619>
- **NASA-SP125-COOLANT-SIDE**: D. K. Huzel and D. H. Huang, *Design of Liquid Propellant Rocket Engines*, NASA SP-125, 2nd ed., 1971. Used: coolant-side heat transfer, Fig. 4-27 and Eqs. (4-23)-(4-25), printed pp. 105-106 (PDF p. 115); Eq. (4-24) not implemented (coolant constant C2 not given). Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- **NASA-SP-8087-1972**: NASA, *Liquid Rocket Engine Fluid-Cooled Combustion Chambers*, NASA SP-8087, April 1972. Used: sec. 2.1.1.5 (Bartz and Hatch-Papell as the accepted methods, burnout data scatter, Titan margins; pp. 15-16), sec. 3.1.1.5.2 (80 % burnout-flux rule, curvature guidance; p. 61); Table XI (film-cooling procedure, p. 92) read but not implemented. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19730022965>
- **NASA-TN-D-130-HATCH-PAPELL-1959**: J. E. Hatch and S. S. Papell, *Use of a Theoretical Flow Model to Correlate Data for Film Cooling or Heating an Adiabatic Wall by Tangential Injection of Gases of Different Fluid Properties*, NASA TN D-130, November 1959. Used: Summary data ranges (p. 1), assumption (6) (p. 2), Eqs. (9)-(12a) (pp. 5-6), accuracy (pp. 6-8, 10), Eq. (13) and worked example (pp. 8-10), symbols (pp. 11-12). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19890068390>

## Thermal radiation (verified 2026-10-25)

- **NIST-CODATA-2022**: NIST Reference on Constants, Units and Uncertainty, CODATA 2022 recommended values: Planck constant, speed of light, Boltzmann constant (exact), Stefan-Boltzmann constant 5.670 374 419... e-8 W m^-2 K^-4, Wien wavelength displacement constant 2.897 771 955... e-3 m K, second radiation constant 1.438 776 877... e-2 m K. Accessed 2026-10-25. <https://physics.nist.gov/cuu/Constants/index.html>, <https://physics.nist.gov/cgi-bin/cuu/Value?sigma>
- **DLMF-25.5**: NIST Digital Library of Mathematical Functions: Eq. (25.5.1) (zeta(s) Gamma(s) as the Bose integral), Eq. (25.6.1) (zeta(4) = pi^4/90), Eq. (8.2.2) (upper incomplete gamma), Eqs. (8.4.8), (8.4.11) (Gamma(n+1, z) = n! e^-z e_n(z)), Eq. (24.2.1) and Table 24.2.1 (Bernoulli numbers). Used to derive the black-body fraction series. Accessed 2026-10-25. <https://dlmf.nist.gov/25.5>, <https://dlmf.nist.gov/25.6>, <https://dlmf.nist.gov/8.2>, <https://dlmf.nist.gov/8.4>, <https://dlmf.nist.gov/24.2>
- **HOWELL-CATALOG-C**: J. R. Howell, *A Catalog of Radiation Heat Transfer Configuration Factors* (web catalog, cited by Lienhard as ref. [10.6], 2nd ed. 2001), Section C factors C-11 (identical parallel rectangles), C-14 (perpendicular rectangles with a common edge), C-40 and C-41 (coaxial parallel disks), C-78 (cylinder interior to itself), C-111 (interior of a conical frustum to its base). Accessed 2026-10-25. <https://www.thermalradiation.net/indexCat.html>, <https://www.thermalradiation.net/calc/sectionc/C-41.html>, <https://www.thermalradiation.net/calc/sectionc/C-111.html>
- **BORDBAR-FRAGA-HOSTIKKA-2020**: H. Bordbar, G. C. Fraga and S. Hostikka, "An extended weighted-sum-of-gray-gases model to account for all CO2-H2O molar fraction ratios in thermal radiation", *Int. Commun. Heat Mass Transfer* 110 (2020) 104400; accepted manuscript (CC BY-NC-ND) in the Aalto repository. Used: Eqs. (1)-(3), Table 2 (pure CO2) and Table 3 (pure H2O) WSGG coefficients at 1 atm, T = 300-2400 K, L = 0.01-60 m (manuscript pp. 4-5), Table 4 double-integration errors (p. 10). Accessed 2026-10-25. <https://doi.org/10.1016/j.icheatmasstransfer.2019.104400>, <https://acris.aalto.fi/ws/portalfiles/portal/39751900/ENG_Bordbar_et_al_An_extended_weighted_sum_of_gray_gases_model_ICHMT.pdf>
- **SADEGHI-2021-FSJ**: H. Sadeghi, S. Hostikka, G. C. Fraga and H. Bordbar, "Weighted-sum-of-gray-gases models for non-gray thermal radiation of hydrocarbon fuel vapors, CH4, CO and soot", *Fire Safety Journal* 125 (2021) 103420, CC BY. Used: Eqs. (9)-(14) (WSGG form and superposition of species), Eq. (24) (soot Planck-mean fit), Table B9 (soot WSGG coefficients, T = 300-3000 K). Table B8 (CO) was read but not adopted. Accessed 2026-10-25. <https://doi.org/10.1016/j.firesaf.2021.103420>, <https://acris.aalto.fi/ws/portalfiles/portal/66772583/1_s2.0_S0379711221001612_main_1.pdf>

## Real-fluid transport properties (verified 2026-10-25)

- **LEMMON-JACOBSEN-2004-TRANSPORT**: E. W. Lemmon and R. T Jacobsen, *Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air*, Int. J. Thermophys. 25(1), 21-69, 2004, doi:10.1023/B:IJOT.0000022327.04529.f3. Used: abstract (validity over all liquid and vapor states), p. 22 (not valid when the mean free path is comparable to the container), Eqs. (1)-(3) and Table I p. 23-24, Eqs. (4)-(5) and Table II p. 24, Tables III-IV pp. 25-26, Eqs. (7)-(11) with R0 = 1.01, nu = 0.63, gamma = 1.2415, k = 1.380658e-23 J/K and the lambda_c = 0 rule p. 27, Table V check values p. 28, uncertainty summary p. 62. Born-digital PDF hosted by NIST TRC; text layer read and Eqs. (1)-(11), Tables I-II checked against the rendered pages. Accessed 2026-10-25. <https://trc.nist.gov/refprop/FAQ/NAO.PDF>
- **ASSAEL-2011-H2-CONDUCTIVITY**: M. J. Assael, J.-A. M. Assael, M. L. Huber, R. A. Perkins and Y. Takata, *Correlation of the Thermal Conductivity of Normal and Parahydrogen from the Triple Point to 1000 K and up to 100 MPa*, J. Phys. Chem. Ref. Data 40(3), 033101, 2011, doi:10.1063/1.3606499. Used: critical constants of the equations of state (033101-5 and -9), Eq. (2) and Table 2 (033101-5), Eqs. (3)-(7) with Gamma = 0.052, xi_0 = 1.5e-10 m, q_D^-1 = 4.0e-10 m and T_ref = 49.7175 K (033101-6), parahydrogen q_D^-1 = 5.0e-10 m and T_ref = 49.407 K, Table 5 (033101-10), Table 7 check values (033101-11). Tables 2 and 5 read from the rendered page images (the text layer drops signs and exponents). Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.3606499.pdf>
- **FRIEND-ELY-INGHAM-1989-CH4**: D. G. Friend, J. F. Ely and H. Ingham, *Thermophysical Properties of Methane*, J. Phys. Chem. Ref. Data 18(2), 583-638, 1989, doi:10.1063/1.555828. Used: abstract (ranges), p. 591 (data ranges), Table 1 p. 590, Eqs. (9)-(15a) and Tables 8-9 p. 596, Eqs. (16)-(22) and Table 10 p. 597, Eqs. (23)-(26) and Table 11 p. 598 (scaled region, not implemented), Table A1 p. 630 and Table A3 p. 632 (check values). Scanned reprint; coefficients read from rendered page images and compared with the OCR text. Eq. (16) prints delta*_sigma = 11 outside the vapor branch; only 1 reproduces Table A3, so 1 is used. Accessed 2026-10-25. <https://srd.nist.gov/jpcrdreprint/1.555828.pdf>
- **MUZNY-2022-H2-VISCOSITY-ERRATUM**: C. D. Muzny, M. L. Huber and A. F. Kazakov, *Erratum: Correlation for the Viscosity of Normal Hydrogen Obtained from Symbolic Regression* [J. Chem. Eng. Data 58, 969-979 (2013), doi:10.1021/je301273j], NIST author copy. Read: the three corrections (Avogadro constant 6.022137e23 /mol in Eq. (6), exponent -i in Eq. (7), scale density 90.909090909 kg/m3 in Eq. (9)) and the three check values (40 K at 0, 50 and 100 kg/m3). THE 2013 PAPER ITSELF WAS NOT OPENED (publisher page blocked; the NIST publication record links only the DOI). The implemented hydrogen viscosity coefficients are transcribed from COOLPROP-TRANSPORT-C00300E, and correctness rests on that transcription as verified against 32 NIST WebBook rows (16 per hydrogen form; largest relative deviation 4.6e-12) and the erratum check values (6.1e-6). Range (triple point to 1000 K, up to 200 MPa) from the NIST publication record <https://www.nist.gov/node/574776>. Accessed 2026-10-25. <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934278>
- **ARP-MCCARTY-FRIEND-1998-HE**: V. D. Arp, R. D. McCarty and D. G. Friend, *Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa*, NIST Technical Note 1334 (revised), September 1998. Used: Sec. 1 (helium I code valid 2.5-1500 K to 100 MPa, p. 1), Sec. 4.1 (conductivity from Hands and Arp, validity unknown below 3.5 K or above 100 MPa, p. 4), Sec. 4.2 Eqs. (1)-(9) with the 100-110 K linear blend of Steward and Maitland-Smith dilute values and the 300 K freeze of the excess function (pp. 4-5). Scanned report: the OCR text is garbled in places, so Eqs. (1)-(9) were read from the rendered page images; they agree with CoolProp's `viscosity_helium_hardcoded` except that CoolProp truncates the last coefficient of Eq. (2) (0.012520841 versus 0.0125208416) and omits the 100-110 K blend; with the printed values the WebBook is matched to 5.5e-10. Accessed 2026-10-25. <https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf>
- **NIST-WEBBOOK-FLUID-TRANSPORT-2026**: E. W. Lemmon, M. O. McLinden and D. G. Friend, *Thermophysical Properties of Fluid Systems*, NIST Chemistry WebBook SRD 69 — `fluid.cgi` Action=Data, Digits=12, Type=IsoBar single-state rows with the viscosity and thermal-conductivity columns for parahydrogen (B5000001), normal hydrogen (C1333740), oxygen (C7782447), nitrogen (C7727379), methane (C74828) and helium (C7440597): 76 rows, each URL stored in `tests/reference_data/nist_real_fluid_transport.json`; transport-model citations read from the "References and Notes" of each fluid page (Action=Load). Accessed 2026-10-25. <https://webbook.nist.gov/chemistry/fluid/> <https://webbook.nist.gov/cgi/fluid.cgi>
- **COOLPROP-TRANSPORT-C00300E**: CoolProp developers, commit c00300e97ad4cdb30032635d7379209d5e45a121: `dev/fluids/Hydrogen.json` and `dev/fluids/ParaHydrogen.json` (`TRANSPORT.viscosity`, BibTeX key Muzny-JCED-2013: dilute collision-integral coefficients, Rainwater-Friend initial-density coefficients, epsilon/k = 30.41 K, sigma = 0.297 nm) and `src/Backends/Helmholtz/TransportRoutines.cpp` (`viscosity_hydrogen_higher_order_hardcoded`, `viscosity_dilute_collision_integral`, `viscosity_initial_density_dependence_Rainwater_Friend`, `conductivity_hardcoded_helium`, `viscosity_helium_hardcoded`) and `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` (`calc_viscosity_background`, how the terms combine). PRIMARY SOURCE for the hydrogen viscosity (Muzny et al. 2013) and the helium conductivity (Hands and Arp, Cryogenics 21, 697-703, 1981), whose papers were not opened; correctness rests on these transcriptions as verified against NIST WebBook rows (hydrogen viscosity 4.6e-12, helium conductivity 9.0e-7). Cross-check only for the helium viscosity (TN 1334). Accessed 2026-10-25. <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids/Hydrogen.json> <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids/ParaHydrogen.json> <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp> <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp>

## Finite-rate chemical kinetics (verified 2026-10-25)

- **NASA-RP-1328**: K. Radhakrishnan, "LSENS -- A General Chemical Kinetics and Sensitivity Analysis Code for Homogeneous Gas-Phase Reactions. I. Theory and Numerical Solution Procedures", NASA RP-1328, January 1994. Used: modified Arrhenius eq. (2.4), detailed balance and concentration equilibrium constant eqs. (2.6)-(2.11) (pp. 6-7), rates of progress, third bodies and species production eqs. (2.12)-(2.33) (pp. 8-11), one-dimensional steady reactive flow eqs. (2.57)-(2.85) (section 2.4, pp. 16-22). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950008642>
- **MFIX-25.2-CHEMKIN**: NETL, MFiX 25.2 documentation, section 4.14.4 "Implementation of CHEMKIN mechanism in MFiX": mass-action rates eqs. (4.4)-(4.7), Arrhenius eq. (4.8), reverse rates from K_c and K_p eqs. (4.9)-(4.13), third bodies eqs. (4.14)-(4.17), falloff and Troe form eqs. (4.18)-(4.20), CHEMKIN default units (cm-mole-sec-K, cal/mole) and DUPLICATE handling. Accessed 2026-10-25. <https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html>
- **CANTERA-3.2-RATES**: Cantera 3.2.0 documentation, "Rate Constant Parameterizations" (Troe falloff function with base-10 logarithms) and "Reaction Rates" (three-body and collider-specific reactions). Accessed 2026-10-25. <https://cantera.org/stable/reference/kinetics/rate-constants.html>, <https://cantera.org/stable/reference/kinetics/reaction-rates.html>
- **GRI-MECH-3.0**: G. P. Smith, D. M. Golden, M. Frenklach, N. W. Moriarty, B. Eiteneer, M. Goldenberg, C. T. Bowman, R. K. Hanson, S. Song, W. C. Gardiner Jr., V. V. Lissianski and Z. Qin, "GRI-Mech 3.0": documentation page (design target "natural gas combustion", 325 reactions, 53 species) and mechanism file grimech30.dat (SHA-256 a7726208619d27b50e43b9fb6b27dc0a737003c049cc79e0a849bc75301cd763), H2-O2 reactions on lines 22-28, 58-76, 151-161, 191-194, 394-395. Accessed 2026-10-25. <http://combustion.berkeley.edu/gri-mech/version30/text30.html>, <http://combustion.berkeley.edu/gri-mech/version30/files30/grimech30.dat>
- **HAIRER-RODAS**: E. Hairer and G. Wanner, Fortran code rodas.f (embedded Rosenbrock method of order 4(3), version of October 28, 1996; subroutines RODAS, ROSCOR, ROCOE with METH=1 coefficients on lines 895-941) and dc_decsol.f (DECOMR, SLVROD). Accessed 2026-10-25. <https://www.unige.ch/~hairer/prog/stiff/rodas.f>, <https://www.unige.ch/~hairer/prog/stiff/dc_decsol.f>
- **HAIRER-TESTSET**: E. Hairer, "Testset and Drivers" for stiff ODEs: ROBER (Robertson) and HIRES equations (equation.f), RODAS drivers (driver_rodas.f: initial values, output points) and reference solutions (res_exact_pic). Accessed 2026-10-25. <https://www.unige.ch/~hairer/testset/testset.html>
- **BRAY-1960-THESIS**: K. N. C. Bray, "Atomic recombination in a hypersonic wind-tunnel nozzle", PhD thesis, University of Southampton, 1960; section 5 "Approximate solutions", conditions (34)-(35) and the sudden-freezing criterion eq. (36), K of order one, experimental range 0.42-1.25 (Wegener), pp. 35-41. Accessed 2026-10-25. <https://eprints.soton.ac.uk/458375/1/38084.pdf>
- **NASA-RP-1311-PART1-EQ**: NASA RP-1311 Part I (S. Gordon and B. J. McBride, 1994), equations used by the finite-rate nozzle model: mixture enthalpy (2.14), entropy (2.16)-(2.17), specific impulse and vacuum specific impulse (6.8)-(6.9), area ratio (6.14). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-CEA-PY-RP1311-EX12**: NASA Glenn, cea 3.3.4 Python package sample samples/rp1311/example12.py ("frozen composition from the throat", n_frz = 2), used to set the frozen-station argument of the CEA bound generator. Accessed 2026-10-25. <https://pypi.org/project/cea/>

## Combustion efficiency and finite-area combustor (verified 2026-10-25)

- **NASA-RP1311-P1-FAC**: S. Gordon and B. J. McBride, NASA RP-1311 Part I, *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications: Analysis* (1994), chapter 6, sections 6.1-6.4 (printed pp. 25-31, PDF pp. 31-37): FAC station definitions (Fig. 6.1), assumptions (6.1), conservation equations (6.1)-(6.4), velocity (6.5), Isp (6.7)-(6.9), c* (6.11), momentum balance (6.26)-(6.28), tolerance (6.29), initial estimate (6.30)-(6.31), updates (6.32)-(6.33); frozen performance section 6.5. PDF SHA-256 6241ca2eece86c36edf82bc400039bc83f137223ebf1db7e3eb0186e6b23d097. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-RP1311-P1-EQRATIO**: Gordon and McBride, NASA RP-1311 Part I (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43 (PDF pp. 48-49): total-reactant molecular weight and density (9.11)-(9.13), oxidation-state totals (9.14)-(9.17), chemical equivalence ratio r (9.18), weight equivalence ratio phi (9.19) and the statement that r and phi coincide when all positive valences are in the fuel and all negative ones in the oxidizer. Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19950013764/downloads/19950013764.pdf>
- **NASA-RP1311-P1-ENTHALPY**: Gordon and McBride, NASA RP-1311 Part I (1994), section 4.1, equations (4.1)-(4.2), printed p. 19 (PDF p. 25): assigned enthalpy H(T) = Delta_f H(298.15) + [H(T) - H(298.15)]. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-CEA-SRC-VALENCE**: NASA Glenn Research Center, CEA source files `source/atomic_data.f90` (SHA-256 a959e5fbf16637ba48b7058ca1eb0d4234ae16eaf7d0dc200c35b4340b8be87b; valences H 1, C 4, N 0, O -2 in the `AtomicData` table) and `source/mixture.f90` (SHA-256 920b01cd478ff56587f9e7cdff064fea533941dc03b11c8d39fc41ce3c2f97c6; functions `mixture_chem_eq_ratio_to_of_ratio`, `mixture_of_ratio_to_chem_eq_ratio`, `mixture_weight_eq_ratio_to_of_ratio`, `mixture_get_valence`), commit 2f79a647e85737742fb08a127da00d1ea1d069da. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/mixture.f90>
- **NASA-CEA-SRC-FAC**: NASA Glenn Research Center, CEA source file `source/rocket.f90` (SHA-256 0c8bc2d820d7decbeb942fcb751f67c86396d8f3a4d15dca05a47b5e738ca16e), subroutine `RocketSolver_solve_fac`, commit 2f79a647e85737742fb08a127da00d1ea1d069da: station order injector/infinity/combustor end/throat, eq. (6.30) start, eq. (6.32) update and the mdot/Ac option. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/rocket.f90>
- **NASA-CEA2-EXAMPLE9-10-OUT**: NASA CEA2 regression outputs of RP-1311 examples 9 (FAC, Ac/At = 1.58) and 10 (FAC, mdot/Ac = 1333.9 kg/(s m^2)), `test/main_interface/reference_output/example9.out` (SHA-256 0fabc39640cadb5edead23275dcd269ea4f8fa9cc9971521f50a9ef784157c27, lines 126-160) and `example10.out` (SHA-256 9d3a2274529f01cee33cf2f88419173efb66496a3f265296b653a85d04fe2663, lines 160-194), same commit. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/test/main_interface/reference_output/example9.out>
- **NASA-CEA-PY-FAC**: NASA Glenn Research Center, `cea` Python package 3.3.4, `RocketSolver.solve(..., iac=False, ac_at=..., mdot=...)` and `Mixture.chem_eq_ratio_to_of_ratio` / `weight_eq_ratio_to_of_ratio`; packaged `thermo.lib` SHA-256 8e5df1cca92d4a48663d1ee5a1372e6508c59cddc2247ceeca32f041a03ec52a. Used only to generate `tests/reference_data/cea_fac_points.json` and `tests/reference_data/combustion_efficiency.json`. Accessed 2026-10-25. <https://pypi.org/project/cea/>
- **NACA-TN4219**: M. F. Heidmann and R. J. Priem, *Propellant Vaporization as a Criterion for Rocket Engine Design; Relation Between Percentage of Propellant Vaporized and Engine Performance*, NACA TN 4219, March 1958, equations (1)-(7), printed pp. 3-4 (PDF pp. 4-5). PDF SHA-256 b5f35aeb7211dd63856c902a7bbc2473323ab01ccf528332a35fbc9d2a300049. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930085136>
- **NASA-CR187109-ROCCID**: J. A. Muss, T. V. Nguyen and C. W. Johnson, *User's Manual for Rocket Combustor Interactive Design (ROCCID) and Analysis Computer Program, Volume I*, NASA CR-187109, May 1991, section 2.2, equations (2.4a)-(2.11), printed pp. 16-18 (PDF pp. 24-26): two-streamtube Rupe mixing split (Ito and Calhoon), eta_m(E_m) fit, vaporization efficiency, delivered c*, c* efficiency and Isp-based energy-release efficiency of the JANNAF simplified procedure. Also equations (3.16)-(3.19), printed p. 40 (PDF p. 48), not implemented. PDF SHA-256 9e6b85487a9cf723fdb4f938bc2fa10a3a9671c91e817c44a9d627646200b3cf. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19910014917>
- **NASA-SP8089-MIXING**: NASA SP-8089, *Liquid Rocket Engine Injectors* (NASA Space Vehicle Design Criteria, March 1976), equation (1), printed p. 19 (PDF p. 28), and Table IV mixing factors, printed p. 20 (PDF p. 29). Scanned PDF (no text layer) read as page images; SHA-256 5e9db6a6f7860ee51acbf5cbe137005ef259ec3fd87ff0ed5f1ad466d6227c44 (14,239,464 bytes). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19760023196>
- **NASA-SP125-LSTAR**: D. K. Huzel and D. H. Huang, *Design of Liquid Propellant Rocket Engines*, NASA SP-125, Rocketdyne / NASA, 1967, equations (4-3) and (4-4), printed p. 86 (PDF p. 95), and Table 4-1 recommended L* with stay-time range 0.002-0.040 s, printed p. 87 (PDF p. 96). PDF SHA-256 471e590a0c8e99c36a98132104a3ddcc7097b0b7c561fb2638ce9d3c80e26c87 (47,789,998 bytes, equal to the NTRS Content-Length). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19710019929>
- **NASA-SP125-DENSITY-IMPULSE**: Huzel and Huang, NASA SP-125, equation (1-45) density impulse = I_s d, printed p. 20 (PDF p. 29); Table 1-7 rows of I_s, d and I_sd, printed p. 25 (PDF p. 34). Equation (1-46) prints `(r_w - 1)`, which contradicts RP-1311 (9.13); not used. Accessed 2026-10-25. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- **NASA-TP1999-208522**: R. S. Jankovsky, T. D. Smith and A. J. Pavli, *High-Area-Ratio Rocket Nozzle at High Combustion Chamber Pressure - Experimental and Analytical Validation*, NASA/TP-1999-208522, June 1999: equations (3)-(4) effective chamber pressure and momentum pressure loss, printed p. 12 (PDF p. 14); equations (5)-(11) c*, C_F, Isp and their efficiencies, printed p. 13 (PDF p. 15); JANNAF loss chain ODE -> ODK -> MOC -> TDK/MABL -> x eta_c*, printed p. 25 with Figure 20 (PDF p. 27) and Figure 24, printed p. 27 (PDF p. 29). PDF SHA-256 ac4a890c582de559405c93c936f1b7e24d979c342ecf565535e5999762f8e0f2. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19990051032>
- **NASA-TM89814-JANNAF**: K. J. Davidian, *Comparison of Two Procedures for Predicting Rocket Engine Nozzle Performance*, NASA TM-89814 (AIAA-87-2071), 1987, printed pp. 1-4: JANNAF ODE/ODK/BLIMP and TDK procedures; boundary-layer thrust decrement converted to an impulse decrement and subtracted from the ideal Isp (p. 4). PDF SHA-256 08635ebed711ebe80c9d15d2189191c107396adfd8f4ff8ab457db0d61074114. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19870010946>

## Condensed-phase equilibrium (verified 2026-10-25)

- **NASA-RP1311-P1-CONDENSED**: S. Gordon and B. J. McBride, NASA RP-1311 Part I, *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications: Analysis* (1994); PDF SHA-256 6241ca2eece86c36edf82bc400039bc83f137223ebf1db7e3eb0186e6b23d097. Read: section 2.1 (R = 8314.51 J/(kg-mole K)) and section 2.2, eqs. (2.1a)-(2.4b), condensed species of negligible volume, printed pp. 3-4 (PDF pp. 9-10); section 2.3, eqs. (2.5)-(2.17), printed pp. 4-6 (PDF pp. 10-12); section 2.3.1, eqs. (2.18)-(2.23), printed p. 6; section 2.3.2, reduced Gibbs iteration eqs. (2.24)-(2.28) and the tp/hp/sp equation table, printed p. 7 (PDF p. 13); section 2.6, eqs. (2.67)-(2.75), printed p. 12 (PDF p. 18); sections 3.1-3.3, initial estimates, SIZE = 18.420681, control factor eqs. (3.1)-(3.3), updates (3.4), convergence tests (3.5)-(3.6d), printed pp. 13-15 (PDF pp. 19-21); section 3.4, condensed inclusion test eq. (3.7) and vapor-pressure form (3.8), one new species per convergence, removal of negative amounts, printed pp. 15-16 (PDF pp. 21-22); section 3.5, phase-transition criteria (50 K rule) and eqs. (3.9)-(3.10), printed pp. 16-17 (PDF pp. 22-23); section 3.6, singularities and the n_j = 1e-6 reset, printed p. 17 (PDF p. 23); section 4.1, eqs. (4.1)-(4.2), printed p. 19 (PDF p. 25); chapter 6, sections 6.1-6.2.4, assumptions (zero temperature and velocity lags) and eqs. (6.1)-(6.7), printed pp. 25-26 (PDF pp. 31-32); section 6.3.3, eqs. (6.15)-(6.17), printed p. 28 (PDF p. 34). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-CEA-THERMO-CONDENSED**: NASA Glenn Research Center, CEA thermodynamic database `data/thermo.inp`, repository commit 2f79a647e85737742fb08a127da00d1ea1d069da, file SHA-256 fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363, Apache-2.0. Records extracted by `scripts/sync_nasa_condensed.py` into `thermochemistry/data/nasa_cea_condensed.json`: 186 non-ionic gas products and 16 condensed products made only of Al, C, Cl, H, N, O (gases already bundled excluded), and the reactants `AL(cr)` (0 J/mol) and `NH4CLO4(I)` (-295767 J/mol) at 298.15 K. Verified in the script: 196 heat-of-formation identities within 0.088 J/mol, Gibbs equality at six phase transitions (alumina at 2327 K to 3.1e-8 RT, ice/water at 273.15 K to 2.7e-4 RT), interval steps in H/RT and S/R at most 1.24e-3 except a source step of 3.03 at 2700 K in `ALN(L)`. Accessed 2026-10-25. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp>
- **NASA-CEA-PY-CONDENSED**: NASA Glenn Research Center, `cea` Python package 3.3.4, `RocketSolver` (infinite-area combustor, `iac=True`, `supar` area ratios) and `EqSolver` (`TP`, `HP`) with `Mixture(..., products_from_reactants=True)`; packaged `thermo.lib` SHA-256 8e5df1cca92d4a48663d1ee5a1372e6508c59cddc2247ceeca32f041a03ec52a; enthalpy and entropy reported in kJ/kg and kJ/(kg K). Used only to generate `tests/reference_data/cea_condensed_points.json` (`scripts/generate_cea_condensed_points.py`). Accessed 2026-10-25. <https://pypi.org/project/cea/>
- **NASA-CEA-REPO-CONDENSED**: NASA Glenn Research Center, CEA source repository (Apache-2.0), home of `thermo.inp` and of the `cea` Python package; the RP-1311 sample installed with package 3.3.4 as `cea/samples/rp1311/example5.py` uses the reactant names `AL(cr)` and `NH4CLO4(I)` adopted here. Accessed 2026-10-25. <https://github.com/nasa/cea>

## Transpiration cooling and hot-fire heat transfer (verified 2026-10-25)

- **NACA-TN-3208-MICKLEY-1954**: H. S. Mickley, R. C. Ross, A. L. Squyers and W. E. Stewart, *Heat, Mass, and Momentum Transfer for Flow Over a Flat Plate With Blowing or Suction*, NACA TN 3208, July 1954. Used: film-theory Eqs. (21)-(30) (pp. 14-18), footnote 2 (p. 17), Table VII(a) theoretical and measured coefficients (pp. 78-79), porous-plate friction remark (pp. 53-54). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930084056>
- **NACA-REPORT-1182-ECKERT-LIVINGOOD-1954**: E. R. G. Eckert and J. N. B. Livingood, *Comparison of Effectiveness of Convection-, Transpiration-, and Film-Cooling Methods With Air as Coolant*, NACA Report 1182, 1954. Used: Eqs. (8), (20)-(31), Fig. 7 and the numerical-results discussion (pp. 8-11). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930092205>
- **NASA-CR-76367-KELLEY-LECUYER-1966**: J. B. Kelley and M. R. L'Ecuyer, *Transpiration Cooling - Its Theory and Application*, Purdue University JPC TM-66-5 / NASA CR-76367, June 1966. Used: Eqs. (2.14)-(2.17) (pp. 16-17, Rannie's analysis as restated), notation (Appendix A, pp. 143-146), sample calculations and Table B.1 (Appendix B, pp. 147-152). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19660021566>
- **NASA-CR-107268-PWA-1969**: Pratt & Whitney Aircraft, *Investigation of Methods for Transpiration Cooling Liquid Rocket Chambers*, PWA FR-3390 / NASA CR-107268, 1969. Used: Eq. (III-9) (p. III-12) and the two-temperature porous-wall solution Eqs. (III-11)-(III-29) (pp. III-13 to III-16). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19700004126>
- **NASA-TN-D-6396-SCHACHT-QUENTMEYER-1971**: R. L. Schacht and R. J. Quentmeyer, *Axial and Circumferential Variations of Hot-Gas-Side Heat-Transfer Rates in a Hydrogen-Oxygen Rocket*, NASA TN D-6396, July 1971. Used: Fig. 1 geometry (p. 3), Table I measured coefficients, wall and adiabatic-wall temperatures (pp. 10-11). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19710021262>
- **NASA-TN-D-2832-SCHACHT-1965**: R. L. Schacht, R. J. Quentmeyer and W. L. Jones, *Experimental Investigation of Hot-Gas Side Heat-Transfer Rates for a Hydrogen-Oxygen Rocket*, NASA TN D-2832, 1965. Used: correlation constants by station (table, p. 17). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19650016811>

## Combustion instability and propellant vaporization (verified 2026-10-25)

- **NASA-SP-194-1972**: D. T. Harrje and F. H. Reardon (eds.), *Liquid Propellant Rocket Combustion Instability*, NASA SP-194, 1972 (NTRS 19720026079). Used: Table 3.5.1 and Eqs. (3.5.1-18)-(3.5.1-19) (pp. 143-144), Fig. 4.2.2a and Eq. (4.2.2-23) (pp. 180-182), the 1L mean-flow factor (p. 182), the response factor (pp. 208-210), Eqs. (5.3.1-17)-(5.3.1-19) (p. 239), Example 5.5.1.1 and Eq. (5.5.1-6) (pp. 252-253), Eqs. (8.3.1-1a)-(8.3.1-2) (p. 401), (8.3.2-1)-(8.3.2-2) (pp. 410-411), (8.3.3-1) (p. 415), Sect. 9.6.6.2 (pp. 523-524). The 1T2R entry of Table 3.5.1 is misprinted (8.5263 for 8.5363). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19720026079>
- **NACA-TN-3985-PRIEM-1957**: R. J. Priem, *Propellant Vaporization as a Criterion for Rocket Engine Design; Calculations of Chamber Length to Vaporize a Single n-Heptane Drop*, NACA TN 3985, July 1957. Used: single-drop model Eqs. (1)-(12) (pp. 2-4). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930084710>
- **NACA-TN-3883-PRIEM-1958**: R. J. Priem, *Propellant Vaporization as a Criterion for Rocket-Engine Design; Calculations of Chamber Length to Vaporize Various Propellants*, NACA TN 3883, September 1958. Used: effective length (p. 5), Table I. Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19930084722>
- **NASA-CR-168066-AEROJET-1983**: J. Fang, *LOX/Hydrocarbon Fuel Carbon Formation and Mixing Data Analysis*, Aerojet Liquid Rocket Company, NASA CR-168066, 1983. Used: generalized-length eq. (1) (p. 38, attributed to Priem and Heidmann, NASA TR R-67, 1960, which is not on NTRS and was not read), Fig. 26 (p. 39), program listing with `DATA LG/EV`, the `CVO` statement and `SUBROUTINE SINTP` (pp. C-6 to C-8). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19830023615>

## Burnout and further hot-fire data (verified 2026-10-25)

- **NASA-TM-1998-206612-MEYER**: M. L. Meyer, D. L. Linne and D. C. Rousar, *Forced Convection Boiling and Critical Heat Flux of Ethanol in Electrically Heated Tube Tests*, NASA/TM-1998-206612 (AIAA-98-1055), January 1998. Used: Eq. (1) and Table III (pp. 5-6), Eqs. (2)-(3) and the high-pressure transition remarks (pp. 5 and 8). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19980017166>
- **NASA-CR-134282-TOBIN-1974**: R. D. Tobin, *Space Shuttle Orbit Maneuvering Engine Reusable Thrust Chamber: Adverse Operating Conditions Test Report*, Rocketdyne ASR-74-131, NASA CR-134282, 1974. Used: Eq. (A-1), Fig. A1 and p. A3 (Appendix). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19740018134>
- **NASA-TN-D-3532-PAVLI-1966**: A. J. Pavli, J. K. Curley, P. A. Masters and R. M. Schwartz, *Design and Cooling Performance of a Dump-Cooled Rocket Engine*, NASA TN D-3532, August 1966. Used: Fig. 3 contour (p. 3), Table I (p. 9), Eq. (1) (p. 11), Fig. 15 (p. 20), Table II (p. 44). Accessed 2026-10-25. <https://ntrs.nasa.gov/citations/19660022898>
- **NASA-SP-8087-BURNOUT-POINTER**: NASA SP-8087 (1972): the pages read (pp. 16-17, 61, 97) give no burnout correlation; Table IV (p. 17) grades the burnout-correlation quality per coolant and refers to Van Huff and Rousar, CPIA Publ. 121 (1966), ref. 34 on p. 97, marked Confidential and not on NTRS. <https://ntrs.nasa.gov/citations/19730022965>

- **NIST-DLMF-4.13-LAMBERT-W**: NIST Digital Library of Mathematical Functions, §4.13, Eq. (4.13.1), defining `W exp(W)=z` and its complex branches. Used for the independent delay-characteristic root witness in `test_feed_coupled_instability.py`; the root transformation is derived in the review rather than attributed to NIST. Accessed 2026-10-27. <https://dlmf.nist.gov/4.13>

- **CANTERA-32-IDEAL-MIXTURE**: IdealGasPhase, frozen ideal gas standard/partial molar thermodynamic relations. Cross-implementation comparison shares NASA coefficients and is not empirical validation. <https://cantera.org/3.2/cxx/d7/dfa/classCantera_1_1IdealGasPhase.html>
