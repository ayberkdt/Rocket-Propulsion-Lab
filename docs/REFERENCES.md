# Engineering references

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
- NASA SP-8041, **Liquid Rocket Engine Nozzles** — nozzle design and empirical
  separation context: <https://ntrs.nasa.gov/citations/19710021390>
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

- **SCIPY-1.14.1-IVP**: SciPy, *solve_ivp*. Event sign-change detection, terminal and direction attributes; multiple roots within one step may be missed. Accessed 2026-10-18. <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.solve_ivp.html>

- **PYTHON-3.12-FRACTION**: Python, *fractions*. Float inputs retain their exact binary value as a rational. Accessed 2026-10-18. <https://docs.python.org/3.12/library/fractions.html>
- **PYTHON-3.12-NEXTAFTER**: Python, *math.nextafter*. Adjacent representable floating-point values for the direct inventory helper's safe time horizon. Accessed 2026-10-19. <https://docs.python.org/3.12/library/math.html#math.nextafter>
- **PYTHON-3.12-DECIMAL**: Python, *decimal*. Precision contexts and exact float conversion for independent numerical checks. Accessed 2026-10-18. <https://docs.python.org/3.12/library/decimal.html>

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
