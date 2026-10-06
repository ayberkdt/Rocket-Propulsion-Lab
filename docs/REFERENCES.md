# Engineering references

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

## Finite-burn, transient, and blowdown references

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

