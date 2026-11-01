# Sentinel-1 / Sentinel-2: nine-slide English transcript

V6 · 1 November 2026. Source notes belong to the preceding slide.

## Slide 1: Mission Overview

- Sentinel-1 and Sentinel-2 | A/B/C/D comparison
- Category
- Sentinel-1
- Sentinel-2
- Purpose
- All-weather, day/night radar imaging
- Land and coastal multispectral imaging
- Operator / prime
- ESA / Thales Alenia Space
- ESA / Airbus Defence and Space
- Payload
- C-band SAR; AIS added on C/D
- MSI: 13 bands, 10 / 20 / 60 m sampling
- Reference orbit
- 693 km, Sun-synchronous; 12-day cycle
- 786 km, Sun-synchronous; 10-day cycle
- Design life
- 7 years minimum
- 7 years minimum
- C-vehicle calculation basis
- 2,186 kg launch; 154 kg hydrazine
- 1,143 kg launch; 133 kg hydrazine
- Series status¹
- A retired; B failed; C and D in service
- A supports calibration; B and C in service; D stored
- Sources: ESA mission facts; VV25 / VV24 launch kits; ESA fuelling records. ¹Status cut: 1 November 2026.

### Source notes

C flight-kit masses are used consistently with C fuelling records. Generic family fact sheets give 2,185 and 1,140 kg; these are not exact confirmations of the C values. Status references are retained in the research source archive. S1D resumed nominal operations on 29 October 2026 after an anomaly; the public cause is not established.

Sources (locators and limits are in PROPULSION_CATALOGUE):
FACT1: Sentinel-1 Facts and figures
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures

FACT2: Sentinel-2 Facts and figures
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures

KIT1C: VV25 launch kit
https://sky-brokers.com/wp-content/uploads/2025/07/Arianespace-Launchkit-VV25-ESA-Copernicus-Sentinel-1c.pdf

KIT2C: VV24 launch kit
https://newsroom.arianespace.com/download/?n=b694039895faf86c0a7a5bcab504d9a1.pdf&id=39971

F1C: Sentinel-1C fuelled ahead of liftoff, 18 November 2024
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff

F2C: Sentinel-2C fully loaded, 19 August 2024
https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/


## Slide 2: Mission Profile & Architecture

- The launch vehicle injects; onboard thrusters trim the orbit
- Flight chain
- Sentinel-1
- Sentinel-2
- A / B
- Soyuz-ST + Fregat; 2014 / 2016
- Vega + AVUM; 2015 / 2017
- C
- Vega-C + AVUM+; 5 Dec 2024
- Vega + AVUM; 5 Sep 2024
- D
- Ariane 6 + Vinci upper stage; 4 Nov 2025
- Ground storage; no launch to report
- All launched vehicles departed Kourou. Images: ESA launch-campaign archive. Launcher propulsion is separate from onboard hydrazine.

### Source notes

Flight architecture includes the launch vehicle and upper stage. Sentinel-1 and Sentinel-2 are Earth-orbiting spacecraft: no cruise stage, aeroshell or descent stage. Launcher dates and stage types are traceable in the accompanying series reports and launch kits. Images and their original URLs are in source_receipts.json.

Sources (locators and limits are in PROPULSION_CATALOGUE):
FACT1: Sentinel-1 Facts and figures
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures

FACT2: Sentinel-2 Facts and figures
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures

KIT1C: VV25 launch kit
https://sky-brokers.com/wp-content/uploads/2025/07/Arianespace-Launchkit-VV25-ESA-Copernicus-Sentinel-1c.pdf

KIT2C: VV24 launch kit
https://newsroom.arianespace.com/download/?n=b694039895faf86c0a7a5bcab504d9a1.pdf&id=39971


## Slide 3: System Budgets & Constraints

- Published system limits; complete flight budgets remain unavailable
- Budget / constraint
- Sentinel-1
- Sentinel-2
- Non-hydrazine mass
- 2,032 kg = 2,186 − 154
- 1,010 kg = 1,143 − 133
- Power, family design
- 5,900 W at end of life
- 1,700 W at end of life
- Pointing / orbit
- SAR repeat-track control: 100 m tube
- Knowledge: 10 µrad/axis (95.5%, 2014)
- Payload storage
- 1,410 Gbit, family design
- 2.4 Tbit, family design
- Payload downlink
- 520 Mbit/s aggregate, two channels
- 560 Mbit/s line rate; 520 effective in SP appendix
- Thermal / power
- Motor interface: +5 to +90°C (Rafael)
- Flight heater and valve power budget not published
- Unclosed allocations
- Subsystem mass, link margin, thermal loads
- Subsystem mass, link margin, thermal loads
- Sources: ESA SP-1322/1 and /2; Airbus AOCS 2014; Rafael p46. Subtracted masses are not certified dry masses. Older design values are labelled.

### Source notes

These are published constraints, not a reconstructed complete budget. The non-hydrazine remainder includes all other mass, including pressurant and residuals. Family power/memory/rate requirements do not certify the C/D configuration. S2 SP main text and appendix give different line/effective rates. Do not silently resolve them to one number. The full pointing confidence levels, frequency bands and thermal context are in the series reports.

Sources (locators and limits are in PROPULSION_CATALOGUE):
S1SP: SP-1322/1, 2012
https://sentinels.copernicus.eu/documents/247904/349449/S1_SP-1322_1.pdf

S2SP: SP-1322/2, 2012
https://sentinels.copernicus.eu/documents/247904/349490/S2_SP-1322_2.pdf

S214: Sentinel-2 AOCS: Challenges and Solutions, GNC 2014
https://www.researchgate.net/publication/271584776_THE_SENTINEL-2_SATELLITE_ATTITUDE_CONTROL_SYSTEM_-_CHALLENGES_AND_SOLUTIONS

R24: Space Propulsion Technical Catalogue, 2024
https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf


## Slide 4: Ayberk Demirkanat

- Concept of Operations
- Orbit control continues from acquisition to disposal
- Phase
- Required action
- Propulsion accounting
- Launch and separation
- Upper stage supplies the main injection
- Keep launcher impulse outside spacecraft Δv
- LEOP and acquisition
- Deploy; acquire attitude; trim injection errors
- Sum each in-plane and out-of-plane manoeuvre
- Routine observation
- Control ground track, eccentricity and inclination
- Track actual pressure, pulse count and impulse
- Collision avoidance
- Execute avoidance; restore the reference orbit
- Add both burns by magnitude
- End of life
- Lower orbit; vent or isolate stored energy
- Reserve fuel for disposal and passivation
- Sentinel-1B case: final 580 × 610 km orbit; partial passivation in September 2024.
- Sources: ESA flight dynamics (2017); SpaceOps 2025 #226, #360 and #617. Burn direction changes do not cancel fuel use.

### Source notes

Sentinel-1B had already ceased its imaging mission in 2022. End-of-life operations followed in 2024. The residual-fuel estimate (~70 kg) is uncertain, not a measured reserve for C/D. There was no dedicated dump valve; full depressurisation was not accomplished. Historical 25-year re-entry criteria must not be described as a new five-year compliance certificate.

Sources (locators and limits are in PROPULSION_CATALOGUE):
S1FD: Sentinel-1B LEOP and Reference Orbit Acquisition, 2017
https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf

S1226: Sentinel-1B End-Of-Life Operations, SpaceOps 2025 #226
https://publications.spaceops.org/2025/download_by_id.php?id=0226

S1360: Stop OCM Burn FDIR, SpaceOps 2025 #360
https://publications.spaceops.org/2025/download_by_id.php?id=0360

S2617: Optimising Tandem Orbit Control, SpaceOps 2025 #617
https://publications.spaceops.org/2025/download_by_id.php?id=0617


## Slide 5: Propulsion Design Drivers

- Precision and end-of-life pressure drive the design
- Driver
- Engineering consequence
- Evidence / limit
- Small orbit corrections
- 1 N class thrusters; short pulses and repeatable impulse
- Minimum impulse bit depends on pressure and duty
- Pointing during burns
- Thruster geometry and attitude control reduce useful impulse
- Installed count is not simultaneous axial thrust
- Blow-down feed
- Thrust and Isp change as tank pressure falls
- Tank pressure is not motor inlet pressure
- Plume impingement
- Use flight-specific efficiency and burn limits
- S1B study reports about 30% effective thrust loss
- Disposal and residuals
- Protect EOL pressure and manoeuvre reserves
- Expulsion efficiency alone does not prove usable fuel
- Qualification and heritage
- Respect pulse, impulse and thermal limits
- Flight heritage does not qualify a new student design
- Sources: Rafael catalogue; ArianeGroup brochure; ESA S1B studies. S1B limits cannot be assigned to all C/D burns.

### Source notes

S1B acquisition paper reports 35 s for plume-affected in-plane firings. A special 90-degree slew aligns OCT3 for large in-plane burns up to 600 s; six hours in nominal pointing restores thermal state and battery charge afterwards. These are 2016 procedure conditions, not a universal out-of-plane or thruster limit. Exact C/D cant angles, thermal recovery rules and efficiency maps were not recovered. No numerical preflight TRL is invented.

Sources (locators and limits are in PROPULSION_CATALOGUE):
R24: Space Propulsion Technical Catalogue, 2024
https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf

AGPDF: Hydrazine Thrusters brochure
https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf

S1FD: Sentinel-1B LEOP and Reference Orbit Acquisition, 2017
https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf

S1360: Stop OCM Burn FDIR, SpaceOps 2025 #360
https://publications.spaceops.org/2025/download_by_id.php?id=0360


## Slide 6: Propulsion Architecture

- Both use hydrazine; branch layouts and published hardware differ
- Item
- Sentinel-1
- Sentinel-2
- Propulsion type
- Hydrazine monopropellant; blow-down
- Helium-pressurised hydrazine; blow-down
- Installed thrusters
- 14 × 1 N; 2 redundant branches × 7
- 8 × 1 N; 2 cold-redundant sets × 4
- Manufacturer / model
- Rafael 1 N; Type 1/2 mix unverified
- EADS ST CHT1N-6 in A/B design paper¹
- Functions
- Per branch: 3 orbit + 4 attitude thrusters
- Tilted array for orbit control and burn attitude control
- Tank
- PEPT 590GB: 204 L; 17.5 kg dry (catalogue)
- Flight tank model / dry mass unavailable
- C loading
- 154 kg hydrazine; initial pressure unknown
- 133 kg hydrazine; He at 21 bar²
- System placement
- OCT mounts: +X / −X / −Y, A/B design
- Separate module in launcher-interface ring
- ¹C/D exact subtype not verified. ²ESA does not state absolute/gauge basis. Sources: ESA operations, Airbus AOCS, Rafael, ESA fuelling.

### Source notes

Independent authority checks: ESA operational architecture agrees with Rafael deliveries on the S1 14-thruster count and 1 N class. ESA and Airbus agree on S2 hydrazine architecture. The Airbus 2014 paper names CHT1N-6 in Table 5 for A/B design. Current ArianeGroup product figures describe the related 1 N family; they are not a C/D acceptance record. S1 PEPT specifications are product data linked by manufacturer heritage, not an as-built C tank certificate. S2C 21 bar is campaign tank pressure, not inlet pressure.

Sources (locators and limits are in PROPULSION_CATALOGUE):
R24: Space Propulsion Technical Catalogue, 2024
https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf

S1360: Stop OCM Burn FDIR, SpaceOps 2025 #360
https://publications.spaceops.org/2025/download_by_id.php?id=0360

S1FD: Sentinel-1B LEOP and Reference Orbit Acquisition, 2017
https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf

S214: Sentinel-2 AOCS: Challenges and Solutions, GNC 2014
https://www.researchgate.net/publication/271584776_THE_SENTINEL-2_SATELLITE_ATTITUDE_CONTROL_SYSTEM_-_CHALLENGES_AND_SOLUTIONS

S209: Sentinel-2 AOCS Model-Based Design, DGLR 2009
https://www.fzt.haw-hamburg.de/pers/Scholz/dglr/dlrk2009_mitReview/Papers/121179.pdf

F1C: Sentinel-1C fuelled ahead of liftoff, 18 November 2024
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff

F2C: Sentinel-2C fully loaded, 19 August 2024
https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/


## Slide 7: Propulsion Performance

- Manufacturer specifications | conditions must stay attached
- Parameter
- Rafael 1 N family
- ArianeGroup 1 N family
- Inlet pressure
- Type 1: 9–24.5; Type 2: 5.5–24.5 bar abs
- 5.5–22 bar; reference basis unspecified
- Specific impulse
- >214 s at 22 bar; >205 s at 5.5 bar
- 200–223 s range; 220 s nominal
- Thrust / minimum impulse
- BOL 1.3–1.4 N (T1), 1.0–1.1 N (T2); 0.008 N·s¹
- 0.320–1.1 N; 0.010–0.043 N·s
- Nozzle expansion / mass
- 130; ≤310 g including 1 m leads
- 80; 290 g
- Total impulse / pulses
- T1: 100 kN·s / 100,000; T2: 60 kN·s / 58,000
- 135 kN·s / 59,000 cycles
- Valve and heater
- Normally closed, dual coil/seat; 9.2 W valve at 28 V
- Two normally closed valves in series; redundant heater
- Thermal / life conditions
- 180°C recommended catalyst preheat; 15-year total life
- 12 h single / 50 h accumulated firing, qualified limits
- ¹0.008 N·s at 5.5 bar and 0.02 s / 1 s duty. Sources: Rafael p46; ArianeGroup pp5–6. These are not C/D flight acceptance maps.

### Source notes

All conditions and additional dimensions/leakage/power data are retained in the catalogue. ArianeGroup nominal 1 N, 220 s and 0.44 g/s are rounded product descriptors and do not form one exact matched operating point. NASA Table 4-2 repeats supplier data: agreement is not independent test verification. Rafael 180°C preheat concerns catalyst temperature, whereas +5 to +90°C concerns the operating interface; do not compare these as the same location.

Sources (locators and limits are in PROPULSION_CATALOGUE):
R24: Space Propulsion Technical Catalogue, 2024
https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf

AGPDF: Hydrazine Thrusters brochure
https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf

AGWEB: 1 N Monopropellant Hydrazine Thruster
https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html

NASA26: 2026 State of the Art, In-Space Propulsion
https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/


## Slide 8: Delta-V Analysis

- Separate published manoeuvres from ideal capability
- Δv = g₀ Isp η ln(m₀ / mf)      |      mp = m₀ [1 − exp(−Δv / (g₀ Isp η))]
- Calculation / record
- Sentinel-1
- Sentinel-2
- Ideal loaded-fuel capability²
- 146.86 m/s at Isp = 205 s
- 266.89 m/s at Isp = 220 s
- Geometric disposal comparison
- 693 km → 580 × 610 km: 52.58 m/s
- 786 km → 653 km: 70.22 m/s
- Historical flight / design context
- S1B 2024 final orbit; value above is not flown burn sum
- 2014 plan: 73 m/s, 47 kg, 116 × 600 s
- S2C published acquisition subset
- —
- Approximately 7.97 m/s, magnitude accounting
- S2C subset: ~6.33 + 0.012 + ~0.57 m/s, plus Table 3 manoeuvre magnitudes.
- It is not a complete lifetime Δv budget. No current C/D EOL allocation is published.
- ²η = 1; all loaded fuel idealised as usable; zero reserve. Sources: SpaceOps #617; Airbus 2014; reproducible calculation appendix.

### Source notes

g0 = 9.80665 m/s². Matched C masses and fills: S1C 2186/154 kg, S2C 1143/133 kg. η includes direction/plume efficiency only in explicitly labelled scenarios. Ideal capability assumes constant Isp and does not model residual pressure. S2C Table 3: vector magnitude for 6 November is sqrt(0.9344²+0.0537²); opposite-signed burns on 23 November cost 0.0207+0.0206 m/s. Total including summary phases is 7.9684418 m/s, with approximate source phases. Hohmann geometry is a two-impulse two-body comparison, not a reconstruction of finite burns.

Sources (locators and limits are in PROPULSION_CATALOGUE):
F1C: Sentinel-1C fuelled ahead of liftoff, 18 November 2024
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff

F2C: Sentinel-2C fully loaded, 19 August 2024
https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/

KIT1C: VV25 launch kit
https://sky-brokers.com/wp-content/uploads/2025/07/Arianespace-Launchkit-VV25-ESA-Copernicus-Sentinel-1c.pdf

KIT2C: VV24 launch kit
https://newsroom.arianespace.com/download/?n=b694039895faf86c0a7a5bcab504d9a1.pdf&id=39971

S2617: Optimising Tandem Orbit Control, SpaceOps 2025 #617
https://publications.spaceops.org/2025/download_by_id.php?id=0617

S214: Sentinel-2 AOCS: Challenges and Solutions, GNC 2014
https://www.researchgate.net/publication/271584776_THE_SENTINEL-2_SATELLITE_ATTITUDE_CONTROL_SYSTEM_-_CHALLENGES_AND_SOLUTIONS

S1226: Sentinel-1B End-Of-Life Operations, SpaceOps 2025 #226
https://publications.spaceops.org/2025/download_by_id.php?id=0226


## Slide 9: Propellant Budget

- One-metre-per-second comparison | explicit ideal assumptions
- Quantity
- Sentinel-1C example
- Sentinel-2C example
- Initial mass / Isp
- 2,186 kg / 205 s
- 1,143 kg / 220 s
- Fuel for Δv = 1 m/s
- 1.08710 kg
- 0.52967 kg
- Flow per assumed 1 N
- 0.49742 g/s
- 0.46351 g/s
- 1 N axial burn time
- 2,185.46 s
- 1,142.74 s
- 4 × 1 N axial burn time
- Geometry not assumed
- 285.68 s; illustrative, zero cant
- Fuel carried at launch
- 154 kg, ESA campaign record
- 133 kg, ESA campaign record
- Extra operating life
- Not determinable from loading alone
- Not determinable from loading alone
- S2C acquisition subset at 220 s, η = 1: 4.21 kg calculated fuel; 128.79 kg model remainder.
- Flight reserve needs pressure, residuals, future demand and disposal allocation.
- ṁ = F / (g₀ Isp); t = mp / ṁ. Constant thrust/Isp, no cant, no losses. Model remainder is not measured onboard fuel.

### Source notes

These are comparison scenarios, not commanded burn durations. At equal Isp and η, parallel thrusters shorten time but do not reduce ideal fuel. Cant reduces axial efficiency; 20° and 30° sensitivities and finite-burn integration are in CALCULATIONS_ENGLISH. A remaining mass cannot be converted to extra years without pressure access, an operational demand forecast and protected disposal/passivation allocations.

Sources (locators and limits are in PROPULSION_CATALOGUE):
F1C: Sentinel-1C fuelled ahead of liftoff, 18 November 2024
https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff

F2C: Sentinel-2C fully loaded, 19 August 2024
https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/

S2617: Optimising Tandem Orbit Control, SpaceOps 2025 #617
https://publications.spaceops.org/2025/download_by_id.php?id=0617

R24: Space Propulsion Technical Catalogue, 2024
https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf

AGPDF: Hydrazine Thrusters brochure
https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf
