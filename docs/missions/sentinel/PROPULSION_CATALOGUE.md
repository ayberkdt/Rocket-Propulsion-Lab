# Sentinel-1 / Sentinel-2 propulsion catalogue

**V6 · English · Source cut: 1 November 2026.**

Published Sentinel-1/2 first-generation A/B/C/D propulsion details. Missing flight data remain null. No exhaustive access to private acceptance documents is claimed.

Different official documents are compared by exact fact, vehicle, date, unit and condition. A repeated press release, supplier-derived NASA row or an older design value is not independent flight verification.

## Verification labels

- MULTI_AUTHORITY_AGREEMENT: the stated fact agrees across distinct official authorities; applicability and conditions still apply.
- SINGLE_PRIMARY_SOURCE: one primary authority/data chain was recovered. Repetition on that authority’s website does not make it independent.
- OFFICIAL_COMPILATION_CROSS_CHECK: a second official publisher reproduces supplier data; this is a consistency check, not an independent measurement.
- REPORTED_IN_MULTIPLE_DOCUMENTS: multiple documents agree, without an established independent measurement chain.
- SOURCE_DISAGREEMENT_OR_VERSION: numbers differ by vehicle, date or definition; the selected basis is stated.
- NOT_PUBLICLY_ESTABLISHED: no sufficient public source was recovered. Null is not zero.

## Critical-source checks

| ID | Fact / applicability | Value | Verification | Evidence and limit |
|---|---|---|---|---|
| S1-001 | Propulsion type / Sentinel-1 family; detailed operational evidence for B | Catalytic hydrazine monopropellant  | MULTI_AUTHORITY_AGREEMENT | S1FD: PDF p2; R24: printed pp43–46 / PDF pp46–49. ESA flight-dynamics description and manufacturer product-family identification agree. |
| S1-002 | Installed thruster count / A/B operational architecture; C/D manufacturer deliveries | 14 thrusters | MULTI_AUTHORITY_AGREEMENT | S1360: p2: two branches of seven; S1FD: p2: three OCT pairs and four RCT pairs; R24: printed p10: 28 for A/B and 28 delivered for C/D. Manufacturer delivery totals are consistent with 14 per spacecraft; operator paper states two times seven directly. This does not identify the Type 1/2 distribution. |
| S1-003 | Nominal thruster class / Sentinel-1 operational architecture | 1 N | MULTI_AUTHORITY_AGREEMENT | S1360: p2; S1FD: p2; R24: printed p43. Operator and manufacturer agree on the 1 N class. Actual thrust varies with pressure and condition. |
| S1-004 | Thruster supplier / Sentinel-1 A/B/C/D program selection | Rafael  | SINGLE_PRIMARY_SOURCE | R24: printed p43; heritage printed p10. Direct manufacturer statement. A separate official operator source explicitly naming Rafael was not recovered. |
| S1-016 | Type 1 inlet pressure / Rafael 1 N product family; flight unit subtypes not established | [9, 24.5] bar abs | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-017 | Type 2 inlet pressure / Rafael 1 N product family; flight unit subtypes not established | [5.5, 24.5] bar abs | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-021 | Steady-state Isp at 22 bar / Rafael 1 N product family; flight unit subtypes not established | >214 s | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-022 | Steady-state Isp at 5.5 bar / Rafael 1 N product family; flight unit subtypes not established | >205 s | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-023 | Minimum impulse bit / Rafael 1 N product family; flight unit subtypes not established | 0.008 N s | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-027 | Type 1 total delivered impulse / Rafael 1 N product family; flight unit subtypes not established | 100000 N s | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-028 | Type 2 total delivered impulse / Rafael 1 N product family; flight unit subtypes not established | 60000 N s | SINGLE_PRIMARY_SOURCE | R24: printed p46 / PDF p49. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S1-056 | Tank model / Sentinel-1 A/B/C/D manufacturer heritage | PEPT-590GB  | SINGLE_PRIMARY_SOURCE | R24: printed p9; p41 / PDF pp12,44. Direct manufacturer heritage; not a C/D signed as-built ICD. |
| S1C-080 | Loaded propellant / Sentinel-1C, November 2024 campaign | 154 kg | SINGLE_PRIMARY_SOURCE | F1C: fuel-loading paragraph; R24: printed p41 nominal product capacity, not the loading event. ESA directly confirms the flight loading. Rafael independently confirms a compatible nominal tank capacity, not the flight fill. No second independent fill measurement found. |
| S1C-081 | Propellant identity / Sentinel-1C | Hydrazine  | MULTI_AUTHORITY_AGREEMENT | R24: printed p43 A/B/C/D selection; S1SP: hydrazine propulsion description. Manufacturer family selection agrees with ESA mission architecture. |
| S1C-082 | Launch mass used for calculations / S1C VV25 launch-kit value | 2186 kg | SOURCE_DISAGREEMENT_OR_VERSION | KIT1C: p3; FACT1: generic family mass 2185 kg. The vehicle-specific kit says 2186 kg and was recovered on a third-party mirror. ESA general facts says 2185 kg. Exact two-authority agreement is not claimed. |
| S2-099 | Propulsion type / Sentinel-2 family | Catalytic hydrazine monopropellant  | MULTI_AUTHORITY_AGREEMENT | S2SP: printed pp15–16 / PDF pp23–24; S214: pp5,12; F2C: loading paragraph. ESA mission design, Airbus design paper and C flight loading agree on hydrazine monopropulsion. |
| S2-100 | Installed thruster count / S2 A/B design; C/D detailed acceptance configuration not public | 8 thrusters | REPORTED_IN_MULTIPLE_DOCUMENTS | S214: pp5,12: 2 cold-redundant branches of four; S209: PDF p2: two times four thrusters. Both primary design documents are from the Astrium/Airbus development chain. They are not independent measurements. |
| S2-101 | Specific thruster model / 2014 A/B AOCS equipment selection | CHT1N-6  | SINGLE_PRIMARY_SOURCE | S214: p9 Table 5. Primary author copy names EADS ST, Germany. Same C/D sub-version has not been independently established. |
| S2-102 | Supplier/product-family association / A/B design and manufacturer heritage | EADS ST / later ArianeGroup 1 N family  | REPORTED_IN_MULTIPLE_DOCUMENTS | S214: p9 Table 5; AGPDF: p6: S2A and S2B heritage. The A/B design and later manufacturer brochure identify the same supplier lineage. This is not evidence for an undisclosed C/D subtype. |
| S2-114 | Nominal Isp / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | 220 s | SINGLE_PRIMARY_SOURCE | AGPDF: p6; AGWEB: characteristics table. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S2-115 | Isp range / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | [200, 223] s | OFFICIAL_COMPILATION_CROSS_CHECK | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row. NASA and supplier values agree. NASA explicitly compiles manufacturer-provided data; independent experimental verification is not claimed. |
| S2-118 | Inlet pressure range / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | [5.5, 22] bar | SINGLE_PRIMARY_SOURCE | AGPDF: p6; AGWEB: characteristics table. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S2-119 | Minimum impulse bit / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | [0.01, 0.043] N s | SINGLE_PRIMARY_SOURCE | AGPDF: p6; AGWEB: characteristics table. Manufacturer specification. Independent flight acceptance/performance data for this parameter were not recovered. |
| S2-121 | Thruster mass with valves / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | 290 g | OFFICIAL_COMPILATION_CROSS_CHECK | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row. NASA and supplier values agree. NASA explicitly compiles manufacturer-provided data; independent experimental verification is not claimed. |
| S2-123 | Total impulse qualification / ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map | 135000 N s | OFFICIAL_COMPILATION_CROSS_CHECK | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row. NASA and supplier values agree. NASA explicitly compiles manufacturer-provided data; independent experimental verification is not claimed. |
| S2C-141 | Loaded hydrazine / 16 August 2024 flight campaign | 133 kg | SINGLE_PRIMARY_SOURCE | F2C: loading paragraph. ESA flight-campaign record is direct. Multiple ESA image captions repeat the same campaign statement, not independent confirmation. |
| S2C-142 | Pressurant / 16 August 2024 flight campaign | Helium  | SINGLE_PRIMARY_SOURCE | F2C: loading paragraph. Flight loading confirmed by ESA. A second independent acceptance record was not recovered. |
| S2C-143 | Tank loading pressure / 16 August 2024 flight campaign | 21 bar | SINGLE_PRIMARY_SOURCE | F2C: loading paragraph.  |
| S2C-144 | Launch mass used for calculations / VV24 vehicle-specific launch-kit value | 1143 kg | SOURCE_DISAGREEMENT_OR_VERSION | KIT2C: p3; FACT2: generic 1140 kg including 123 kg fuel. The flight kit is vehicle-specific. The general ESA fact-sheet value is not substituted for it. No exact two-source mass confirmation claimed. |
| S1-182 | Feed mode / Sentinel-1 family | Pressure-fed blow-down  | SINGLE_PRIMARY_SOURCE | S1FD: p2; R24: printed p41.  |
| S2-183 | Feed mode / A/B design; C helium loading confirmed | Helium-pressurised blow-down  | MULTI_AUTHORITY_AGREEMENT | S214: p12; F2C: loading paragraph. Airbus design explicitly specifies blow-down; ESA C campaign confirms helium pressurisation. A C pressure-performance map is unavailable. |

## S1: Architecture

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-001 | Propulsion type | Catalytic hydrazine monopropellant  | Sentinel-1 family; detailed operational evidence for B.  | S1FD: PDF p2; R24: printed pp43–46 / PDF pp46–49 |
| S1-002 | Installed thruster count | 14 thrusters | A/B operational architecture; C/D manufacturer deliveries.  | S1360: p2: two branches of seven; S1FD: p2: three OCT pairs and four RCT pairs; R24: printed p10: 28 for A/B and 28 delivered for C/D |
| S1-003 | Nominal thruster class | 1 N | Sentinel-1 operational architecture.  | S1360: p2; S1FD: p2; R24: printed p43 |
| S1-004 | Thruster supplier | Rafael  | Sentinel-1 A/B/C/D program selection.  | R24: printed p43; heritage printed p10 |
| S1-005 | Redundant branches | 2 branches | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-006 | Thrusters per branch | 7 thrusters | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-007 | Orbit-control thrusters per branch | 3 thrusters | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-008 | Attitude thrusters per branch | 4 thrusters | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-009 | RCT1 function | In-plane retrograde / orbit lowering  | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-010 | RCT2 function | In-plane prograde / orbit raising  | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-011 | RCT3 function | Out-of-plane control  | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-012 | OCT1 mounting side | +X spacecraft frame | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-013 | OCT2 mounting side | −X spacecraft frame | S1B documented architecture.  | S1360: p2; S1FD: p2 |
| S1-014 | OCT3 mounting side | −Y spacecraft frame | S1B documented architecture.  | S1360: p2; S1FD: p2 |

## S1: Thruster performance

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-015 | Propellant | N2H4  | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-016 | Type 1 inlet pressure | [9, 24.5] bar abs | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-017 | Type 2 inlet pressure | [5.5, 24.5] bar abs | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-018 | Type 1 BOL steady thrust | [1.3, 1.4] N | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-019 | Type 2 BOL steady thrust | [1, 1.1] N | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-020 | EOL steady thrust | >0.2 N | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-021 | Steady-state Isp at 22 bar | >214 s | Rafael 1 N product family; flight unit subtypes not established. Catalogued pressure point, not a flight map | R24: printed p46 / PDF p49 |
| S1-022 | Steady-state Isp at 5.5 bar | >205 s | Rafael 1 N product family; flight unit subtypes not established. Applicable inlet range must include this pressure | R24: printed p46 / PDF p49 |
| S1-023 | Minimum impulse bit | 0.008 N s | Rafael 1 N product family; flight unit subtypes not established. 5.5 bar; 0.02 s on / 1 s cycle | R24: printed p46 / PDF p49 |
| S1-024 | Nominal duty specification | 0.1 s on / 1 s cycle  | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-025 | Hot-pulse rise time | <200 ms | Rafael 1 N product family; flight unit subtypes not established. At nominal duty | R24: printed p46 / PDF p49 |
| S1-026 | Hot-pulse decay time | <300 ms | Rafael 1 N product family; flight unit subtypes not established. At nominal duty | R24: printed p46 / PDF p49 |
| S1-027 | Type 1 total delivered impulse | 100000 N s | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-028 | Type 2 total delivered impulse | 60000 N s | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-029 | Type 1 pulse qualification | 100000 pulses | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-030 | Type 2 pulse qualification | 58000 pulses | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-031 | Nozzle area ratio | 130 Ae/At | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-032 | Mass including 1000 mm leads | <=310 g | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |
| S1-033 | Combined storage and flight life | 15 years | Rafael 1 N product family; flight unit subtypes not established.  | R24: printed p46 / PDF p49 |

## S1: Valve and interfaces

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-034 | Valve type | Solenoid FCV, normally closed, two coils and two seats  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-035 | Closure mechanism | Springs at both seats; de-energise to close  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-036 | FCV construction | All-welded body, two springs and two plungers  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-037 | Inlet filter | 15 micrometre absolute | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-038 | FCV supply voltage | [23, 36] V DC | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-039 | FCV power | 9.2 W | Rafael 1 N product family. At 28 V DC | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-040 | Internal He leakage | <1e-5 standard cm3/s | Rafael 1 N product family. 5.5 and 24.5 bar as printed | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-041 | External He leakage | <1e-6 standard cm3/s | Rafael 1 N product family. 24.5 bar | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-042 | Threaded inlet | MS 33656-4  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-043 | Optional welded inlet tubing | ["1/4", "3/8"] inch | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-044 | Nozzle geometry | Bell nozzle welded to decomposition chamber  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-045 | Thermal barrier | Perforated Hayness-alloy tube between TCA and FCV  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |
| S1-046 | TCA-to-FCV seal | Metallic  | Rafael 1 N product family.  | R24: printed pp43,45–46 / PDF pp46,48–49 |

## S1: Thermal control

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-047 | Specified operating interface range | [5, 90] degC | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-048 | Non-operating range | [-10, 95] degC | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-049 | Catalyst-bed heater | Two resistance coils  | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-050 | Heater voltage | [23, 36] V DC | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-051 | Resistance per heater element | 257 ohm | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-052 | Recommended catalyst preheat temperature | 180 degC | Rafael 1 N product family; component limits have different boundaries. Heated catalyst bed, not the 5–90 degC interface envelope | R24: printed pp45–46 / PDF pp48–49 |
| S1-053 | Typical preheat time from +25 degC | 20 min | Rafael 1 N product family; component limits have different boundaries. Manufacturer example, not a spacecraft procedure | R24: printed pp45–46 / PDF pp48–49 |
| S1-054 | Heat-soak-back control | Thermal barrier limits conduction to FCV  | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |
| S1-055 | Additional FCV thermal fittings | External free areas allow heater/thermistor bonding  | Rafael 1 N product family; component limits have different boundaries.  | R24: printed pp45–46 / PDF pp48–49 |

## S1: Tank

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-056 | Tank model | PEPT-590GB  | Sentinel-1 A/B/C/D manufacturer heritage.  | R24: printed p9; p41 / PDF pp12,44 |
| S1-057 | Nominal internal volume | 204 L | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-058 | Nominal propellant capacity | 154 kg | PEPT-590GB product family. At stated 4:1 blow-down ratio | R24: printed pp40–41 / PDF pp43–44 |
| S1-059 | Nominal outer diameter | 587 mm | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-060 | Nominal overall length | 1170 mm | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-061 | Shape | Cylindrically extended hemispheres  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-062 | Dry tank mass | 17.5 kg | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-063 | Tank material | Ti-6Al-4V  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-064 | Diaphragm | EPN-40, silica free  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-065 | Compatible media | ["Hydrazine", "IPA", "distilled water", "GN2", "GHe"]  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-066 | Compatible pressurants | ["GN2", "GHe"]  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-067 | Liquid volumetric expulsion | >99 percent | PEPT-590GB product family. Does not include motor-pressure cutoff or unusable propellant from other causes | R24: printed pp40–41 / PDF pp43–44 |
| S1-068 | Operating pressure | [5.5, 27.4] bar | PEPT-590GB product family. Pressure reference not explicit in this table | R24: printed pp40–41 / PDF pp43–44 |
| S1-069 | Proof pressure | 41.1 bar | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-070 | Specified burst pressure | 54.8 bar | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-071 | Reported actual rupture pressure | 72 bar | PEPT-590GB product family. Single catalogue statement, not a flight allowable | R24: printed pp40–41 / PDF pp43–44 |
| S1-072 | Internal He leakage | <1e-2 standard cm3/s | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-073 | External He leakage | <1e-6 standard cm3/s | PEPT-590GB product family. At 27.4 bar | R24: printed pp40–41 / PDF pp43–44 |
| S1-074 | Qualification operating temperature | [4, 50] degC | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-075 | Qualification non-operating temperature | [-20, 60] degC | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-076 | Storage life | 7 years | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-077 | Flight life | 15 years | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-078 | Mounting | Polar  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |
| S1-079 | Gas and liquid ports | Weldable 1/4-inch titanium tubes  | PEPT-590GB product family.  | R24: printed pp40–41 / PDF pp43–44 |

## S1: Feed system

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1-182 | Feed mode | Pressure-fed blow-down  | Sentinel-1 family.  | S1FD: p2; R24: printed p41 |

## S1C: Flight loading

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1C-080 | Loaded propellant | 154 kg | Sentinel-1C, November 2024 campaign.  | F1C: fuel-loading paragraph; R24: printed p41 nominal product capacity, not the loading event |
| S1C-081 | Propellant identity | Hydrazine  | Sentinel-1C.  | R24: printed p43 A/B/C/D selection; S1SP: hydrazine propulsion description |
| S1C-082 | Launch mass used for calculations | 2186 kg | S1C VV25 launch-kit value.  | KIT1C: p3; FACT1: generic family mass 2185 kg |

## S1C: Open flight data

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1C-159 | C/D as-built thruster subtype distribution | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-160 | C/D thruster cant angles and moment arms | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-161 | Flight thrust/Isp/mass-flow map versus P/T/duty | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-162 | Flight residual P/T and pressurant mass | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-163 | Motor-by-motor impulse/pulse/cold-start ledger | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-164 | Current C/D EOL propellant allocation | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-165 | Complete propulsion dry mass | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-166 | Flight valve/heater switching power profile | Not publicly established  | S1C.  | No public acceptance record recovered |
| S1C-181 | Actual loaded pressurant species and initial pressure | Not publicly established  | S1C.  | No public acceptance record recovered |

## S1B: Operational constraints

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1B-083 | Affected in-plane burn duration | 35 s | Specific 2016 acquisition conditions, not all vehicles or epochs. Reduced from 300 s because of plume constraints | S1FD: PDF pp2–3 |
| S1B-084 | Reported effective-thrust loss | approximately 30 percent | Specific 2016 acquisition conditions, not all vehicles or epochs. Plume interaction; not a universal nozzle-efficiency factor | S1FD: PDF pp2–3 |
| S1B-085 | OCT3 large in-plane burn duration | 600 s | Specific 2016 acquisition conditions, not all vehicles or epochs. Special 90-degree slew procedure; not a universal out-of-plane limit | S1FD: PDF pp2–3 |
| S1B-086 | Recovery after OCT3 slew/burn/return sequence | 6 h | Specific 2016 acquisition conditions, not all vehicles or epochs. Nominal pointing restores thermal state and battery charge | S1FD: PDF pp2–3 |

## S1B: End-of-life operation

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S1B-087 | Simulated RCT4A / 5A / 6A / 7A forces | [0.57, 0.56, 0.47, 0.58] N | Specific 2023–2024 disposal campaign. Simulator inputs representing a degraded flight condition | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-088 | Residual-fuel uncertainty | approximately 5 kg | Specific 2023–2024 disposal campaign. Campaign estimate, not a C/D requirement | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-089 | Stop-Burn wheel speed threshold | 250 rad/s magnitude | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-090 | Stop-Burn attitude-error threshold | 0.01 rad magnitude | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-091 | Stop-Burn rate-error threshold | 0.001 rad/s magnitude | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-092 | Double-thruster failure validation | Not performed in the reported study  | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-093 | Final reported orbit | [580, 610] km perigee/apogee | Specific 2023–2024 disposal campaign. Approximate | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-094 | Remaining fuel at campaign stop | approximately 70 kg | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-095 | Predicted remaining orbital life | approximately 24 years | Specific 2023–2024 disposal campaign. Model prediction, not observed re-entry | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-096 | Dedicated in-orbit fuel vent | Not present  | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-097 | Electrical passivation | 11–12 September 2024 UTC dates | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |
| S1B-098 | Complete fuel/pressure passivation | Not achieved  | Specific 2023–2024 disposal campaign.  | S1360: pp3,7,10,12; S1226: pp16–21 |

## S2: Architecture

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-099 | Propulsion type | Catalytic hydrazine monopropellant  | Sentinel-2 family.  | S2SP: printed pp15–16 / PDF pp23–24; S214: pp5,12; F2C: loading paragraph |
| S2-100 | Installed thruster count | 8 thrusters | S2 A/B design; C/D detailed acceptance configuration not public.  | S214: pp5,12: 2 cold-redundant branches of four; S209: PDF p2: two times four thrusters |
| S2-101 | Specific thruster model | CHT1N-6  | 2014 A/B AOCS equipment selection.  | S214: p9 Table 5 |
| S2-102 | Supplier/product-family association | EADS ST / later ArianeGroup 1 N family  | A/B design and manufacturer heritage.  | S214: p9 Table 5; AGPDF: p6: S2A and S2B heritage |
| S2-103 | Branch redundancy | Two cold-redundant branches  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-104 | Thrusters per active branch | 4 thrusters | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-105 | Thruster installation | Tilted relative to spacecraft x-axis  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-106 | Function | Orbit manoeuvres and three-axis torque  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-107 | Integration | Self-standing tank/pipe/thruster module through launcher interface ring  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-108 | OCM sequence | ["Slew", "Stabilise", "Delta-V", "Back slew"]  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-109 | OCM attitude actuation | Thruster off-modulation  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-110 | Normal attitude actuation | Reaction wheels with magnetic momentum unloading  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |
| S2-111 | Safe-mode actuation | Thruster pulses supported by magnetic torquers  | 2014 A/B design.  | S214: pp5,8,12; S2SP: printed p16 / PDF p24 |

## S2: Thruster performance

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-112 | Nominal thrust | 1 N | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-113 | Thrust range | [0.32, 1.1] N | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-114 | Nominal Isp | 220 s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-115 | Isp range | [200, 223] s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row |
| S2-116 | Nominal mass flow | 0.44 g/s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-117 | Mass-flow range | [0.142, 0.447] g/s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-118 | Inlet pressure range | [5.5, 22] bar | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map. Pressure reference is not explicitly labelled in the cited table | AGPDF: p6; AGWEB: characteristics table |
| S2-119 | Minimum impulse bit | [0.01, 0.043] N s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map. PDF unit; web page incorrectly prints N | AGPDF: p6; AGWEB: characteristics table |
| S2-120 | Nozzle area ratio | 80 Ae/At | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-121 | Thruster mass with valves | 290 g | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row |
| S2-122 | Propellant grade | High-purity hydrazine  | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-123 | Total impulse qualification | 135000 N s | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table; NASA26: Table 4-2, ArianeGroup 1 N row |
| S2-124 | Cycle-life qualification | 59000 cycles | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-125 | Propellant throughput qualification | 67 kg | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-126 | Single-burn life qualification | 12 h | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-127 | Accumulated burn-life qualification | 50 h | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |
| S2-128 | Qualified cold starts below 10 degC | 10 starts | ArianeGroup general 1 N catalogue; A/B heritage, not C/D acceptance map.  | AGPDF: p6; AGWEB: characteristics table |

## S2: Valve and thermal design

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-129 | Flow-control valve | Two identical monostable normally-closed valves in series, single housing  | ArianeGroup general 1 N product family.  | AGPDF: p5; AGWEB: general description |
| S2-130 | Catalyst heater | Internally redundant catalyst-bed heater  | ArianeGroup general 1 N product family.  | AGPDF: p5; AGWEB: general description |
| S2-131 | Thermal protection | Insulation and heat barrier to protect valve/structure  | ArianeGroup general 1 N product family.  | AGPDF: p5; AGWEB: general description |
| S2-132 | Cold-start capability | Multiple cold starts qualified  | ArianeGroup general 1 N product family.  | AGPDF: p5; AGWEB: general description |
| S2-133 | Valve and chamber materials | Selected for hydrazine compatibility; exact alloys not stated  | ArianeGroup general 1 N product family.  | AGPDF: p5; AGWEB: general description |

## S2: Dimensions

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-134 | Valve-side reference length | 97 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-135 | Mounting-plane to hot-end reference length | 75 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-136 | Heater assembly reference length | 42 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-137 | Hot-end drawing reference diameter | 30 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-138 | Drawing reference diameter | 41 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-139 | Mounting hole diameter | 4.3 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |
| S2-140 | Mounting feature reference dimension | 23 mm | Manufacturer general 1 N drawing; not a C/D flight ICD.  | AGWEB: Elevation / Plan / Cut-out drawing links |

## S2: End-of-life plan

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-145 | Initial circular altitude | 786 km | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-146 | Disposal circular altitude | 653 km | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-147 | Planned delta-V | 73 m/s | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-148 | Planned propellant | 47 kg | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-149 | Planned manoeuvre count | 116 burns | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-150 | Planned burn duration | 600 s | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-151 | Planned frequency | One burn per orbit  | 2014 family design; not an executed C/D campaign.  | S214: p24 |
| S2-152 | Post-passivation re-entry target | 25 years | 2014 family design; not an executed C/D campaign. Historical design requirement | S214: p24 |

## S2: Feed system

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2-183 | Feed mode | Helium-pressurised blow-down  | A/B design; C helium loading confirmed.  | S214: p12; F2C: loading paragraph |

## S2C: Flight loading

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2C-141 | Loaded hydrazine | 133 kg | 16 August 2024 flight campaign.  | F2C: loading paragraph |
| S2C-142 | Pressurant | Helium  | 16 August 2024 flight campaign.  | F2C: loading paragraph |
| S2C-143 | Tank loading pressure | 21 bar | 16 August 2024 flight campaign. Absolute/gauge reference not stated; tank pressure, not motor inlet | F2C: loading paragraph |
| S2C-144 | Launch mass used for calculations | 1143 kg | VV24 vehicle-specific launch-kit value.  | KIT2C: p3; FACT2: generic 1140 kg including 123 kg fuel |

## S2C: Published manoeuvres

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2C-153 | Large in-plane acquisition total | approximately 6.33 m/s | 2024 acquisition/tandem summary; not a complete telemetry ledger.  | S2617: pp15–16; Table 3 |
| S2C-154 | 29 October fine acquisition | 0.012 m/s | 2024 acquisition/tandem summary; not a complete telemetry ledger.  | S2617: pp15–16; Table 3 |
| S2C-155 | Initial out-of-plane total | approximately 0.57 m/s | 2024 acquisition/tandem summary; not a complete telemetry ledger. Source also uses planning language | S2617: pp15–16; Table 3 |
| S2C-156 | 6 November calibrated OOP/IP components | [0.9344, 0.0537] m/s | 2024 acquisition/tandem summary; not a complete telemetry ledger.  | S2617: pp15–16; Table 3 |
| S2C-157 | 6 November duration | 314.1 s | 2024 acquisition/tandem summary; not a complete telemetry ledger.  | S2617: pp15–16; Table 3 |
| S2C-158 | 23 November CAM components in separate burns | [-0.0207, 0.0206] m/s | 2024 acquisition/tandem summary; not a complete telemetry ledger. Fuel cost sums magnitudes, not signed cancellation | S2617: pp15–16; Table 3 |

## S2C: Open flight data

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2C-167 | C/D as-built thruster subtype distribution | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-168 | C/D thruster cant angles and moment arms | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-169 | Flight thrust/Isp/mass-flow map versus P/T/duty | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-170 | Flight residual P/T and pressurant mass | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-171 | Motor-by-motor impulse/pulse/cold-start ledger | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-172 | Current C/D EOL propellant allocation | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-173 | Complete propulsion dry mass | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-174 | Flight valve/heater switching power profile | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-175 | Tank supplier/model | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-176 | Tank volume and dry mass | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-177 | Tank material/diaphragm specification | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-178 | Tank proof/burst pressures | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-179 | Flight piping/filter pressure losses | Not publicly established  | S2C.  | No public acceptance record recovered |
| S2C-180 | Valve power and heater resistance | Not publicly established  | S2C.  | No public acceptance record recovered |

## S2C: Feed system

| ID | Parameter | Value | Applicability / condition | Evidence |
|---|---|---|---|---|
| S2C-184 | Latch-valve isolation | Closed and verified during launch campaign  | 16 August 2024 campaign.  | F2C: post-loading spacecraft checks |

## Source registry

- **R24** — [Space Propulsion Technical Catalogue, 2024](https://www.rafael.co.il/wp-content/uploads/2024/05/SPACE-Propulsion-Technical-doc_compressed.pdf). Rafael. official publisher. Local file: `sources/Rafael_propulsion_2024.pdf`. SHA-256: `f621c47cf06ed8f7064b7ec41217d950d10e562928169f1a8d80c3a44123ed1b`.
- **AGPDF** — [Hydrazine Thrusters brochure](https://www.space-propulsion.com/brochures/hydrazine-thrusters/hydrazine-thrusters.pdf). ArianeGroup. official publisher. Local file: `sources/ArianeGroup_hydrazine_thrusters.pdf`. SHA-256: `90068e721a853845c46f559ddb79cdc7c6276961bafdcc16d7152cd5d52d969b`.
- **AGWEB** — [1 N Monopropellant Hydrazine Thruster](https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html). ArianeGroup. official publisher. Local file: `sources/ArianeGroup_1N_product.html`. SHA-256: `20a5e656b2ad7fedad82bef52ecffadc273f6ed8874bcc625ee3c6eaf3f77a49`.
- **S1360** — [Stop OCM Burn FDIR, SpaceOps 2025 #360](https://publications.spaceops.org/2025/download_by_id.php?id=0360). ESA / Serco mission operations. official publisher. Local file: `sources/SpaceOps_2025_0360_S1B_StopBurn.pdf`. SHA-256: `153df80ffa07cf69e91b855cd1a34883bfd8a040570a832334350a420285ab89`.
- **S1226** — [Sentinel-1B End-Of-Life Operations, SpaceOps 2025 #226](https://publications.spaceops.org/2025/download_by_id.php?id=0226). ESA mission operations. official publisher. Local file: `sources/SpaceOps_2025_0226_S1B_EndOfLife.pdf`. SHA-256: `c03b7153b8c98f9f3d289e6edf6414f5f26499a0fe6a5871deb80cf6facb6afd`.
- **S2617** — [Optimising Tandem Orbit Control, SpaceOps 2025 #617](https://publications.spaceops.org/2025/download_by_id.php?id=0617). ESA / CS Group mission operations. official publisher. Local file: `sources/SpaceOps_2025_0617_S2_Tandem.pdf`. SHA-256: `53b1f906bc44082904d5192851a6631f526d4531b127d1dd36c1ccbf3e7f51de`.
- **S1FD** — [Sentinel-1B LEOP and Reference Orbit Acquisition, 2017](https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-029__ISSFD-2017-029.pdf). ESA / ESOC flight dynamics. official publisher. Local file: `sources/ISSFD_2017_S1_OrbitAcquisition.pdf`. SHA-256: `f57ded32775361b2154c6a1a7387a72ee383a11f67044a71540091cb63d14370`.
- **S1SP** — [SP-1322/1, 2012](https://sentinels.copernicus.eu/documents/247904/349449/S1_SP-1322_1.pdf). ESA. official publisher. Local file: `sources/Sentinel1_ESA_SP1322_2012.pdf`. SHA-256: `b31d601a798a8cbf37614b8d5e78cb96499acaae9648c202e3b47ad61f9c5b8e`.
- **S2SP** — [SP-1322/2, 2012](https://sentinels.copernicus.eu/documents/247904/349490/S2_SP-1322_2.pdf). ESA. official publisher. Local file: `sources/Sentinel2_ESA_SP1322_2012.pdf`. SHA-256: `517366ed447a21e7e2f00188abc9b8e4c4228e133f6c9f07c473f1d894288fb4`.
- **S209** — [Sentinel-2 AOCS Model-Based Design, DGLR 2009](https://www.fzt.haw-hamburg.de/pers/Scholz/dglr/dlrk2009_mitReview/Papers/121179.pdf). EADS Astrium. primary authors; university-hosted copy. Local file: `sources/Sentinel2_AOCS_ModelBasedDesign_DGLR2009.pdf`. SHA-256: `a0e68403ad4894c3b6528dde54a7aa9d4743820f15b56528c9024030d41923e4`.
- **S214** — [Sentinel-2 AOCS: Challenges and Solutions, GNC 2014](https://www.researchgate.net/publication/271584776_THE_SENTINEL-2_SATELLITE_ATTITUDE_CONTROL_SYSTEM_-_CHALLENGES_AND_SOLUTIONS). Airbus Defence and Space. primary author copy on ResearchGate; full-text read; automated download unavailable. Local file: `None`. SHA-256: `not downloaded`.
- **F1C** — [Sentinel-1C fuelled ahead of liftoff, 18 November 2024](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Sentinel-1C_fuelled_ahead_of_liftoff). ESA. official publisher. Local file: `sources/ESA_S1C_Fuel.html`. SHA-256: `261a3a955f9104541bccae8d900036207964a8f532a5cd1cccb5fecf8de7dcea`.
- **F2C** — [Sentinel-2C fully loaded, 19 August 2024](https://blogs.esa.int/eolaunches/2024/08/19/sentinel-2c-fully-loaded/). ESA. official publisher. Local file: `sources/ESA_S2C_Fuel.html`. SHA-256: `bd8e7d9e607093ce45f458f90f9f1935770615e1937b4802261f820486250a9f`.
- **FACT1** — [Sentinel-1 Facts and figures](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures). ESA. official publisher. Local file: `sources/ESA_S1_Facts.html`. SHA-256: `6faf70a76f55a1f53aee76757ccc267c7163ee2033062d306a318d3fd57ff81f`.
- **FACT2** — [Sentinel-2 Facts and figures](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures). ESA. official publisher. Local file: `sources/ESA_S2_Facts.html`. SHA-256: `1b32fe6901d1b90875b09e89f2f829addb70c3bb95b0e1ec724dda26663fa3ea`.
- **KIT1C** — [VV25 launch kit](https://sky-brokers.com/wp-content/uploads/2025/07/Arianespace-Launchkit-VV25-ESA-Copernicus-Sentinel-1c.pdf). Arianespace. primary-authored third-party mirror; official-host copy not recovered. Local file: `sources/S1C_VV25_LaunchKit_mirror.pdf`. SHA-256: `5f387399c33840c042a6f353b7b4870ca4059305f7f153f8aea77c0afe5a2367`.
- **KIT2C** — [VV24 launch kit](https://newsroom.arianespace.com/download/?n=b694039895faf86c0a7a5bcab504d9a1.pdf&id=39971). Arianespace. official publisher. Local file: `sources/S2C_VV24_Official_LaunchKit.pdf`. SHA-256: `fd96c6f27b335d95e122fe48b2be51a02e70679939cc76c43d4c330f254dedd1`.
- **NASA26** — [2026 State of the Art, In-Space Propulsion](https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/). NASA. official compilation; table values supplied by manufacturers. Local file: `sources/NASA_propulsion_2026.html`. SHA-256: `751ab86c4d9b85a464eb680d7a98eef7db1f6777aae3198bc28fd8b84fe7d89f`.
- **N2H4** — [NASA/TM-2013-216560, Recommended Figures of Merit for Green Monopropellants](https://ntrs.nasa.gov/api/citations/20140009582/downloads/20140009582.pdf). NASA. official publisher. Local file: `sources/NASA_TM_2013_216560.pdf`. SHA-256: `3758f8500534c4bd82b88c44cc7180df3b15e75d99a2e99b636494162434d6b8`.
- **GUIDE** — [Guidelines for the computation of Delta-V and propellant budget, Issue 1 Rev 0, 11 November 2019](https://eopro.esa.int/wp-content/uploads/2020/10/Guidelines-for-DeltaV-and-propellant-budget-computation.pdf). ESA EOP. official publisher. Local file: `sources/ESA_DeltaV_Propellant_Guidelines_2019.pdf`. SHA-256: `b1ea2d4cf311f825150b6590424514b55c37aae5d6a4574ba7612c7642cce7c2`.

## Calculations

The [English calculation appendix](CALCULATIONS_ENGLISH.md) contains the derivations. Machine-readable results are in SERI_SONUCLARI.json and DERIN_ITKI_SONUCLARI.json. Scenario Isp, cant angle, line loss and annual demand are not catalogue specifications. The 9-slide presentation labels its ideal calculations. ESA guideline GUIDE §§4–5 provides a separate official method check; its Annex A example is not Sentinel flight data.
