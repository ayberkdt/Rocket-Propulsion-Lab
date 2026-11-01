# Sentinel-1 / Sentinel-2 propulsion calculations

**V6 · English · 1 November 2026.** This appendix supports the nine-slide presentation.
The [propulsion catalogue](PROPULSION_CATALOGUE.md) contains the source locations,
vehicle applicability and verification status of the hardware facts. Calculations
below are analytical estimates and labelled scenarios, not measured fuel consumption.

## 1. Matched inputs and units

| Input | Sentinel-1C | Sentinel-2C | Evidence / use |
|---|---:|---:|---|
| Launch mass, m₀ | 2,186 kg | 1,143 kg | Vehicle-specific VV25 / VV24 kits, p3 |
| Loaded hydrazine | 154 kg | 133 kg | ESA flight-campaign records, F1C / F2C |
| Mass after subtracting hydrazine | 2,032 kg | 1,010 kg | Arithmetic remainder, not certified dry mass |
| Reference altitude | 693 km | 786 km | ESA family mission design |
| Illustrative constant Isp | 205 s | 220 s | Chosen calculation points, not C flight maps |
| Isp sensitivity | 205 / 214 / 220 s | 200 / 220 / 223 s | Scenarios informed by product data |

Use g₀ = 9.80665 m/s², Earth μ = 3.986004418 × 10¹⁴ m³/s² and Earth equatorial
radius R = 6,378,137 m. The last two constants are the model convention in the
reproducible scripts. Heights and radii must not be interchanged. Pressure units in
the catalogue retain the source's absolute/gauge qualification. A pressure with an
unspecified reference must not enter a gas-law calculation as known absolute pressure.

The general ESA fact sheets report 2,185 kg and 1,140 kg. The latter includes 123 kg
fuel, while the S2C campaign reports 133 kg. These records refer to different
vehicle/design summaries. This study does not mix them into a fictitious mass budget.
The S1C kit was recovered through a third-party mirror; the S2C kit is archived from
the official publisher. Exact independent two-authority confirmation of either
C launch mass was not recovered.

## 2. Scalar thrust, useful thrust and fuel

For a constant operating point:

```text
c = g₀ Isp                                      [m/s]
ṁ = ΣFᵢ / c                                     [kg/s]
Faxis = Σ(Fᵢ cos θᵢ)                             [N]
η = Faxis / ΣFᵢ                                 [dimensionless]
Δv = η c ln(m₀ / mf)                            [m/s]
mp = m₀ [1 − exp(−Δv / (η c))]                  [kg]
t = mp / ṁ                                     [s]
```

These expressions assume common Isp, constant thrust directions and a fixed useful
axis. Attitude off-modulation, torque generation and plume losses require additional
accounting. For arbitrary flight histories, integrate thrust vectors and each
thruster's mass flow separately. Installed thruster count is not the number firing
at once. A redundant branch is not an additional continuously available axial force.

| Ideal Δv = 1 m/s, η = 1 | S1C at 205 s | S2C at 220 s |
|---|---:|---:|
| Consumed propellant | 1.087095 kg | 0.529666 kg |
| Flow for one assumed 1 N | 0.497423 g/s | 0.463507 g/s |
| Duration for one axial 1 N | 2,185.456 s | 1,142.735 s |
| Duration for four axial 1 N | Not assigned to the flight layout | 285.684 s, illustrative zero-cant case |

Four ideal parallel thrusters shorten the burn by four. They do not reduce the
ideal fuel at unchanged Isp and η. At a common cant angle θ, η = cos θ: 20° gives
0.939693 and 30° gives 0.866025. Flight cant angles were not established, so these
are sensitivities. Do not infer a flight angle from a calibrated mean force alone.

The ArianeGroup nominal labels 1 N, 220 s and 0.44 g/s are not an exact matched
point: 1/(220g₀) = 0.463507 g/s, whereas 0.44 g/s × 220g₀ = 0.949284 N.
The catalogue's range endpoints must not be paired arbitrarily to construct a
flight performance map. This is a consistency observation, not a claim of faulty hardware.

## 3. Ideal loaded-fuel capability

Set mf = m₀ − loaded hydrazine and η = 1. This idealises all loaded fuel as usable;
it includes no disposal reserve, inaccessible residual or pressure cutoff.

| Vehicle | Assumed Isp | Ideal capability |
|---|---:|---:|
| S1C | 205 s | 146.862785 m/s |
| S1C | 214 s | 153.310420 m/s |
| S1C | 220 s | 157.608843 m/s |
| S2C | 200 s | 242.628395 m/s |
| S2C | 220 s | 266.891234 m/s |
| S2C | 223 s | 270.530660 m/s |

These are neither mission requirements nor certified available Δv. Higher constant
Isp produces a larger ideal estimate; actual blow-down Isp varies with conditions.

## 4. Sentinel-2C published manoeuvre subset

ESA/CS Group SpaceOps 2025 paper #617, pp15–16, describes acquisition and tandem
control. Its summary includes approximate/planning language for some phases.
Table 3 contains calibrated small manoeuvres; it is not a complete lifetime ledger.

| Event | Fuel-cost Δv basis, m/s | Published duration |
|---|---:|---:|
| Four large in-plane manoeuvres, summary | approximately 6.33 | Not specified in the summary |
| 29 October fine acquisition | 0.012 | Not specified here |
| Initial out-of-plane summary | approximately 0.57 | Source also uses planning language |
| 31 October | 0.0146 | 5.3 s |
| 6 November | √(0.9344² + 0.0537²) = 0.935941798 | 314.1 s |
| 8 November | 0.0025 | 0.9 s |
| 23 November avoidance and return | 0.0207 + 0.0206 = 0.0413 | 7.3 + 7.3 s |
| 27 November | 0.0309 | 11 s |
| 29 November | 0.0029 | 1 s |
| 18 December | 0.0283 | 10 s |
| Sum of the listed magnitudes | approximately 7.97 | Not a single burn |

The arithmetic sum before presentation rounding is 7.968441798 m/s. Source precision does not justify
reporting that many digits as flight accuracy. Avoidance and return have opposite
signed components but both consume fuel. Their fuel cost is not the signed sum.

Update the mass after each event:

```text
m[j+1] = m[j] exp(−Δv[j] / (g₀ Isp η))
fuel[j] = m[j] − m[j+1]
```

At m₀ = 1,143 kg, Isp = 220 s and η = 1, the subset consumes **about 4.21 kg**;
the model remainder is **128.79 kg**. These are not onboard measurements. The
27 combinations of mass 1,140 / 1,143 / 1,146 kg, Isp 200 / 220 / 223 s and cant
0 / 20 / 30° give **4.1463–5.3636 kg**. This bounded scenario envelope is not a
statistical confidence interval. The mass alternatives are assumptions, not three
documented S2C launch masses.

For 6 November, the model mass before firing is 1,139.336 kg. The rocket equation
and 314.1 s duration imply mean useful force of about **3.394 N**. This cannot
uniquely identify cant, instantaneous pressure, valve duty or individual thrust.

## 5. Disposal geometry and mission allocation

For a coplanar two-impulse transfer between circular radii r₁ and r₂:

```text
a = (r₁ + r₂)/2
v₁ = √(μ/r₁);  v₂ = √(μ/r₂)
vt1 = √[μ(2/r₁ − 1/a)];  vt2 = √[μ(2/r₂ − 1/a)]
Δv = |vt1 − v₁| + |v₂ − vt2|
```

The ideal 786→653 km circular comparison is **70.2157 m/s**. Airbus's 2014
Sentinel-2 design plan gives **73 m/s, 47 kg and 116 burns of 600 s**. The difference
can include design assumptions and operational margins; it is not evidence that
the executed C campaign consumed 47 kg. A current C/D disposal allocation is unavailable.

For S1, 693→530 km circular gives **88.0576 m/s** under this model. The 530 km
value is an older target, not S1B's achieved final orbit. A geometric comparison
from 693 km circular to **580 × 610 km** elliptical requires about **52.5822 m/s**.
It is not a sum reconstructed from S1B's flight burns. The 2025 ESA EOL paper
reports approximately 70 kg remaining and a predicted roughly 24-year orbital
lifetime, with partial passivation. Both are campaign-specific estimates.

The complete spacecraft budget must separately allocate injection correction,
ground-track/eccentricity control, plane corrections, collision avoidance and return,
disposal, uncertainty and operational contingency. The published subset cannot
close these categories over the entire mission.

## 6. Blow-down sensitivity and pressure-limited residuals

Combine the **S1 family product tank** volume 204 L with the S1C 154 kg fill only
as a sensitivity case. Density ρ = 1,004 kg/m³ at 300 K is from NASA/TM-2013-216560,
Table 1. A diaphragm separates pressurant from liquid; effective gas volume and
thermal behaviour require a flight tank model.

```text
Vg0 = Vtank − mp0/ρ = 50.61355 L
P Vg^n = P0 Vg0^n
Pcutoff = motor minimum inlet + assumed feed loss = 5.5 + 0.5 = 6 bar abs
Vg,cutoff = min[Vtank, Vg0 (P0/Pcutoff)^(1/n)]
mp,residual = max[0, ρ(Vtank − Vg,cutoff)]
```

| Assumed initial pressure, bar abs | n | Model residual, kg |
|---:|---:|---:|
| 22 | 1 | 18.49 |
| 22 | 1.2 | 54.77 |
| 22 | 5/3 | 94.01 |
| 24.5 | 1 | 0.00 |
| 24.5 | 1.2 | 40.69 |
| 24.5 | 5/3 | 86.62 |

Initial pressure, n, feed loss, pressurant behaviour and chosen Type 2 cutoff are
assumed. These are not predicted S1C residuals. Type 1 has a different minimum inlet.
S2C's reported 21 bar has no stated absolute/gauge basis and must not be inserted
as a known absolute pressure. S2 tank volume is not inferred from the S1 tank.
The S1 tank's 27.4 bar operating limit is not the motor's 24.5 bar inlet limit;
assuming only 0.5 bar loss at 27.4 bar would exceed that motor inlet limit.
Volumetric expulsion >99% does not guarantee >99% usable fuel at the required inlet pressure.

## 7. Impulsive and finite-burn models

The finite-burn example uses m₀ = 2,100 kg, altitude 693 km, Isp = 205 s and a
0.7 m/s scalar retrograde manoeuvre. Thrust follows the instantaneous anti-velocity
direction. RK4 integrates position, velocity, mass, thrust work and angular impulse.
J₂, drag, plume interaction, attitude dynamics and pressure variation are excluded.

| Assumed force | Duration | Final perigee | Final apogee |
|---:|---:|---:|---:|
| 1 N | 1,469.744 s | 690.49328 km | 692.87021 km |
| 2.18 N | 674.195 s | 690.39150 km | 692.97203 km |
| 4 N | 367.436 s | 690.37188 km | 692.99165 km |
| Instantaneous impulse | Zero | 690.36354 km | 693.00000 km |

All finite cases consume 0.731084 kg under the constant-Isp scalar specification.
At 2.18 N the burn occupies about 11.4% of an orbit; its apsides differ by about
28 m from the instantaneous case. That is an orbital-shape difference, not a
28% fuel penalty. The half-step calculation and work/angular-momentum/mass
residual checks are retained in DERIN_ITKI_SONUCLARI.json and the numerical tests.

## 8. Impulse, pulse life, thermal duty and power

```text
I_total = ∫ΣFᵢ dt ≈ mp g₀ Isp                 [N·s]
I_unit = I_total / N                         [equal-sharing screen only]
Δv_min ≈ I_bit / m                           [small pulse approximation]
```

Burning 154 kg at 205 s corresponds to 309,596 N·s scalar impulse. Dividing that
equally among four S1 Type 2 units gives 77,399 N·s per unit, above the 60,000 N·s
product figure. This does not prove a flight violation: the actual subtype,
allocation, operating point and impulse history are unknown. It demonstrates why
a whole-tank fuel number cannot replace motor-by-motor life accounting.

For S2, 133 kg at 220 s corresponds to 286,943 N·s, or 71,736 N·s each under
four-unit equal sharing. This is below the 135,000 N·s catalogue figure, but pulse
count, throughput, cold starts and firing duration still require separate checks.
The historical 116 × 600 s disposal plan totals 19.333 hours, not a single
19.333-hour firing. Catalogue single-burn and accumulated limits must stay separate.

At launch mass, the S1 0.008 N·s bit corresponds to about 3.66 µm/s; the S2
0.010–0.043 N·s range gives about 8.75–37.62 µm/s. These ideal velocity quanta
are not orbit-determination accuracy or commanded repeatability guarantees.

The Rafael valve is specified at 9.2 W at 28 V. A 257 Ω heater element would draw
28²/257 = 3.05 W at 28 V under an ideal resistive calculation. Element count,
switching logic, warm-up temperature and duty determine system power. Recommended
180°C catalyst preheat and +5 to +90°C interface operation refer to different
thermal locations. No S2 flight heater wattage or complete thermal budget is invented.

The 2017 S1 acquisition paper limits plume-affected in-plane firings to 35 s.
Its special 90° slew procedure aligns OCT3 for large in-plane burns up to 600 s,
then requires six hours in nominal pointing for thermal/battery recovery. This is
not a universal 600 s out-of-plane thruster limit.

## 9. Drag, inclination and extra operating life

```text
a_drag = 0.5 ρ_atm Cd A v² / m
Δv_drag,year ≈ a_drag × (365.25 × 86400)
Δv_plane = 2v sin(Δi/2)
```

Drag sensitivities in the numerical file use atmospheric density 10⁻¹⁴ / 10⁻¹³ /
10⁻¹² kg/m³, Cd = 2.2 and area 10 / 20 m². These are assumptions. Solar activity,
attitude, density history and geometry are needed for a mission forecast. A small
inclination change is not an altitude change; use the appropriate vector geometry.

If current mass is m, genuinely accessible remaining fuel is mr, withheld fuel
for disposal/uncertainty is mh and annual useful demand is d:

```text
mavailable = max(0, mr − mh)
Δvavailable = g₀ Isp η ln[m / (m − mavailable)]
extra years ≈ Δvavailable / d
```

An explicitly assumed example m = 2,100 kg, mr = 70 kg, mh = 60 kg, Isp = 205 s,
η = 1 and d = 1 m/s/year gives **9.596 ideal years**. It is not an S1B/C/D lifetime
prediction. Changing pressure accessibility, reserve or annual demand changes it
directly. Public loading records alone cannot establish additional operating life.

## 10. Reproduction and evidence boundary

The method was also checked against ESA's [2019 EOP budget guideline](https://eopro.esa.int/wp-content/uploads/2020/10/Guidelines-for-DeltaV-and-propellant-budget-computation.pdf),
§§4–5. It builds fuel requirements backward from disposal, with separate mission
phases, residuals, gauging uncertainty and manoeuvre-dependent efficiency.
That edition uses 1% of maximum tank load for residuals and 2% for gauging
uncertainty; these are guideline allocations, not measured Sentinel residuals.
Its Annex A uses an example mission and must not supply Sentinel flight inputs.
Mission-specific requirements take precedence. No arbitrary blanket percentage
is added to the displayed ideal calculations. A design budget would instead use
phase-specific Isp/efficiency and explicitly justified allocations.

Run from the repository root with Python 3.11 or later:

```text
python scripts/sentinel_series_budget.py --output docs/missions/sentinel/SERI_SONUCLARI.json
python scripts/sentinel_propulsion_deep_dive.py --output docs/missions/sentinel/DERIN_ITKI_SONUCLARI.json
python -m unittest discover -s tests -p "test_sentinel*.py"
```

The scripts use the standard library. Tests compare independent high-precision
equations, conservation checks, event accounting and finite-step convergence;
they do not test the spacecraft hardware. Qualification tables establish only
the stated product envelope. Private C/D acceptance maps, complete fuel ledgers,
current EOL allocations and numerical preflight TRL certificates were not recovered.
Flight heritage is evidence of flown hardware, not qualification of our future design.
