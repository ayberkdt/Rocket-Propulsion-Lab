# Propellant vaporization (Priem-Heidmann)

`rocket_propulsion.thermochemistry.vaporization` estimates how much of each
propellant vaporizes before the throat, and the c* efficiency that follows.
Priem computed drop histories (Ranz-Marshall heat and mass transfer, drag
`C_D = 27 Re^-0.84`, gas velocity proportional to the mass already
vaporized) for log-normal sprays and collapsed them onto one curve of percent
vaporized against a "generalized length".

Tests: `tests/test_vaporization.py` (9 tests). HTTP route:
`POST /api/v1/propellant-vaporization` (LOX/CH₄ and LOX/LH₂, with T_c, H_v and
M from the reference equations of state); workbench: **Yanma → Ürünler ve
verim → Buharlaşma**.

## Public API

| Function | What it gives | Source |
|---|---|---|
| `generalized_vaporization_length(...)` | `L_gen` of one propellant (SI in and out) | NASA CR-168066 eq. (1) and program statement `CVO`, after Priem & Heidmann, NASA TR R-67 |
| `percent_vaporized(L_gen)` | percent of mass vaporized, the 35-point TR R-67 average curve, linear interpolation | CR-168066 `DATA LG/EV` and `SINTP` |
| `vaporization_cstar_efficiency(...)` | both fractions and the c* efficiency | NACA TN 4219 Eq. (5) |
| `priem_effective_length(...)` | the earlier, propellant-independent effective length on the TN 3883 scale | NACA TN 3883 p. 5 |

In the source's units (L and r in inches, p_c in psia, H_v in Btu/lb, V in
ft/s):

```
L_gen = L (p_c/300)^0.66 / CR^0.44
        / [(1 - T/T_c)^0.4 (H_v/140)^0.8 (M/100)^0.35 (r/0.003)^1.45 (V/100)^0.75]
```

The pressure, drop-radius, velocity and `(1 - T_r)` exponents are the same
as in Priem's TN 3883 effective length; R-67 replaced the final-gas-velocity
term with the property and contraction-ratio factors.

## Verification

| Check | Result |
|---|---|
| Curve | the 35 nodes match an independent transcription of CR-168066 p. C-6; interpolation is linear as `SINTP` does; outside 0.01-100 in. it raises, as the program stops |
| Reference conditions (300 psia, 140 Btu/lb, M = 100, r = 0.003 in., V = 100 ft/s, CR = 1) | `L_gen` equals the chamber length to 1e-8 |
| Each exponent | doubling each input changes `L_gen` by exactly `2^exponent` |
| Exponents shared with TN 3883 | identical |
| c* efficiency | equals TN 4219 Eq. (5) to 1e-12; tends to 1 at full vaporization |

## Limits

- NASA TR R-67 (1960) is not on NTRS and was not read. The generalized length
  and the curve are quoted from CR-168066 (J. Fang, Aerojet, 1983), a
  secondary source. Its printed SI equation is dimensionally inconsistent with
  its own program (`V/328` with V in m/s, `H_v` "in J/kg" over 325.64), so the
  program's English-unit form is used.
- Only the cylindrical chamber length counts; CR-168066's eq. (1) adds a
  convergent-section term that its program leaves out.
- One average curve; the spray spread it belongs to is not stated (Priem's
  reports use a geometric standard deviation of 2.3).
- The mass-median drop radius is an input; no atomization correlation is
  bundled.
- Vaporized propellant is assumed to burn at once (TN 4219's model).

## References

- **NACA-TN-3985-PRIEM-1957**: <https://ntrs.nasa.gov/citations/19930084710>
- **NACA-TN-3883-PRIEM-1958**: <https://ntrs.nasa.gov/citations/19930084722>
- **NASA-CR-168066-AEROJET-1983**: <https://ntrs.nasa.gov/citations/19830023615>
- **NACA-TN4219**: <https://ntrs.nasa.gov/citations/19930085136>

Full records are in [REFERENCES.md](REFERENCES.md).
