"""Sourced reaction mechanisms for :mod:`.kinetics`.

The only mechanism provided is the hydrogen-oxygen subset of GRI-Mech 3.0
[1]_: every reaction of ``grimech30.dat`` whose species all belong to
``{H2, H, O, O2, OH, H2O, HO2, H2O2}``, with the rate parameters copied from
the file in its CHEMKIN-II default units (cm, mol, s, cal/mol) [2]_ and
converted to SI by :meth:`.kinetics.ArrheniusRate.from_chemkin`. Each reaction
records the line of ``grimech30.dat`` it was read from (file SHA-256
``a7726208619d27b50e43b9fb6b27dc0a737003c049cc79e0a849bc75301cd763``,
downloaded 2026-10-25).

Adaptations, all recorded in :attr:`.kinetics.Mechanism.notes`:

* GRI-Mech's own thermochemistry (``thermo30.dat``) is not used; reverse rates
  come from the repository's NASA Glenn polynomials [3]_.
* Collision efficiencies of species outside the set (CH4, CO, CO2, C2H6, N2,
  AR) are dropped, because those species have zero concentration here; the
  collider-specific reactions ``H+O2+N2`` (line 62) and ``H+O2+AR`` (line 63)
  are omitted for the same reason.
* GRI-Mech 3.0 was optimized for natural-gas combustion [1]_; its use for
  rocket-nozzle recombination above about 2500 K and at tens of bar is outside
  the target set the authors list, so results carry that model-form
  uncertainty.

References
----------
.. [1] GRI-MECH-3.0. G. P. Smith, D. M. Golden, M. Frenklach, N. W. Moriarty,
       B. Eiteneer, M. Goldenberg, C. T. Bowman, R. K. Hanson, S. Song,
       W. C. Gardiner Jr., V. V. Lissianski and Z. Qin, "GRI-Mech 3.0",
       mechanism file grimech30.dat and documentation page.
       http://combustion.berkeley.edu/gri-mech/version30/text30.html
.. [2] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, section 4.14.4.
       https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
.. [3] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from functools import lru_cache

from .kinetics import ArrheniusRate, Mechanism, Reaction, ThirdBody, TroeFalloff

GRI30_URL = "http://combustion.berkeley.edu/gri-mech/version30/files30/grimech30.dat"
GRI30_SHA256 = "a7726208619d27b50e43b9fb6b27dc0a737003c049cc79e0a849bc75301cd763"
H2O2_SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2")

# (line, equation, reactants, products, A [cm-mol-s], b, E [cal/mol],
#  third-body efficiencies or None, duplicate). Values as printed in
#  grimech30.dat; efficiencies of species outside H2O2_SPECIES are dropped.
_GRI30_H2O2_ROWS: tuple[tuple, ...] = (
    (22, "2O+M<=>O2+M", (("O", 2),), (("O2", 1),), 1.200e17, -1.000, 0.00,
     (("H2", 2.40), ("H2O", 15.40)), False),
    (24, "O+H+M<=>OH+M", (("O", 1), ("H", 1)), (("OH", 1),), 5.000e17, -1.000, 0.00,
     (("H2", 2.00), ("H2O", 6.00)), False),
    (26, "O+H2<=>H+OH", (("O", 1), ("H2", 1)), (("H", 1), ("OH", 1)),
     3.870e04, 2.700, 6260.00, None, False),
    (27, "O+HO2<=>OH+O2", (("O", 1), ("HO2", 1)), (("OH", 1), ("O2", 1)),
     2.000e13, 0.000, 0.00, None, False),
    (28, "O+H2O2<=>OH+HO2", (("O", 1), ("H2O2", 1)), (("OH", 1), ("HO2", 1)),
     9.630e06, 2.000, 4000.00, None, False),
    (58, "H+O2+M<=>HO2+M", (("H", 1), ("O2", 1)), (("HO2", 1),), 2.800e18, -0.860, 0.00,
     (("O2", 0.00), ("H2O", 0.00)), False),
    (60, "H+2O2<=>HO2+O2", (("H", 1), ("O2", 2)), (("HO2", 1), ("O2", 1)),
     2.080e19, -1.240, 0.00, None, False),
    (61, "H+O2+H2O<=>HO2+H2O", (("H", 1), ("O2", 1), ("H2O", 1)), (("HO2", 1), ("H2O", 1)),
     11.26e18, -0.760, 0.00, None, False),
    (64, "H+O2<=>O+OH", (("H", 1), ("O2", 1)), (("O", 1), ("OH", 1)),
     2.650e16, -0.6707, 17041.00, None, False),
    (65, "2H+M<=>H2+M", (("H", 2),), (("H2", 1),), 1.000e18, -1.000, 0.00,
     (("H2", 0.00), ("H2O", 0.00)), False),
    (67, "2H+H2<=>2H2", (("H", 2), ("H2", 1)), (("H2", 2),), 9.000e16, -0.600, 0.00,
     None, False),
    (68, "2H+H2O<=>H2+H2O", (("H", 2), ("H2O", 1)), (("H2", 1), ("H2O", 1)),
     6.000e19, -1.250, 0.00, None, False),
    (70, "H+OH+M<=>H2O+M", (("H", 1), ("OH", 1)), (("H2O", 1),), 2.200e22, -2.000, 0.00,
     (("H2", 0.73), ("H2O", 3.65)), False),
    (72, "H+HO2<=>O+H2O", (("H", 1), ("HO2", 1)), (("O", 1), ("H2O", 1)),
     3.970e12, 0.000, 671.00, None, False),
    (73, "H+HO2<=>O2+H2", (("H", 1), ("HO2", 1)), (("O2", 1), ("H2", 1)),
     4.480e13, 0.000, 1068.00, None, False),
    (74, "H+HO2<=>2OH", (("H", 1), ("HO2", 1)), (("OH", 2),), 0.840e14, 0.000, 635.00,
     None, False),
    (75, "H+H2O2<=>HO2+H2", (("H", 1), ("H2O2", 1)), (("HO2", 1), ("H2", 1)),
     1.210e07, 2.000, 5200.00, None, False),
    (76, "H+H2O2<=>OH+H2O", (("H", 1), ("H2O2", 1)), (("OH", 1), ("H2O", 1)),
     1.000e13, 0.000, 3600.00, None, False),
    (151, "OH+H2<=>H+H2O", (("OH", 1), ("H2", 1)), (("H", 1), ("H2O", 1)),
     2.160e08, 1.510, 3430.00, None, False),
    (156, "2OH<=>O+H2O", (("OH", 2),), (("O", 1), ("H2O", 1)), 3.570e04, 2.400, -2110.00,
     None, False),
    (157, "OH+HO2<=>O2+H2O", (("OH", 1), ("HO2", 1)), (("O2", 1), ("H2O", 1)),
     1.450e13, 0.000, -500.00, None, True),
    (159, "OH+H2O2<=>HO2+H2O", (("OH", 1), ("H2O2", 1)), (("HO2", 1), ("H2O", 1)),
     2.000e12, 0.000, 427.00, None, True),
    (161, "OH+H2O2<=>HO2+H2O", (("OH", 1), ("H2O2", 1)), (("HO2", 1), ("H2O", 1)),
     1.700e18, 0.000, 29410.00, None, True),
    (191, "2HO2<=>O2+H2O2", (("HO2", 2),), (("O2", 1), ("H2O2", 1)),
     1.300e11, 0.000, -1630.00, None, True),
    (193, "2HO2<=>O2+H2O2", (("HO2", 2),), (("O2", 1), ("H2O2", 1)),
     4.200e14, 0.000, 12000.00, None, True),
    (394, "OH+HO2<=>O2+H2O", (("OH", 1), ("HO2", 1)), (("O2", 1), ("H2O", 1)),
     0.500e16, 0.000, 17330.00, None, True),
)

# Falloff reaction, lines 152-155: 2OH(+M)<=>H2O2(+M); high-pressure limit,
# LOW, TROE (alpha, T***, T*, T**) and efficiencies.
_GRI30_FALLOFF = {
    "lines": "152-155",
    "equation": "2OH(+M)<=>H2O2(+M)",
    "high": (7.400e13, -0.370, 0.00),
    "low": (2.300e18, -0.900, -1700.00),
    "troe": (0.7346, 94.00, 1756.00, 5182.00),
    "efficiencies": (("H2", 2.00), ("H2O", 6.00)),
}


def _order(reactants: tuple[tuple[str, int], ...], third_body: bool) -> int:
    return sum(count for _, count in reactants) + (1 if third_body else 0)


@lru_cache(maxsize=1)
def gri30_h2o2_mechanism() -> Mechanism:
    """Return the H2-O2 subset of GRI-Mech 3.0 (27 reactions, 8 species) [1]_.

    Returns
    -------
    Mechanism
        Reactions in file order with SI rate parameters; each ``source`` names
        the ``grimech30.dat`` line(s). The ``+M`` reactions use the file's
        efficiencies for species in the set (default 1.0); ``H+2O2`` and
        ``H+O2+H2O`` are the file's collider-specific companions of
        ``H+O2+M`` (whose O2 and H2O efficiencies are zero). Duplicate
        reactions add.

    Notes
    -----
    CHEMKIN default units are cm-mol-s for ``A`` and cal/mol for ``E`` [2]_;
    the conversion uses the reaction order counting ``M`` once. The falloff
    reaction ``2OH(+M)<=>H2O2(+M)`` uses the Troe form with the file's four
    parameters [2]_.

    References
    ----------
    .. [1] GRI-MECH-3.0. G. P. Smith et al., "GRI-Mech 3.0", grimech30.dat.
           http://combustion.berkeley.edu/gri-mech/version30/text30.html
    .. [2] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, section 4.14.4.
           https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
    """

    reactions: list[Reaction] = []
    for (line, equation, reactants, products, a, b, e, efficiencies,
         duplicate) in _GRI30_H2O2_ROWS:
        third = efficiencies is not None
        reactions.append(Reaction(
            equation=equation,
            reactants=reactants,
            products=products,
            rate=ArrheniusRate.from_chemkin(a, b, e, order=_order(reactants, third)),
            third_body=ThirdBody(efficiencies) if third else None,
            duplicate=duplicate,
            source=f"GRI-Mech 3.0 grimech30.dat line {line}",
        ))
        if line == 151:
            data = _GRI30_FALLOFF
            reactions.append(Reaction(
                equation=data["equation"],
                reactants=(("OH", 2),),
                products=(("H2O2", 1),),
                rate=ArrheniusRate.from_chemkin(*data["high"], order=2),
                third_body=ThirdBody(data["efficiencies"]),
                falloff=TroeFalloff(ArrheniusRate.from_chemkin(*data["low"], order=3),
                                    *data["troe"]),
                source=f"GRI-Mech 3.0 grimech30.dat lines {data['lines']}",
            ))
    return Mechanism(
        name="GRI-Mech 3.0 H2-O2 subset",
        species=H2O2_SPECIES,
        reactions=tuple(reactions),
        source="GRI-Mech 3.0, grimech30.dat (sha256 " + GRI30_SHA256 + ")",
        source_url=GRI30_URL,
        notes=(
            "Reverse rates from the NASA Glenn polynomials, not GRI-Mech thermo30.dat.",
            "Efficiencies of CH4, CO, CO2, C2H6, N2 and AR dropped (species absent).",
            "Lines 62 (H+O2+N2) and 63 (H+O2+AR) omitted (colliders absent).",
            "GRI-Mech 3.0 targets natural-gas combustion; rocket-nozzle use is an "
            "extrapolation beyond its optimization targets.",
        ),
    )
