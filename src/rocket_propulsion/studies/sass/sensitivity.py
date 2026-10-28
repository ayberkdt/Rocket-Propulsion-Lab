"""Deterministic design sensitivity, keeping catalog data unchanged.

References
----------
.. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; 6.2.1.2.4 and 6.8, sensitivity to requirements and assumptions.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import records, text

from .evidence import Citation, CitedScalar, derived, quantity
from .project import ProjectStudy, project_from_dict, scalar_from_dict
from .rational_numerics import binary_ledger


@dataclass(frozen=True, slots=True)
class DesignCase:
    """A sourced change of mission requirement or bus mass; not a hardware extrapolation.

    Attributes
    ----------
    name : str
        Unique case name.
    delta_v_scale : CitedScalar
        Multiplier of every nominal delta-v allocation; non-negative.
    bus_dry_mass_scale : CitedScalar
        Multiplier of the bus dry mass; positive.

    Raises
    ------
    DomainError
        If the name is empty or a multiplier is out of domain.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.2.1.2.4, functional/sensitivity analysis of changed
           requirements.
    """

    name: str
    delta_v_scale: CitedScalar
    bus_dry_mass_scale: CitedScalar

    def __post_init__(self) -> None:
        """Require identified and evidence-bearing scenario multipliers.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; 6.8, assumptions and sensitivity in decision analysis.
        """
        text(self.name, "sensitivity case name")
        quantity(self.delta_v_scale, "1", minimum=0)
        quantity(self.bus_dry_mass_scale, "1", strictly_positive=True)


def _scaled_quantity(value: CitedScalar, factor: CitedScalar, methods: tuple[Citation, ...],
                     note: str) -> CitedScalar:
    """Preserve positive requirement changes or reject unrepresentable outputs.

    An explicit zero factor is a valid hypothetical allocation change. A
    positive exact product must not silently become the same zero requirement.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.2/6.8, traceable hypothetical input changes.
    .. [2] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           declared numerical domain and finite representation limits.
    """
    scaled = Fraction(value.value) * Fraction(factor.value)
    binary = binary_ledger((("scaled_quantity", scaled, value.unit),), "Design sensitivity")[0][1]
    return derived(binary, value.unit, (value, factor), methods, note)


def design_sensitivity(request: dict[str, object], cases: tuple[DesignCase, ...]
                        ) -> tuple[ProjectStudy, ...]:
    """Rerun a validated request under explicit scenario changes.

    Parameters
    ----------
    request : dict
        Valid project request (see :func:`project_from_dict`).
    cases : tuple of DesignCase
        Non-empty, uniquely named scenarios.

    Returns
    -------
    tuple of ProjectStudy
        One study per case, in input order.

    Raises
    ------
    DomainError
        If the request or a case is invalid or a scaled quantity is not
        representable.

    Notes
    -----
    All cases start from the same validated request. Catalogue thrust and Isp,
    selection anchors and margins are preserved; scaling is a hypothetical
    requirement change [1]_, never a performance-map extrapolation. No probability or
    convergence claim is made.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.2.1.2.4 and 6.8.1.2.4, requirement and decision sensitivities.
    """
    project_from_dict(request)
    scenarios = records(cases, DesignCase, "sensitivity cases")
    if not scenarios or len({x.name for x in scenarios}) != len(scenarios):
        raise DomainError("Use non-empty, uniquely identified design sensitivity cases.")
    output = []
    methods = (Citation("NASA-SE-6", "6.2/6.8: hypothetical requirement/input sensitivity"),)
    for case in scenarios:
        variant = deepcopy(request)
        mass = scalar_from_dict(variant["bus_dry_mass_excluding_propulsion"])
        variant["bus_dry_mass_excluding_propulsion"] = asdict(_scaled_quantity(
            mass, case.bus_dry_mass_scale,
            methods, f"Design sensitivity case {case.name}; not measured hardware scaling.",
        ))
        items = variant["items"]
        if not isinstance(items, list):  # already enforced by project_from_dict above
            raise DomainError("items must be a JSON array.")
        for item in items:
            value =scalar_from_dict(item["nominal_delta_v"])
            item["nominal_delta_v"] = asdict(_scaled_quantity(
                value, case.delta_v_scale,
                methods, f"Design sensitivity case {case.name}; mission requirement change.",
            ))
        variant["title"] = f"{request['title']} / {case.name}"
        output.append(project_from_dict(variant))
    return tuple(output)
