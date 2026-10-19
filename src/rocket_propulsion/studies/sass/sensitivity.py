"""Deterministic design sensitivity, keeping catalog data unchanged.

References
----------
NASA-SE-6: 6.2.1.2.4 and 6.8, sensitivity to requirements and assumptions.
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
    """A sourced change of mission requirement or bus mass, not a hardware extrapolation.

    References
    ----------
    NASA-SE-6: 6.2.1.2.4, functional/sensitivity analysis of changed requirements.
    """

    name: str
    delta_v_scale: CitedScalar
    bus_dry_mass_scale: CitedScalar

    def __post_init__(self) -> None:
        """Require identified and evidence-bearing scenario multipliers.

        References
        ----------
        NASA-SE-6: 6.8, assumptions and sensitivity in decision analysis.
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
    NASA-SE-6: 6.2/6.8, traceable hypothetical input changes.
    NASA-7009B: declared numerical domain and finite representation limits.
    """
    scaled = Fraction(value.value) * Fraction(factor.value)
    binary = binary_ledger((("scaled_quantity", scaled, value.unit),), "Design sensitivity")[0][1]
    return derived(binary, value.unit, (value, factor), methods, note)


def design_sensitivity(request: dict[str, object], cases: tuple[DesignCase, ...]
                        ) -> tuple[ProjectStudy, ...]:
    """Rerun explicit scenario changes; no probability or convergence claim is invented.

    All cases start from the same validated request. Catalogue thrust/Isp,
    selection anchors and margins are preserved. Scaling the mission request
    is a hypothetical requirement change, never a performance-map extrapolation.

    References
    ----------
    NASA-SE-6: 6.2.1.2.4 and 6.8.1.2.4, requirement and decision sensitivities.
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
        for item in variant["items"]:
            value = scalar_from_dict(item["nominal_delta_v"])
            item["nominal_delta_v"] = asdict(_scaled_quantity(
                value, case.delta_v_scale,
                methods, f"Design sensitivity case {case.name}; mission requirement change.",
            ))
        variant["title"] = f"{request['title']} / {case.name}"
        output.append(project_from_dict(variant))
    return tuple(output)
