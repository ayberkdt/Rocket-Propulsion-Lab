"""NASA CEA input, subprocess, parser, and deterministic cache adapter.

NASA CEA itself is not distributed with this project.  Users point this
adapter at their own installation; the ideal-gas and local experimental
providers remain fully usable when CEA is absent.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from math import isfinite
from pathlib import Path
from typing import Literal

from rocket_propulsion.core.errors import (
    ConvergenceError,
    DomainError,
    FeatureUnavailableError,
    InputError,
)
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import (
    EquilibriumConvergence,
    EquilibriumProblem,
    EquilibriumResiduals,
    EquilibriumResult,
    EquilibriumSpeciesAmount,
)
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import get_species

CEA_MANUAL_URL = "https://ntrs.nasa.gov/citations/19960044559"
_FORTRAN_EXPONENT = re.compile(r"(?<![EeDd])([+-]\d+)$")


@dataclass(frozen=True, slots=True)
class CeaReactant:
    """One fuel, oxidizer, or unclassified CEA reactant record."""

    name: str
    role: Literal["fuel", "oxid", "name"] = "name"
    amount: float = 1.0
    basis: Literal["moles", "wt%"] = "moles"
    temperature_k: float | None = None


@dataclass(frozen=True, slots=True)
class CeaRocketProblem:
    """Rocket-mode CEA request with explicit chemistry and exit stations."""

    fuel: CeaReactant
    oxidizer: CeaReactant
    oxidizer_fuel_ratio: float
    chamber_pressure_pa: float
    area_ratios: tuple[float, ...] = (10.0, 25.0, 50.0, 75.0)
    chemistry: Literal["equilibrium", "frozen", "both"] = "equilibrium"
    freezing_station: Literal[1, 2] = 2
    case_name: str = "rplab"


@dataclass(frozen=True, slots=True)
class CeaPlotTable:
    """Numeric rows from a CEA ``.plt`` file with caller-known columns."""

    columns: tuple[str, ...]
    rows: tuple[dict[str, float], ...]


@dataclass(frozen=True, slots=True)
class CeaRunArtifact:
    """Traceable raw and parsed result from one external CEA execution."""

    schema_version: str
    provider: str
    provider_version: str
    cache_key: str
    cached: bool
    input_text: str
    output_text: str
    plot_text: str
    plot: CeaPlotTable
    warnings: tuple[str, ...]
    source: str
    source_url: str


def parse_fortran_number(token: str) -> float:
    """Parse CEA's E/D notation and compact exponents such as ``9.18-05``."""

    normalized = token.strip().replace("D", "E").replace("d", "e")
    if not normalized:
        raise ValueError("Empty numeric token.")
    if "e" not in normalized.lower():
        normalized = _FORTRAN_EXPONENT.sub(r"E\1", normalized)
    value = float(normalized)
    if not isfinite(value):
        raise ValueError(f"Non-finite CEA value: {token!r}")
    return value


def parse_cea_plot(text: str, columns: Iterable[str]) -> CeaPlotTable:
    """Parse the numeric plot file emitted for an explicitly ordered schema."""

    schema = tuple(str(column).strip() for column in columns)
    if not schema or any(not column for column in schema):
        raise InputError("CEA plot columns must be non-empty.", field="columns")
    rows: list[dict[str, float]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "!")):
            continue
        tokens = stripped.replace(",", " ").split()
        if len(tokens) != len(schema):
            raise InputError(
                "CEA plot row does not match the requested schema.",
                code="cea_plot_shape",
                details={
                    "line": line_number,
                    "expected_columns": len(schema),
                    "actual_columns": len(tokens),
                },
            )
        try:
            row = {
                column: parse_fortran_number(token)
                for column, token in zip(schema, tokens, strict=True)
            }
        except ValueError as error:
            raise InputError(
                "CEA plot contains an invalid numeric value.",
                code="cea_plot_number",
                details={"line": line_number},
            ) from error
        rows.append(row)
    if not rows:
        raise InputError("CEA plot file contains no data rows.", code="cea_plot_empty")
    return CeaPlotTable(columns=schema, rows=tuple(rows))


def cea_output_warnings(text: str) -> tuple[str, ...]:
    """Extract CEA warnings while treating non-convergence as a hard error."""

    warnings = tuple(
        line.strip()
        for line in text.splitlines()
        if "WARNING" in line.upper() or "ERROR" in line.upper()
    )
    failure = next(
        (
            line
            for line in warnings
            if "NO CONVERGENCE" in line.upper()
            or "ANSWERS NOT RELIABLE" in line.upper()
            or "FATAL" in line.upper()
        ),
        None,
    )
    if failure is not None:
        raise ConvergenceError(
            "NASA CEA reported a non-converged solution.",
            code="cea_non_convergence",
            details={"cea_message": failure},
        )
    return warnings


def _reactant_line(reactant: CeaReactant) -> str:
    role = reactant.role.strip().lower()
    if role not in {"fuel", "oxid", "name"}:
        raise InputError("CEA reactant role must be fuel, oxid, or name.", field="role")
    basis = reactant.basis.strip().lower()
    if basis not in {"moles", "wt%"}:
        raise InputError("CEA amount basis must be moles or wt%.", field="basis")
    if not reactant.name.strip() or any(character.isspace() for character in reactant.name):
        raise InputError(
            "CEA reactant names must be non-empty and contain no spaces.", field="name"
        )
    if not isfinite(reactant.amount) or reactant.amount <= 0.0:
        raise DomainError("CEA reactant amount must be positive.", field="amount")
    line = f"  {role} = {reactant.name} {basis}={reactant.amount:.12g}"
    if reactant.temperature_k is not None:
        if not isfinite(reactant.temperature_k) or reactant.temperature_k <= 0.0:
            raise DomainError("CEA reactant temperature must be positive.")
        line += f" t,k={reactant.temperature_k:.12g}"
    return line


def render_cea_equilibrium_input(
    problem: EquilibriumProblem,
    *,
    plot_columns: tuple[str, ...] = ("t", "p", "mw", "gam", "cp", "h", "s"),
) -> str:
    """Render a generic TP, HP, or UV input deck from the common model."""

    kind = problem.problem_type.strip().lower()
    if kind not in {"tp", "hp", "uv"}:
        raise InputError("CEA supports TP, HP, and UV through this adapter.")
    reactant_temperature = problem.reactant_temperature_k
    reactants = [
        CeaReactant(
            name=reactant.species,
            role="name",
            amount=reactant.moles,
            basis="moles",
            temperature_k=reactant_temperature if kind in {"hp", "uv"} else None,
        )
        for reactant in problem.reactants
    ]
    problem_parts = [f"problem {kind}"]
    if kind in {"tp", "hp"}:
        if problem.pressure_pa is None:
            raise InputError("CEA TP/HP requires pressure.", field="pressure_pa")
        problem_parts.append(f"p,bar={problem.pressure_pa / 100_000.0:.12g}")
    if kind == "tp":
        if problem.temperature_k is None:
            raise InputError("CEA TP requires temperature.", field="temperature_k")
        problem_parts.append(f"t,k={problem.temperature_k:.12g}")
    if kind == "hp" and problem.target_enthalpy_j is not None:
        total_mass = sum(
            reactant.moles * get_species(reactant.species).molecular_mass_kg_mol
            for reactant in problem.reactants
        )
        h_over_r = problem.target_enthalpy_j / (
            total_mass * UNIVERSAL_GAS_CONSTANT_J_MOL_K * 1000.0
        )
        problem_parts.append(f"h/r={h_over_r:.12g}")
    if kind == "uv":
        if problem.volume_m3 is None or problem.volume_m3 <= 0.0:
            raise InputError("CEA UV requires volume.", field="volume_m3")
        total_mass = sum(
            reactant.moles * get_species(reactant.species).molecular_mass_kg_mol
            for reactant in problem.reactants
        )
        problem_parts.append(f"rho,kg/m**3={total_mass / problem.volume_m3:.12g}")
        if problem.target_internal_energy_j is not None:
            u_over_r = problem.target_internal_energy_j / (
                total_mass * UNIVERSAL_GAS_CONSTANT_J_MOL_K * 1000.0
            )
            problem_parts.append(f"u/r={u_over_r:.12g}")
    lines = [" ".join(problem_parts), "reactants"]
    lines.extend(_reactant_line(reactant) for reactant in reactants)
    if problem.candidate_species is not None:
        lines.append(f"only {' '.join(problem.candidate_species)}")
    lines.append(f"output siunits trace=1.e-12 plot {' '.join(plot_columns)}")
    lines.append("end")
    return "\n".join(lines) + "\n"


def render_cea_rocket_input(
    problem: CeaRocketProblem,
    *,
    plot_columns: tuple[str, ...] = (
        "aeat",
        "t",
        "p",
        "mw",
        "gam",
        "cstar",
        "cf",
        "ivac",
        "isp",
        "mach",
    ),
) -> str:
    """Render a NASA RP-1311 rocket input deck."""

    if problem.fuel.role != "fuel" or problem.oxidizer.role != "oxid":
        raise InputError("Rocket reactants must use fuel and oxid roles.")
    if not isfinite(problem.oxidizer_fuel_ratio) or problem.oxidizer_fuel_ratio <= 0.0:
        raise DomainError("Rocket O/F ratio must be positive.", field="oxidizer_fuel_ratio")
    if not isfinite(problem.chamber_pressure_pa) or problem.chamber_pressure_pa <= 0.0:
        raise DomainError("Chamber pressure must be positive.", field="chamber_pressure_pa")
    if not problem.area_ratios or any(
        not isfinite(value) or value <= 1.0 for value in problem.area_ratios
    ):
        raise DomainError("Rocket area ratios must all exceed one.", field="area_ratios")
    chemistry = problem.chemistry.strip().lower()
    if chemistry not in {"equilibrium", "frozen", "both"}:
        raise InputError("Rocket chemistry must be equilibrium, frozen, or both.")
    chemistry_text = chemistry
    if chemistry == "both":
        chemistry_text = f"equilibrium frozen nfz={problem.freezing_station}"
    areas = ",".join(f"{value:.12g}" for value in problem.area_ratios)
    lines = [
        "reactants",
        _reactant_line(problem.fuel),
        _reactant_line(problem.oxidizer),
        (
            f"problem rocket case={problem.case_name} {chemistry_text} "
            f"o/f={problem.oxidizer_fuel_ratio:.12g} "
            f"p,bar={problem.chamber_pressure_pa / 100_000.0:.12g} "
            f"supersonic,ae/at={areas}"
        ),
        f"output siunits trace=1.e-12 plot {' '.join(plot_columns)}",
        "end",
    ]
    return "\n".join(lines) + "\n"


class CeaSubprocessAdapter:
    """Run a user-supplied NASA CEA executable in an isolated directory."""

    name = "NASA CEA"

    def __init__(
        self,
        executable: str | Path | None,
        *,
        version: str = "user-installation",
        cache_directory: str | Path | None = None,
        timeout_seconds: float = 60.0,
        invocation: Literal["stdin-prefix", "input-argument"] = "stdin-prefix",
    ) -> None:
        self.executable = executable
        self.version = version
        self.cache_directory = (
            Path(cache_directory).expanduser().resolve()
            if cache_directory is not None
            else None
        )
        self.timeout_seconds = timeout_seconds
        self.invocation = invocation

    def _resolved_executable(self) -> Path:
        if self.executable is None:
            raise FeatureUnavailableError(
                "NASA CEA executable is not configured. Set an explicit CEA executable path.",
                details={"documentation": CEA_MANUAL_URL},
            )
        raw = str(self.executable)
        found = shutil.which(raw)
        path = Path(found if found is not None else raw).expanduser().resolve()
        if not path.is_file():
            raise FeatureUnavailableError(
                "Configured NASA CEA executable was not found.",
                details={"path": str(path), "documentation": CEA_MANUAL_URL},
            )
        return path

    def _cache_key(self, input_text: str, columns: tuple[str, ...]) -> str:
        document = json.dumps(
            {
                "adapter_schema": "1.0",
                "provider_version": self.version,
                "input": input_text,
                "columns": columns,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        return hashlib.sha256(document.encode("utf-8")).hexdigest()

    def _cache_path(self, key: str) -> Path | None:
        if self.cache_directory is None:
            return None
        return self.cache_directory / f"{key}.json"

    def _from_cache(self, path: Path) -> CeaRunArtifact:
        payload = json.loads(path.read_text(encoding="utf-8"))
        plot = CeaPlotTable(
            columns=tuple(payload["plot"]["columns"]),
            rows=tuple(payload["plot"]["rows"]),
        )
        return CeaRunArtifact(
            schema_version=payload["schema_version"],
            provider=payload["provider"],
            provider_version=payload["provider_version"],
            cache_key=payload["cache_key"],
            cached=True,
            input_text=payload["input_text"],
            output_text=payload["output_text"],
            plot_text=payload["plot_text"],
            plot=plot,
            warnings=tuple(payload["warnings"]),
            source=payload["source"],
            source_url=payload["source_url"],
        )

    def run(
        self,
        input_text: str,
        *,
        plot_columns: tuple[str, ...],
    ) -> CeaRunArtifact:
        """Execute CEA, parse its plot file, and cache the complete trace."""

        executable = self._resolved_executable()
        if self.invocation not in {"stdin-prefix", "input-argument"}:
            raise InputError("Unknown CEA invocation mode.", field="invocation")
        key = self._cache_key(input_text, plot_columns)
        cache_path = self._cache_path(key)
        if cache_path is not None and cache_path.is_file():
            return self._from_cache(cache_path)

        with tempfile.TemporaryDirectory(prefix="rplab-cea-") as directory:
            work = Path(directory)
            prefix = "case"
            input_path = work / f"{prefix}.inp"
            input_path.write_text(input_text, encoding="ascii")
            for library_name in ("thermo.lib", "trans.lib"):
                source = executable.parent / library_name
                if source.is_file():
                    shutil.copy2(source, work / library_name)
            command = [str(executable)]
            standard_input = f"{prefix}\n"
            if self.invocation == "input-argument":
                command.append(str(input_path))
                standard_input = None
            try:
                completed = subprocess.run(
                    command,
                    cwd=work,
                    input=standard_input,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout_seconds,
                    check=False,
                )
            except subprocess.TimeoutExpired as error:
                raise ConvergenceError(
                    "NASA CEA execution exceeded the configured timeout.",
                    code="cea_timeout",
                    details={"timeout_seconds": self.timeout_seconds},
                ) from error
            output_path = work / f"{prefix}.out"
            plot_path = work / f"{prefix}.plt"
            output_text = (
                output_path.read_text(encoding="utf-8", errors="replace")
                if output_path.is_file()
                else completed.stdout
            )
            if completed.returncode != 0:
                raise ConvergenceError(
                    "NASA CEA executable returned a failure status.",
                    code="cea_process_failure",
                    details={
                        "returncode": completed.returncode,
                        "stderr": completed.stderr[-2000:],
                    },
                )
            if not plot_path.is_file():
                raise ConvergenceError(
                    "NASA CEA did not produce the requested plot file.",
                    code="cea_missing_plot",
                    details={"stdout": completed.stdout[-2000:]},
                )
            plot_text = plot_path.read_text(encoding="utf-8", errors="replace")

        warnings = cea_output_warnings(output_text)
        plot = parse_cea_plot(plot_text, plot_columns)
        artifact = CeaRunArtifact(
            schema_version="1.0",
            provider=self.name,
            provider_version=self.version,
            cache_key=key,
            cached=False,
            input_text=input_text,
            output_text=output_text,
            plot_text=plot_text,
            plot=plot,
            warnings=warnings,
            source="NASA RP-1311 Part II, CEA User's Manual",
            source_url=CEA_MANUAL_URL,
        )
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(
                json.dumps(asdict(artifact), indent=2, allow_nan=False),
                encoding="utf-8",
            )
        return artifact


class CeaEquilibriumProvider:
    """Adapt a CEA installation to the common equilibrium result contract."""

    name = "NASA CEA"

    def __init__(self, adapter: CeaSubprocessAdapter) -> None:
        self.adapter = adapter
        self.version = adapter.version

    def _candidate_keys(self, problem: EquilibriumProblem) -> tuple[str, ...]:
        elements = {
            element
            for reactant in problem.reactants
            for element, _ in get_species(reactant.species).formula
        }
        if problem.candidate_species is not None:
            return tuple(key.upper() for key in problem.candidate_species)
        from .species import species_catalog

        return tuple(
            species.key
            for species in species_catalog().values()
            if set(dict(species.formula)).issubset(elements)
        )

    def solve(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Execute CEA and convert its numeric plot to versioned Python data."""

        candidates = self._candidate_keys(problem)
        resolved_problem = replace(problem, candidate_species=candidates)
        base_columns = ("t", "p", "mw", "gam", "cp", "h", "s")
        plot_columns = base_columns + candidates
        deck = render_cea_equilibrium_input(
            resolved_problem, plot_columns=plot_columns
        )
        artifact = self.adapter.run(deck, plot_columns=plot_columns)
        row = artifact.plot.rows[-1]
        raw_fractions = {key: max(0.0, row[key]) for key in candidates}
        fraction_sum = sum(raw_fractions.values())
        if fraction_sum <= 0.0:
            raise ConvergenceError(
                "NASA CEA returned no retained product fractions.",
                code="cea_empty_composition",
            )
        mole_fractions = {
            key: value / fraction_sum for key, value in raw_fractions.items() if value > 0.0
        }
        total_mass = sum(
            reactant.moles * get_species(reactant.species).molecular_mass_kg_mol
            for reactant in problem.reactants
        )
        molecular_mass = row["mw"] / 1000.0
        total_moles = total_mass / molecular_mass
        species_moles = {
            key: total_moles * fraction for key, fraction in mole_fractions.items()
        }
        species_masses = {
            key: moles * get_species(key).molecular_mass_kg_mol
            for key, moles in species_moles.items()
        }
        mass_sum = sum(species_masses.values())
        inventory: dict[str, float] = {}
        for reactant in problem.reactants:
            for element, coefficient in get_species(reactant.species).formula:
                inventory[element] = inventory.get(element, 0.0) + (
                    reactant.moles * coefficient
                )
        calculated = {
            element: sum(
                moles * dict(get_species(key).formula).get(element, 0.0)
                for key, moles in species_moles.items()
            )
            for element in inventory
        }
        element_residual = max(
            abs(calculated[element] - target) / target
            for element, target in inventory.items()
        )
        temperature_k = row["t"]
        pressure_pa = row["p"] * 100_000.0
        cp = row["cp"] * 1000.0
        gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / molecular_mass
        cv = cp - gas_constant
        if cv <= 0.0:
            raise DomainError(
                "NASA CEA returned non-positive mixture cv.",
                code="cea_invalid_property",
                details={"cp_j_kg_k": cp, "gas_constant_j_kg_k": gas_constant},
            )
        enthalpy = row["h"] * 1000.0 * total_mass
        internal_energy = enthalpy - (
            total_moles * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature_k
        )
        entropy = row["s"] * 1000.0 * total_mass
        warnings = tuple(
            WarningMessage(code="cea_warning", message=message)
            for message in artifact.warnings
        )
        residuals = EquilibriumResiduals(
            maximum_element_relative=element_residual,
            mole_fraction_sum=abs(fraction_sum - 1.0),
        )
        return EquilibriumResult(
            schema_version="1.0",
            provider=self.name,
            provider_version=self.version,
            problem_type=problem.problem_type,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            total_moles=total_moles,
            total_mass_kg=total_mass,
            molecular_mass_kg_mol=molecular_mass,
            gas_constant_j_kg_k=gas_constant,
            cp_j_kg_k=cp,
            cv_j_kg_k=cv,
            gamma=row["gam"],
            enthalpy_j=enthalpy,
            internal_energy_j=internal_energy,
            entropy_j_k=entropy,
            species=tuple(
                EquilibriumSpeciesAmount(
                    species=key,
                    moles=species_moles[key],
                    mole_fraction=fraction,
                    mass_fraction=species_masses[key] / mass_sum,
                )
                for key, fraction in sorted(
                    mole_fractions.items(), key=lambda item: -item[1]
                )
            ),
            element_moles=tuple(sorted(inventory.items())),
            residuals=residuals,
            convergence=EquilibriumConvergence(
                converged=True,
                composition_iterations=0,
                temperature_iterations=0,
                pressure_iterations=0,
                maximum_residual=max(
                    residuals.maximum_element_relative,
                    residuals.mole_fraction_sum,
                ),
            ),
            assumptions=(
                "NASA CEA equilibrium model",
                "product species restricted to the configured candidate set",
            ),
            warnings=warnings,
            source="NASA RP-1311 Part I and Part II",
            source_url=CEA_MANUAL_URL,
        )

