"""Versioned, deterministic Rocket Propulsion Lab workspace documents."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field, replace
from typing import Any

from rocket_propulsion.core.errors import InputError

CURRENT_SCHEMA_VERSION = 2


@dataclass(frozen=True, slots=True)
class StudyCase:
    """One named input/result snapshot inside a workspace."""

    case_id: str
    name: str
    model: str
    inputs: Mapping[str, Any]
    results: Mapping[str, Any]
    model_versions: Mapping[str, str] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Workspace:
    """Portable workspace schema; large source tables are referenced, not embedded."""

    schema_version: int
    name: str
    cases: tuple[StudyCase, ...]
    dataset_versions: Mapping[str, str] = field(default_factory=dict)
    notes: str = ""


def _validate_case_id(case_id: str) -> None:
    if not case_id or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for character in case_id):
        raise InputError(
            "Case id must contain only letters, numbers, hyphens, and underscores.",
            field="case_id",
        )


def create_workspace(name: str, cases: tuple[StudyCase, ...] = ()) -> Workspace:
    """Create a current-version workspace after enforcing stable case identifiers."""

    if not name.strip():
        raise InputError("Workspace name cannot be empty.", field="name")
    identifiers: set[str] = set()
    for case in cases:
        _validate_case_id(case.case_id)
        if case.case_id in identifiers:
            raise InputError("Workspace case ids must be unique.", field="case_id")
        identifiers.add(case.case_id)
    return Workspace(CURRENT_SCHEMA_VERSION, name.strip(), cases)


def duplicate_case(case: StudyCase, *, case_id: str, name: str | None = None) -> StudyCase:
    """Duplicate a case without sharing mutable input/result dictionaries."""

    _validate_case_id(case_id)
    return replace(
        case,
        case_id=case_id,
        name=name.strip() if name is not None else f"{case.name} copy",
        inputs=dict(case.inputs),
        results=dict(case.results),
        model_versions=dict(case.model_versions),
        warnings=tuple(case.warnings),
    )


def workspace_to_json(workspace: Workspace, *, indent: int = 2) -> str:
    """Serialize with sorted keys so equal workspaces produce byte-identical JSON."""

    if workspace.schema_version != CURRENT_SCHEMA_VERSION:
        raise InputError("Only current-version workspaces can be serialized.")
    return json.dumps(
        asdict(workspace),
        ensure_ascii=False,
        indent=indent,
        sort_keys=True,
        allow_nan=False,
    ) + "\n"


def _migrate_v1(document: dict[str, Any]) -> dict[str, Any]:
    raw_case = document.get("case")
    cases = [] if raw_case is None else [raw_case]
    return {
        "schema_version": 2,
        "name": document.get("name", "Migrated workspace"),
        "cases": cases,
        "dataset_versions": document.get("datasets", {}),
        "notes": document.get("notes", ""),
    }


def migrate_workspace_document(document: Mapping[str, Any]) -> dict[str, Any]:
    """Migrate supported historical shapes to the current schema."""

    migrated = dict(document)
    version = migrated.get("schema_version", migrated.get("version"))
    if version == 1:
        migrated = _migrate_v1(migrated)
        version = 2
    if version != CURRENT_SCHEMA_VERSION:
        raise InputError(
            f"Unsupported workspace schema version: {version}.",
            code="unsupported_workspace_version",
            field="schema_version",
        )
    return migrated


def workspace_from_json(text: str) -> Workspace:
    """Parse, migrate, and validate a workspace JSON document."""

    try:
        raw = json.loads(text)
    except (TypeError, json.JSONDecodeError) as error:
        raise InputError("Workspace is not valid JSON.", code="invalid_workspace_json") from error
    if not isinstance(raw, dict):
        raise InputError("Workspace root must be a JSON object.")
    document = migrate_workspace_document(raw)
    raw_cases = document.get("cases")
    if not isinstance(raw_cases, list):
        raise InputError("Workspace cases must be a list.", field="cases")
    cases: list[StudyCase] = []
    for raw_case in raw_cases:
        if not isinstance(raw_case, dict):
            raise InputError("Each workspace case must be an object.", field="cases")
        try:
            cases.append(
                StudyCase(
                    case_id=str(raw_case["case_id"]),
                    name=str(raw_case["name"]),
                    model=str(raw_case["model"]),
                    inputs=dict(raw_case.get("inputs", {})),
                    results=dict(raw_case.get("results", {})),
                    model_versions=dict(raw_case.get("model_versions", {})),
                    warnings=tuple(str(item) for item in raw_case.get("warnings", [])),
                )
            )
        except (KeyError, TypeError, ValueError) as error:
            raise InputError("Workspace case is incomplete.", field="cases") from error
    workspace = create_workspace(str(document.get("name", "")), tuple(cases))
    return replace(
        workspace,
        dataset_versions=dict(document.get("dataset_versions", {})),
        notes=str(document.get("notes", "")),
    )
