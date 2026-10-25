"""Versioned, deterministic Rocket Propulsion Lab workspace documents."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field, replace
from math import isfinite
from typing import Any

from rocket_propulsion.core.errors import InputError
from rocket_propulsion.core.validation import DOCUMENT_ERRORS

CURRENT_SCHEMA_VERSION = 2


def _validate_json_tree(value: Any, path: str = "workspace") -> None:
    """Reject non-JSON values, non-string keys, NaN and infinities recursively.

    Parameters
    ----------
    value : object
        Candidate JSON tree.
    path : str, default 'workspace'
        Location prefix used in errors.

    Raises
    ------
    InputError
        At the first offending location.
    """

    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not isfinite(value):
            raise InputError(f"{path} must contain only finite numbers.", field=path)
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise InputError(f"{path} keys must be strings.", field=path)
            _validate_json_tree(item, f"{path}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _validate_json_tree(item, f"{path}[{index}]")
        return
    raise InputError(f"{path} contains a non-JSON value.", field=path)


@dataclass(frozen=True, slots=True)
class StudyCase:
    """One named input/result snapshot inside a workspace.

    Attributes
    ----------
    case_id : str
        Stable identifier of letters, digits, hyphens and underscores.
    name : str
        Display name.
    model : str
        Name of the model that produced the results.
    inputs, results : Mapping
        JSON-compatible model inputs and results.
    model_versions : Mapping of str to str
        Versions of the models used.
    warnings : tuple of str
        Warnings attached to the results.
    """

    case_id: str
    name: str
    model: str
    inputs: Mapping[str, Any]
    results: Mapping[str, Any]
    model_versions: Mapping[str, str] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Workspace:
    """Portable workspace document; large source tables are referenced, not embedded.

    Attributes
    ----------
    schema_version : int
        Document schema version.
    name : str
        Workspace name.
    cases : tuple of StudyCase
        Cases with unique identifiers.
    dataset_versions : Mapping of str to str
        Versions of referenced datasets.
    notes : str
        Free-text notes.
    """

    schema_version: int
    name: str
    cases: tuple[StudyCase, ...]
    dataset_versions: Mapping[str, str] = field(default_factory=dict)
    notes: str = ""


def _validate_case_id(case_id: str) -> None:
    """Refuse an empty case id or one with characters outside ``[A-Za-z0-9_-]``.

    Parameters
    ----------
    case_id : str

    Raises
    ------
    InputError
    """
    if not case_id or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"
        for character in case_id
    ):
        raise InputError(
            "Case id must contain only letters, numbers, hyphens, and underscores.",
            field="case_id",
        )


def create_workspace(name: str, cases: tuple[StudyCase, ...] = ()) -> Workspace:
    """Create a current-version workspace after enforcing stable case identifiers.

    Parameters
    ----------
    name : str
        Workspace name; surrounding whitespace is removed and it must not be empty.
    cases : tuple of StudyCase, default ()
        Cases with valid, unique identifiers.

    Returns
    -------
    Workspace

    Raises
    ------
    InputError
        If the name is empty or a case id is invalid or duplicated.
    """

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
    """Copy a case under a new identifier without sharing its mappings.

    Parameters
    ----------
    case : StudyCase
        Case to copy.
    case_id : str
        Identifier of the copy.
    name : str, optional
        Name of the copy; defaults to ``'<name> copy'``.

    Returns
    -------
    StudyCase

    Raises
    ------
    InputError
        If ``case_id`` is invalid.
    """

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
    """Serialize a workspace with sorted keys so equal workspaces give identical JSON.

    Parameters
    ----------
    workspace : Workspace
        Current-version workspace.
    indent : int, default 2
        JSON indentation.

    Returns
    -------
    str
        UTF-8-ready JSON text ending with a newline; NaN and infinities are refused
        [1]_.

    Raises
    ------
    InputError
        If the workspace is not current-version or contains non-JSON values.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format", 2017.
    """

    if workspace.schema_version != CURRENT_SCHEMA_VERSION:
        raise InputError("Only current-version workspaces can be serialized.")
    payload = asdict(workspace)
    _validate_json_tree(payload)
    return (
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )


def _migrate_v1(document: dict[str, Any]) -> dict[str, Any]:
    """Convert a schema-1 document (single ``case``) to the schema-2 layout.

    Parameters
    ----------
    document : dict

    Returns
    -------
    dict
    """
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
    """Migrate a supported earlier document shape to the current schema.

    Parameters
    ----------
    document : Mapping
        Parsed workspace document with ``schema_version`` or ``version``.

    Returns
    -------
    dict
        Document in the current schema layout.

    Raises
    ------
    InputError
        ``unsupported_workspace_version`` for any other version.
    """

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
    """Parse, migrate and validate a workspace JSON document.

    A workspace is a user-editable study file without an integrity hash, so a changed
    document yields a changed workspace; malformed content is refused.

    Parameters
    ----------
    text : str
        JSON document [1]_.

    Returns
    -------
    Workspace

    Raises
    ------
    InputError
        If the text is not JSON, the root is not an object, a value is not JSON, the
        schema version is unsupported, or the content is malformed.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format", 2017.
    """

    if not isinstance(text, str):
        raise InputError("Workspace must be JSON text.", code="invalid_workspace_json")
    try:
        raw = json.loads(text)
    except (ValueError, RecursionError) as error:
        raise InputError("Workspace is not valid JSON.", code="invalid_workspace_json") from error
    if not isinstance(raw, dict):
        raise InputError("Workspace root must be a JSON object.")
    _validate_json_tree(raw)
    try:
        return _workspace_from_document(raw)
    except InputError:
        raise
    except DOCUMENT_ERRORS as error:
        raise InputError("Workspace content is malformed.", field="workspace") from error


def _workspace_from_document(raw: dict) -> Workspace:
    """Build a :class:`Workspace` from a parsed, JSON-validated document.

    Parameters
    ----------
    raw : dict

    Returns
    -------
    Workspace

    Raises
    ------
    InputError
        If the document or a case is incomplete or malformed.
    """
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
