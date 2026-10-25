"""Reject ambiguous input before any evidence envelope is interpreted.

References
----------
.. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; 6.6, controlled and interpretable technical inputs.
"""

import json

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import DOCUMENT_ERRORS, json_text


def unique_object(pairs: list[tuple[str, object]]) -> dict:
    """Build a JSON object, refusing duplicate fields at every nesting level.

    Parameters
    ----------
    pairs : list of (str, object)
        Key/value pairs from the JSON decoder.

    Returns
    -------
    dict

    Raises
    ------
    InputError
        On the first duplicate field.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.6, unambiguous input data; this is an implementation guard.
    """
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f"Duplicate JSON field: {key!r}.")
        result[key] = value
    return result


def invalid_constant(value: str):
    """Refuse the non-standard JSON constants ``NaN``, ``Infinity`` and ``-Infinity``.

    Parameters
    ----------
    value : str
        Constant name from the decoder.

    Raises
    ------
    InputError
        Always.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.6, numerical data integrity; implementation input guard.
    """
    raise InputError(f"Non-finite JSON constant: {value}.")


def read_project_document(value: str, build):
    """Parse a strict JSON request of at most 2 MiB and pass it to ``build``.

    Parameters
    ----------
    value : str
        JSON text.
    build : callable
        ``build(payload) -> object`` applied to the parsed root object.

    Returns
    -------
    object
        Result of ``build``.

    Raises
    ------
    InputError
        If the text is too large, malformed, has duplicate fields or non-finite
        constants, or the root is not an object.
    DomainError
        If ``build`` rejects the content.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.6, controlled technical inputs; no silent last-field override.
    """
    source = json_text(value, "SASS project")
    try:
        if len(source) > 2_097_152 or len(source.encode("utf-8")) > 2_097_152:
            raise InputError("SASS JSON exceeds the implementation's 2 MiB input limit.")
        payload = json.loads(source, object_pairs_hook=unique_object, parse_constant=invalid_constant)
        if not isinstance(payload, dict):
            raise InputError("SASS project root must be a JSON object.")
        return build(payload)
    except (InputError, DomainError):
        raise
    except (*DOCUMENT_ERRORS, UnicodeError) as error:
        raise InputError(f"SASS project is malformed ({type(error).__name__}).") from error
