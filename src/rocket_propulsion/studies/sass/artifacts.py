"""Staged report writes with alias protection and best-effort rollback.

References
----------
.. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; 6.5/6.6 configuration and technical-data integrity.
"""

from pathlib import Path
from tempfile import NamedTemporaryFile

from rocket_propulsion.core.errors import DomainError


def distinct_artifact_paths(paths: tuple[Path, ...]) -> None:
    """Refuse path aliases, existing hard links and non-file destinations.

    Parameters
    ----------
    paths : tuple of pathlib.Path
        Request, JSON output and report paths.

    Raises
    ------
    DomainError
        If an existing path is not a regular file or two paths resolve to the same
        file.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.5/6.6 preserve source and result configurations.
    """
    for i, path in enumerate(paths):
        if path.exists() and not path.is_file():
            raise DomainError("Input and report destinations must be regular files.")
        for other in paths[:i]:
            if path.resolve() == other.resolve() or (
                path.exists() and other.exists() and path.samefile(other)
            ):
                raise DomainError("Input, JSON output and HTML report must be distinct paths.")


def write_artifact_pair(artifacts: tuple[tuple[Path, str], ...]) -> None:
    """Stage all files before replacing any, rolling back on ordinary write errors.

    Parameters
    ----------
    artifacts : tuple of (pathlib.Path, str)
        Destination paths and UTF-8 text contents.

    Raises
    ------
    OSError
        If staging or replacement fails and the rollback succeeds.
    DomainError
        If the rollback itself fails.

    Notes
    -----
    Each replacement is atomic on its filesystem. A power loss between replacements
    can still leave a mixed pair; the embedded study hash detects that case.
    Concurrent writers are outside the permitted use of the command-line tool.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.5/6.6 retain internally identifiable result snapshots.
    """
    staged = []
    committed = []
    try:
        for path, content in artifacts:
            path.parent.mkdir(parents=True, exist_ok=True)
            previous = path.read_bytes() if path.exists() else None
            with NamedTemporaryFile(dir=path.parent, prefix=".sass-", delete=False) as stream:
                temporary = Path(stream.name)
                staged.append((path, temporary, previous))
                stream.write(content.encode("utf-8"))
                stream.flush()
        for path, temporary, previous in staged:
            temporary.replace(path)
            committed.append((path, previous))
    except OSError:
        try:
            for path, previous in reversed(committed):
                if previous is None:
                    path.unlink(missing_ok=True)
                else:
                    with NamedTemporaryFile(dir=path.parent, prefix=".sass-", delete=False) as stream:
                        rollback = Path(stream.name)
                        staged.append((path, rollback, None))
                        stream.write(previous)
                    rollback.replace(path)
        except OSError as error:
            raise DomainError("Report write and rollback failed; inspect both artifact hashes.") from error
        raise
    finally:
        for _, temporary, _ in staged:
            temporary.unlink(missing_ok=True)
