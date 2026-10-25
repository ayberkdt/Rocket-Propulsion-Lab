"""Dependency-free local HTTP server for the Rocket Propulsion Lab UI."""

from __future__ import annotations

import argparse
import json
import sys
import threading
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse
from urllib.request import urlopen

from rocket_propulsion import __version__
from rocket_propulsion.core.errors import (
    ConvergenceError,
    FeatureUnavailableError,
    InputError,
    RocketPropulsionError,
)

from .routes import calculation_metadata, dispatch_calculation, sidera_capabilities_route

WEB_ROOT = Path(__file__).resolve().parents[1] / "web"
# Only these file types are served from the web tree; anything else is 404.
CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".woff2": "font/woff2",
}
DEFAULT_PORT = 8765
FALLBACK_PORTS = (8766, 8767, 8768, 8769, 8770)


def resolve_web_asset(url_path: str) -> Path | None:
    """Map a decoded URL path to a servable file inside ``WEB_ROOT``.

    Parameters
    ----------
    url_path : str
        Percent-decoded URL path; ``""`` and ``"/"`` mean ``/index.html``.

    Returns
    -------
    pathlib.Path or None
        The resolved file, or ``None`` when the path contains a backslash or
        NUL, resolves outside ``WEB_ROOT``, has a suffix not listed in
        ``CONTENT_TYPES``, or is not an existing file.
    """

    if url_path in ("", "/"):
        url_path = "/index.html"
    if "\\" in url_path or "\x00" in url_path:
        return None
    candidate = (WEB_ROOT / url_path.lstrip("/")).resolve()
    root = WEB_ROOT.resolve()
    if candidate.suffix not in CONTENT_TYPES or not candidate.is_relative_to(root):
        return None
    return candidate if candidate.is_file() else None


class RocketPropulsionHandler(BaseHTTPRequestHandler):
    """Serve the single-page application and the versioned calculation API.

    ``GET`` serves ``/api/v1/health``, the Sidera capability snapshot and
    static web assets; ``POST`` dispatches JSON calculation requests. Error
    responses use the envelope ``{"error": {"code", "message", ...}}``.
    """

    server_version = f"RocketPropulsionLab/{__version__}"

    def _write_json(self, status: HTTPStatus, body: dict[str, Any]) -> None:
        """Send ``body`` as strict JSON (no NaN or Infinity) with ``no-store`` caching."""
        payload = json.dumps(body, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def _write_error(self, status: HTTPStatus, error: RocketPropulsionError) -> None:
        """Write ``{"error": error.as_dict()}`` with the given HTTP status."""

        self._write_json(status, {"error": error.as_dict()})

    def do_GET(self) -> None:
        """Answer health and capability requests, otherwise serve a web asset or 404."""
        path = urlparse(self.path).path
        if path == "/api/v1/health":
            self._write_json(
                HTTPStatus.OK,
                {
                    "status": "ok",
                    "version": __version__,
                    "meta": {
                        "service": "rocket-propulsion-lab",
                        "api_version": "v1",
                        "unit_system": "SI",
                    },
                },
            )
            return
        if path == "/api/v1/integrations/sidera/capabilities":
            self._write_json(
                HTTPStatus.OK,
                {
                    "data": sidera_capabilities_route({}),
                    "meta": calculation_metadata(path),
                },
            )
            return

        asset = resolve_web_asset(unquote(path))
        if asset is None:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            payload = asset.read_bytes()
        except OSError:
            self.send_error(HTTPStatus.INTERNAL_SERVER_ERROR, "Web asset unavailable")
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", CONTENT_TYPES[asset.suffix])
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-cache")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(payload)

    def do_POST(self) -> None:
        """Decode a JSON object body, dispatch it and write the result envelope.

        The body must be 1 byte to 1 MB of UTF-8 JSON whose root is an
        object. Unknown routes return 404, convergence failures 422,
        unavailable optional features 501, other library errors 400 and
        unexpected exceptions 500 with a generic message.
        """
        path = urlparse(self.path).path
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 1_000_000:
                raise InputError(
                    "Request body must contain JSON and be smaller than 1 MB.",
                    code="invalid_body_size",
                    field="body",
                )
            try:
                decoded = json.loads(self.rfile.read(length))
            except (ValueError, RecursionError) as error:  # incl. JSON and Unicode decode errors
                raise InputError(
                    "Request body must be valid UTF-8 JSON.",
                    code="invalid_json",
                    field="body",
                ) from error
            if not isinstance(decoded, dict):
                raise InputError(
                    "Request JSON must be an object.",
                    code="invalid_json_type",
                    field="body",
                )
            result = dispatch_calculation(path, decoded)
        except KeyError:
            self._write_error(
                HTTPStatus.NOT_FOUND,
                InputError("Unknown API route.", code="unknown_route", field="path"),
            )
            return
        except ConvergenceError as error:
            self._write_error(HTTPStatus.UNPROCESSABLE_ENTITY, error)
            return
        except FeatureUnavailableError as error:
            self._write_error(HTTPStatus.NOT_IMPLEMENTED, error)
            return
        except RocketPropulsionError as error:
            self._write_error(HTTPStatus.BAD_REQUEST, error)
            return
        except Exception:  # noqa: BLE001 - HTTP boundary must map unexpected faults
            self._write_error(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                RocketPropulsionError(
                    "Unexpected calculation error.", code="internal_error"
                ),
            )
            return
        self._write_json(
            HTTPStatus.OK,
            {"data": result, "meta": calculation_metadata(path)},
        )

    def log_message(self, message_format: str, *args: object) -> None:
        """Print one timestamped access-log line to standard output."""
        print(f"[{self.log_date_time_string()}] {message_format % args}")


class LabServer(ThreadingHTTPServer):
    """Threading HTTP server that refuses to share a port on Windows.

    ``HTTPServer`` enables ``SO_REUSEADDR``, which on Windows lets a second
    process bind a port that is already serving; the port fallback in
    :func:`main` would then never trigger. Handler threads are daemon threads.
    """

    allow_reuse_address = sys.platform != "win32"
    daemon_threads = True
    # The workbench opens with a burst of parallel requests (base solves plus chart
    # sweeps); the socketserver default backlog of 5 refuses some of them.
    request_queue_size = 128


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser (``--host``, ``--port`` and ``--open``)."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1", help="Interface to bind.")
    parser.add_argument(
        "--port", type=int, default=None,
        help=f"TCP port to bind (default {DEFAULT_PORT}, then {FALLBACK_PORTS[0]}–"
             f"{FALLBACK_PORTS[-1]} if it is busy).",
    )
    parser.add_argument(
        "--open", action="store_true",
        help="Open the workbench in the default browser once the server is ready.",
    )
    return parser


def running_lab_url(host: str, port: int) -> str | None:
    """Return the base URL if a Rocket Propulsion Lab server already answers there.

    Parameters
    ----------
    host : str
        Host name or address.
    port : int
        TCP port.

    Returns
    -------
    str or None
        ``"http://host:port"`` when ``/api/v1/health`` answers within one
        second and identifies the service as ``rocket-propulsion-lab``;
        otherwise ``None``.
    """

    url = f"http://{host}:{port}"
    try:
        with urlopen(f"{url}/api/v1/health", timeout=1.0) as response:
            body = json.load(response)
    except (OSError, ValueError):
        return None
    service = body.get("meta", {}).get("service") if isinstance(body, dict) else None
    return url if service == "rocket-propulsion-lab" else None


def main() -> None:
    """Run the local calculation server until interrupted.

    ``rocket-propulsion burn ...`` forwards to the burn command-line
    interface. Otherwise the default port and its fallbacks are tried in
    order; if a Lab server already answers on one, its URL is reported (and
    opened with ``--open``) instead of starting a second server.

    Raises
    ------
    SystemExit
        When no candidate port is free, or with the burn CLI's exit status.
    """

    if len(sys.argv) > 1 and sys.argv[1] == "burn":
        from rocket_propulsion.propulsion.burns.cli import main as burn_main

        raise SystemExit(burn_main(sys.argv[2:]))

    arguments = build_parser().parse_args()
    ports = [arguments.port] if arguments.port is not None else [DEFAULT_PORT, *FALLBACK_PORTS]
    server = None
    for port in ports:
        existing = running_lab_url(arguments.host, port)
        if existing:
            print(f"Rocket Propulsion Lab is already running: {existing}")
            if arguments.open:
                webbrowser.open(existing)
            return
        try:
            server = LabServer((arguments.host, port), RocketPropulsionHandler)
            break
        except OSError:
            print(f"Port {port} is busy.")
    if server is None:
        raise SystemExit(f"No free port among {', '.join(map(str, ports))}; use --port.")

    host, port = server.server_address[:2]
    url = f"http://{host}:{port}"
    print(f"Rocket Propulsion Lab: {url}")
    print("Press Ctrl+C (or close this window) to stop the server.")
    if arguments.open:
        threading.Timer(0.3, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()

