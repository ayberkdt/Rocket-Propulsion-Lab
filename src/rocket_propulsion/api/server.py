"""Dependency-free local HTTP server for the Rocket Propulsion Lab UI."""

from __future__ import annotations

import argparse
import json
import sys
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from rocket_propulsion import __version__
from rocket_propulsion.core.errors import (
    ConvergenceError,
    FeatureUnavailableError,
    InputError,
    RocketPropulsionError,
)

from .routes import calculation_metadata, dispatch_calculation, sidera_capabilities_route

WEB_ROOT = Path(__file__).resolve().parents[1] / "web"
STATIC_FILES = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/styles.css": ("styles.css", "text/css; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
}


class RocketPropulsionHandler(BaseHTTPRequestHandler):
    """Serve the single-page application and versioned calculation API."""

    server_version = f"RocketPropulsionLab/{__version__}"

    def _write_json(self, status: HTTPStatus, body: dict[str, Any]) -> None:
        payload = json.dumps(body, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def _write_error(self, status: HTTPStatus, error: RocketPropulsionError) -> None:
        """Write the stable machine-readable error envelope."""

        self._write_json(status, {"error": error.as_dict()})

    def do_GET(self) -> None:
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

        asset = STATIC_FILES.get(path)
        if asset is None:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        filename, content_type = asset
        try:
            payload = (WEB_ROOT / filename).read_bytes()
        except OSError:
            self.send_error(HTTPStatus.INTERNAL_SERVER_ERROR, "Web asset unavailable")
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_POST(self) -> None:
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
            except (json.JSONDecodeError, UnicodeDecodeError) as error:
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
        print(f"[{self.log_date_time_string()}] {message_format % args}")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for local server configuration."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1", help="Interface to bind.")
    parser.add_argument("--port", default=8765, type=int, help="TCP port to bind.")
    return parser


def main() -> None:
    """Run the local calculation server until interrupted."""

    if len(sys.argv) > 1 and sys.argv[1] == "burn":
        from rocket_propulsion.propulsion.burns.cli import main as burn_main

        raise SystemExit(burn_main(sys.argv[2:]))

    arguments = build_parser().parse_args()
    server = ThreadingHTTPServer((arguments.host, arguments.port), RocketPropulsionHandler)
    print(f"Rocket Propulsion Lab: http://{arguments.host}:{arguments.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()

