"""Web workbench shell: static serving, vendored Sunumatik presets and launcher logic."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import threading
import unittest
from http.server import ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import urlopen

from rocket_propulsion.api import server as server_module
from rocket_propulsion.api.server import WEB_ROOT, RocketPropulsionHandler, resolve_web_asset

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VENDOR = WEB_ROOT / "vendor" / "sunumatik"


class QuietHandler(RocketPropulsionHandler):
    def log_message(self, message_format: str, *args: object) -> None:
        pass


class StaticServingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), QuietHandler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base_url = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls) -> None:
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=2.0)

    def get(self, path: str) -> tuple[int, str, bytes]:
        try:
            with urlopen(self.base_url + path, timeout=2.0) as response:
                return response.status, response.headers["Content-Type"], response.read()
        except HTTPError as error:
            return error.code, error.headers["Content-Type"], b""

    def test_shell_and_vendored_assets_are_served_with_their_media_types(self) -> None:
        expected = {
            "/": "text/html",
            "/styles.css?v=3": "text/css",
            "/visuals.js": "text/javascript",
            "/app.js": "text/javascript",
            "/shell.mjs": "text/javascript",
            "/vendor/sunumatik/pano-odak.mjs": "text/javascript",
            "/vendor/sunumatik/palettes.json": "application/json",
            "/vendor/sunumatik/icons.svg": "image/svg+xml",
            "/vendor/sunumatik/typefaces/source-sans-3.woff2": "font/woff2",
        }
        for path, media_type in expected.items():
            with self.subTest(path=path):
                status, content_type, body = self.get(path)
                self.assertEqual(status, 200)
                self.assertTrue(content_type.startswith(media_type))
                self.assertTrue(body)

    def test_paths_outside_the_web_root_or_unlisted_types_are_not_served(self) -> None:
        for path in ("/../api/server.py", "/..%2f..%2fapi%2fserver.py", "/vendor/../../api/server.py",
                     "/vendor/sunumatik/typefaces/source-sans-3.LICENSE", "/README.md",
                     "/%5c..%5capi%5cserver.py", "/missing.css"):
            with self.subTest(path=path):
                self.assertEqual(self.get(path)[0], 404)

    def test_resolver_never_escapes_the_web_root(self) -> None:
        self.assertEqual(resolve_web_asset("/"), WEB_ROOT.resolve() / "index.html")
        for path in ("/../api/server.py", "/vendor/../../../pyproject.toml", "/a\\..\\b.css",
                     "/index.html\x00.css"):
            with self.subTest(path=path):
                self.assertIsNone(resolve_web_asset(path))


class VendoredPresetTests(unittest.TestCase):
    def test_vendored_files_match_their_recorded_hashes(self) -> None:
        provenance = json.loads((VENDOR / "PROVENANCE.json").read_text(encoding="utf-8"))
        self.assertEqual(provenance["repository"], "https://github.com/ayberkdt/sunumatik")
        self.assertRegex(provenance["commit"], r"^[0-9a-f]{40}$")
        for record in provenance["files"]:
            with self.subTest(path=record["path"]):
                digest = hashlib.sha256((VENDOR / record["path"]).read_bytes()).hexdigest()
                self.assertEqual(digest, record["sha256"])

    def test_check_mode_of_the_sync_script_passes(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "sync_sunumatik_presets.py"), "--check"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_every_local_reference_in_the_shell_exists(self) -> None:
        html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
        references = set(re.findall(r'(?:href|src)="(/[^"#?]+)', html))
        script = (WEB_ROOT / "shell.mjs").read_text(encoding="utf-8")
        references |= set(re.findall(r'"(/vendor/[^"]+)"', script))
        self.assertTrue(references)
        for reference in sorted(references):
            with self.subTest(reference=reference):
                self.assertIsNotNone(resolve_web_asset(reference), reference)

    def test_every_icon_used_by_the_shell_is_in_the_sprite(self) -> None:
        html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
        sprite = (VENDOR / "icons.svg").read_text(encoding="utf-8")
        for icon in sorted(set(re.findall(r'icons\.svg#([\w-]+)', html))):
            with self.subTest(icon=icon):
                self.assertIn(f'<symbol id="{icon}"', sprite)

    def test_default_palettes_exist_with_six_data_series(self) -> None:
        css = (VENDOR / "palette-library.css").read_text(encoding="utf-8")
        palettes = {item["id"]: item for item in
                    json.loads((VENDOR / "palettes.json").read_text(encoding="utf-8"))["palettes"]}
        for palette, dark in (("cerulean-vermilion", False), ("midnight-tangerine", True),
                              ("wada-salvia-indigo", False), ("wada-salvia-indigo-dark", True)):
            with self.subTest(palette=palette):
                self.assertEqual(palettes[palette]["dark"], dark)
                block = css.split(f'[data-palette="{palette}"] {{', 1)[1].split("}", 1)[0]
                self.assertIn("--color-data-6", block)


class LauncherTests(unittest.TestCase):
    def test_launcher_runs_the_server_module_with_browser_opening(self) -> None:
        launcher = (PROJECT_ROOT / "app.cmd").read_bytes()
        self.assertIn(b"\r\n", launcher, "cmd.exe needs CRLF line endings")
        self.assertIn(b"-m rocket_propulsion.api.server --open", launcher)
        self.assertIn(b'set "PYTHONPATH=%~dp0src"', launcher)

    def test_busy_default_port_falls_back_and_reuse_is_disabled_on_windows(self) -> None:
        if sys.platform == "win32":
            self.assertFalse(server_module.LabServer.allow_reuse_address)
        blocker = server_module.LabServer(("127.0.0.1", 0), QuietHandler)
        busy_port = blocker.server_address[1]
        started = {}

        def fake_serve(self) -> None:
            started["port"] = self.server_address[1]
            raise KeyboardInterrupt

        try:
            with patch.object(server_module, "DEFAULT_PORT", busy_port), \
                 patch.object(server_module, "FALLBACK_PORTS", (0,)), \
                 patch.object(server_module, "running_lab_url", return_value=None), \
                 patch.object(server_module.LabServer, "serve_forever", fake_serve), \
                 patch.object(sys, "argv", ["rocket-propulsion"]), \
                 patch("builtins.print"):
                server_module.main()
        finally:
            blocker.server_close()
        self.assertNotEqual(started["port"], busy_port)

    def test_an_already_running_lab_is_reused_instead_of_started_twice(self) -> None:
        with patch.object(server_module, "running_lab_url", return_value="http://127.0.0.1:8765"), \
             patch.object(server_module.webbrowser, "open") as opened, \
             patch.object(server_module, "LabServer") as server_class, \
             patch.object(sys, "argv", ["rocket-propulsion", "--open"]), \
             patch("builtins.print"):
            server_module.main()
        opened.assert_called_once_with("http://127.0.0.1:8765")
        server_class.assert_not_called()


if __name__ == "__main__":
    unittest.main()
