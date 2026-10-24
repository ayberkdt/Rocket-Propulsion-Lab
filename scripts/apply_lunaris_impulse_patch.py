"""Apply the reviewed, byte-bound native fix; preserve unrelated live changes."""

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    options = parser.parse_args()
    manifest_path = ROOT / "docs/assurance/lunaris_impulse_patch.json"
    manifest = json.loads(manifest_path.read_text())
    declared = Path(manifest["live_path"]).absolute()
    target = declared.resolve(strict=True)
    if target.parent != declared.parent.resolve() or target.name != declared.name:
        raise ValueError("Live source resolves outside its explicitly named target directory.")
    candidate = (ROOT / manifest["candidate_source"]).resolve(strict=True)
    if not candidate.is_relative_to(ROOT):
        raise ValueError("Candidate source must remain in this workspace.")
    after = candidate.read_bytes()
    before = target.read_bytes()
    if digest(after) != manifest["candidate_source_sha256"]:
        raise ValueError("Candidate bytes differ from the reviewed manifest.")
    if digest(before) not in {manifest["base_source_sha256"], manifest["candidate_source_sha256"]}:
        raise ValueError("Live source changed; rebase and reverify instead of overwriting it.")
    result = {
        "target": str(target),
        "before_sha256": digest(before),
        "candidate_sha256": digest(after),
        "status": "dry_run",
        "unrelated_files_modified": False,
    }
    if options.apply:
        if before != after:
            backup = ROOT / "output/assurance/native_lunaris_live_before_patch.py"
            if backup.exists() and backup.read_bytes() != before:
                raise ValueError("Existing backup differs; do not overwrite evidence.")
            backup.write_bytes(before)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(
                    dir=target.parent, prefix=".rocketprop-impulse-", suffix=".tmp", delete=False
                ) as stream:
                    temporary = Path(stream.name).resolve()
                    stream.write(after)
                # Both absolute move targets must be inside the explicit live
                # module directory, and the source must still match the review.
                if temporary.parent != target.parent or target.read_bytes() != before:
                    raise ValueError("Target directory or live source changed during staging.")
                os.replace(temporary, target)
                temporary = None
            finally:
                if temporary is not None and temporary.parent == target.parent:
                    temporary.unlink(missing_ok=True)
        if target.read_bytes() != after:
            raise ValueError("Live module does not match the reviewed candidate.")
        manifest["status"] = "applied_to_live_source_sha_verified"
        manifest["live_target_sha256_after"] = digest(after)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        result["status"] = "applied" if before != after else "already_applied"
        result["after_sha256"] = digest(target.read_bytes())
    output = ROOT / "output/assurance/native_lunaris_apply_result.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
