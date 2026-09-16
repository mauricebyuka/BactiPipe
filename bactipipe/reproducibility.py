"""Application source and database identities used in retained outputs."""

from __future__ import annotations

import json
import os
import hashlib
from pathlib import Path
from typing import Any

from bactipipe.__version__ import __version__


def source_identity() -> dict[str, Any]:
    return _source_identity("bactipipe", __version__)


def database_identity(database_root: str | Path | None = None) -> dict[str, Any]:
    root = Path(database_root or os.environ.get("BACTIPIPE_DB_DIR", "/opt/bactipipe/data"))
    return _read_database_identity(root / "database_identity.json", "Unversioned BactiPipe database")


def database_identity_label(database_root: str | Path | None = None) -> str:
    identity = database_identity(database_root)
    return str(identity.get("version") or identity.get("name") or "Unversioned")


def database_log_record(species_database_root: str | Path | None = None) -> str:
    """Compact report headers must not discard detailed database provenance."""
    identities = {"reference_bundle": database_identity()}
    if species_database_root:
        identities["species_identification"] = database_identity(species_database_root)
    return "Database provenance (bundle inventory, not a list of invoked components): " + json.dumps(
        identities, sort_keys=True
    )


def _source_identity(package_name: str, version: str) -> dict[str, Any]:
    payload: dict[str, Any] = {"version": version}
    record = _source_lock_record(package_name)
    if record.get("revision"):
        payload["source_revision"] = record["revision"]
    if record.get("tree_sha256"):
        payload["source_tree_sha256"] = record["tree_sha256"]
    if record.get("repository"):
        payload["repository"] = record["repository"]
    payload.setdefault("source_tree_sha256", _tree_digest(Path(__file__).resolve().parent))
    return payload


def _source_lock_record(package_name: str) -> dict[str, Any]:
    candidates = [
        Path(os.environ["PIPELINE_SOURCE_LOCK"])
        if os.environ.get("PIPELINE_SOURCE_LOCK") else None,
        Path(__file__).resolve().parents[1] / "pipeline-sources.lock.json",
        Path("/app/pipeline-sources.lock.json"),
    ]
    for lock_path in candidates:
        if lock_path is None:
            continue
        try:
            lock = json.loads(lock_path.read_text(encoding="utf-8"))
            record = lock.get("packages", {}).get(package_name, {})
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            continue
        if isinstance(record, dict):
            return record
    return {}


def _tree_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part in {"__pycache__", ".pytest_cache", ".DS_Store"} for part in relative.parts):
            continue
        if not path.is_file() or path.suffix in {".pyc", ".pyo"}:
            continue
        payload = path.read_bytes()
        digest.update(relative.as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(len(payload)).encode("ascii"))
        digest.update(b"\0")
        digest.update(payload)
        digest.update(b"\0")
    return digest.hexdigest()


def _read_database_identity(path: Path, fallback_name: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return {"name": fallback_name, "version": "Unversioned", "identity_file": str(path)}
    if not isinstance(payload, dict):
        return {"name": fallback_name, "version": "Unversioned", "identity_file": str(path)}
    return {**payload, "identity_file": str(path)}
