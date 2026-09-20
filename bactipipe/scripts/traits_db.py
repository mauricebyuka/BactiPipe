#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
traits_db.py
Version collection and lightweight DB checks for traits-detect.

Collects tool and database versions *only* for selected components,
so the report header lists exactly what was used in this run.
"""

from __future__ import annotations

import os
import re
import json
import subprocess
from typing import Dict, Tuple, Optional

# Optional centralized env launcher (recommended).
# env_cmd("amrfinder") returns a command prefix list to run the tool inside the right conda env.
try:
    from bactipipe.scripts.core import env_cmd  # type: ignore
except Exception:
    env_cmd = None


# ----------------------------- public API ----------------------------- #

def ensure_and_collect_versions(
    *,
    want_amr: bool,
    want_vf: bool,
    want_resfinder: bool = False,
    want_card: bool = False,
    want_vfinder: bool = False,
    want_vfdb: bool = False,
    want_ecoli_vf: bool = False,
    want_ncbi: bool = False,
    db_root: Optional[str] = None,
) -> Dict:
    """
    Probe versions for exactly the requested components.

    Returns dict:
      {
        "tools": {
          "amrfinder": "vX.Y.Z",
          "virulencefinder": "vA.B.C",
          "abricate": "v1.0.1",
        },
        "databases": {
          "amrfinder_db": "YYYY-MM-DD.N",
          "virulencefinder_db": "vA.B.C",
          "vfdb": "YYYY-MM-DD",
          "ecoli_vf": "YYYY-MM-DD",
          "resfinder_db": "YYYY-MM-DD",
          "card_db": "YYYY-MM-DD",
        }
      }
    """
    tools: Dict[str, str] = {}
    dbs: Dict[str, str] = {}

    # --- AMRFinderPlus (tool + DB) ---
    if want_amr:
        amr_tool_ver, amr_db_ver = _amrfinder_versions()
        if amr_tool_ver:
            tools["amrfinder"] = amr_tool_ver
        if amr_db_ver:
            dbs["amrfinder_db"] = amr_db_ver

    # --- VirulenceFinder (tool) ---
    if want_vf and want_vfinder:
        vfinder_ver = _virulencefinder_version()
        if vfinder_ver:
            tools["virulencefinder"] = vfinder_ver

    # --- VirulenceFinder DB (from VERSION file) ---
    if want_vf and want_vfinder:
        vdb_ver = _virulencefinder_db_version(db_root)
        if vdb_ver:
            dbs["virulencefinder_db"] = vdb_ver

    # --- ABRicate (tool) + DBs via `abricate --list` ---
    # Only probe if any abricate-backed DB was requested.
    if want_ncbi or want_resfinder or want_card or (want_vf and (want_vfdb or want_ecoli_vf)):
        ab_ver = _abricate_version()
        if ab_ver:
            tools["abricate"] = ab_ver

        # Parse abricate --list once and reuse
        list_map = _abricate_list_versions()
        if want_ncbi and list_map.get("ncbi"):
            dbs["ncbi"] = list_map["ncbi"]

        if want_vfdb:
            ver = list_map.get("vfdb")
            if ver:
                dbs["vfdb"] = ver

        if want_ecoli_vf:
            ver = list_map.get("ecoli_vf")
            if ver:
                dbs["ecoli_vf"] = ver

        if want_resfinder:
            ver = list_map.get("resfinder")
            if ver:
                # Keep key consistent with your header labeling
                dbs["resfinder_db"] = ver

        if want_card:
            ver = list_map.get("card")
            if ver:
                dbs["card_db"] = ver
    return {"tools": tools, "databases": dbs}


# ----------------------------- helpers ----------------------------- #

def _run_cmd(cmd_list: list[str]) -> Tuple[int, str, str]:
    """Run a command, return (rc, stdout, stderr) with text outputs."""
    try:
        p = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return p.returncode, p.stdout or "", p.stderr or ""
    except Exception as e:
        return 127, "", str(e)


def _env_launch(exe: str, *args: str) -> list[str]:
    """
    Build a command list, honoring env_cmd if available.
    """
    if env_cmd is not None:
        return env_cmd(exe) + list(args)
    return [exe, *args]


# ------------------ AMRFinderPlus versions ------------------ #

_AMR_TOOL_RE = re.compile(r"Software version:\s*([0-9][^\s]*)", re.I)
_AMR_DB_RE   = re.compile(r"Database version:\s*([0-9][^\s]*)", re.I)

def _amrfinder_versions() -> Tuple[Optional[str], Optional[str]]:
    """
    Run `amrfinder -V` and parse:
      Software version: X
      Database version: YYYY-MM-DD.N
    """
    cmd = _env_launch("amrfinder", "-V")
    rc, out, err = _run_cmd(cmd)
    txt = f"{out}\n{err}"
    tool_v = None
    db_v = None
    if rc == 0:
        m = _AMR_TOOL_RE.search(txt)
        if m:
            tool_v = _fmt_ver(m.group(1))
        m = _AMR_DB_RE.search(txt)
        if m:
            db_v = m.group(1).strip()
    return tool_v, db_v


# ------------------ VirulenceFinder tool version ------------------ #

def _virulencefinder_version() -> Optional[str]:
    """
    Best-effort discovery of VirulenceFinder version.
    Query the configured viramr environment. Never report a hard-coded version
    when the installed version cannot be verified.
    """
    env_name = os.environ.get("VIRAMR_ENV", "viramr")
    rc, out, _ = _run_cmd(["conda", "list", "-n", env_name, "virulencefinder", "--json"])
    if rc != 0 or not out.strip():
        return None
    try:
        records = json.loads(out)
    except (TypeError, ValueError):
        return None
    for record in records:
        if str(record.get("name", "")).casefold() == "virulencefinder":
            version = str(record.get("version", "")).strip()
            return _fmt_ver(version) if version else None
    return None


def _virulencefinder_db_version(db_root: Optional[str]) -> Optional[str]:
    """
    Read first line of <db_root>/virulencefinder_db/VERSION.
    """
    if not db_root:
        return None
    path = os.path.join(db_root, "virulencefinder_db", "VERSION")
    try:
        with open(path) as f:
            line = f.readline().strip()
            return line if line else None
    except Exception:
        return None


# ------------------ ABRicate tool & DB versions ------------------ #

_AB_VERSION_RE = re.compile(r"\babricate\s+([0-9][^\s]*)", re.I)

def _abricate_version() -> Optional[str]:
    """
    Parse `abricate --version` first line like 'abricate 1.0.1'.
    """
    cmd = _env_launch("abricate", "--version")
    rc, out, err = _run_cmd(cmd)
    txt = (out or err).strip()
    m = _AB_VERSION_RE.search(txt.splitlines()[0] if txt else "")
    if m:
        return _fmt_ver(m.group(1).strip())
    return None


def _abricate_list_versions() -> Dict[str, str]:
    """
    Run `abricate --list` in the abricate environment and parse its table.
    Returns {db_name_lower: date_string} using the DATE column for each DB.
    """
    out_map: Dict[str, str] = {}

    # Build command in the right env (no shell, no piping)
    if env_cmd is not None:
        cmd = env_cmd("abricate") + ["--list"]
    else:
        cmd = ["abricate", "--list"]

    rc, out, err = _run_cmd(cmd)
    if rc != 0:
        # If abricate prints to stderr in some envs, try stderr text too
        text = out or err
        if not text:
            return out_map
    else:
        text = out

    lines = (text or "").strip().splitlines()
    if not lines:
        return out_map

    # Parse header to find the DATE column index (robust to extra whitespace/tabs)
    header = lines[0].strip().split()
    # Expected: ["DATABASE","SEQUENCES","DBTYPE","DATE"], but be defensive
    try:
        date_idx = header.index("DATE")
    except ValueError:
        # Fallback: assume last column is date
        date_idx = -1

    # Iterate rows, skip the header
    for i, line in enumerate(lines):
        parts = line.strip().split()
        if not parts:
            continue
        # Skip header row by name match
        if i == 0 and parts[0].upper() == "DATABASE":
            continue
        # Defensive: need at least 2 columns (DB + DATE)
        if len(parts) < 2:
            continue

        db_raw = parts[0]
        db = db_raw.lstrip("*").strip().lower()  # strip active marker, normalize

        # DATE column (use located index; fallback to last token)
        date_str = parts[date_idx] if (date_idx != -1 and date_idx < len(parts)) else parts[-1]
        date_str = date_str.strip()

        if db and date_str:
            out_map[db] = date_str

    return out_map

# ------------------ formatting helpers ------------------ #

def _fmt_ver(v: str) -> str:
    """Normalize versions to 'vX.Y.Z' style without double 'v'."""
    v = (v or "").strip()
    if not v:
        return v
    return v if v.startswith("v") else f"v{v}"
