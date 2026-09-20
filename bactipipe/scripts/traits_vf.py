#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
traits_vf.py
Virulence detection from assembled genomes.

Sources:
  - VirulenceFinder (DTU python tool; requires virulencefinder_db)
  - ABRicate VFDB
  - ABRicate ecoli_vf

Rules:
  - For organism == "Escherichia": run all three sources.
  - Otherwise: run only the sources enabled by flags from run_traits.

Mapping:
  The virulence map file has 4 columns:
    gene    display_name    functional_category    notes
  Only genes present in the map are kept. Others are dropped.
  For Escherichia:
    - default map: dx_vf_category_map.tsv
    - if extended_vf=True: all_vf_category_map.tsv
  For all other organisms:
    - map: all_vf_category_map.tsv
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
from bactipipe.scripts.traits_organisms import virulence_databases

# Optional centralized env launcher (recommended)
try:
    from bactipipe.scripts.core import env_cmd  # returns list prefix to run tool in the right env
except Exception:
    env_cmd = None

# Logger type hint: callable(msg: str, level?: str)
LogFn = callable


# ----------------------------- public entry ----------------------------- #

def run_vf_for_sample(
    *,
    sample: str,
    fasta: str,
    out_root: str,
    min_id: float,
    min_cov: float,
    enable_vfinder: bool,
    enable_vfdb: bool,
    enable_ecoli_vf: bool,
    threads: int,
    logger_fn,
    db_root: str,
    organism: str = "",
    extended_vf: bool = False,
) -> List[Dict]:
    """
    Run virulence detection for a single sample and return normalized rows.

    Each returned row has:
      {
        "gene":         <display_name>,           # from map
        "base_gene":    <raw symbol as called>,   # original symbol (e.g., espB)
        "category":     <functional_category>,
        "note":         <notes>,
        "source":       "VirulenceFinder" | "VFDB" | "ecoli_vf" | "VirulenceFinder,VFDB" (after de-dupe)
      }
    """

    # --- Determine which sources to run (Escherichia rule)
    if organism.casefold() == "escherichia":
        enable_vfinder = True
        enable_vfdb = True
        enable_ecoli_vf = True

    vdb = os.path.join(db_root or "", "virulencefinder_db") if db_root else ""
    logger_fn(
        f"[{sample}] VF plan | organism={organism or '-'} | "
        f"VirulenceFinder={enable_vfinder} (db={vdb or 'n/a'}) | "
        f"VFDB={enable_vfdb} | ecoli_vf={enable_ecoli_vf}"
    )

    # --- Load the map
    map_path = _select_vf_map_path(organism, extended_vf)
    vf_map = _read_vf_map(map_path)
    logger_fn(f"[{sample}] Using virulence map: {map_path.name} ({len(vf_map)} entries)")

    rows: List[Dict] = []

    # --- VirulenceFinder
    if enable_vfinder:
        if not db_root:
            raise RuntimeError(f"[{sample}] VirulenceFinder database root is not configured")
        elif not os.path.isdir(vdb):
            raise RuntimeError(f"[{sample}] VirulenceFinder database not found: {vdb}")
        else:
            databases = virulence_databases(organism)
            if not databases:
                raise ValueError(f"No VirulenceFinder database configured for {organism}")
            vf_tsv = _run_virulencefinder(sample, fasta, out_root, db_root, logger_fn,
                min_id=min_id / 100, min_cov=min_cov / 100, databases=databases)
            if vf_tsv:
                rows.extend(_parse_vfinder_results(vf_tsv, min_id, min_cov, vf_map, logger_fn))
            else:
                raise RuntimeError(f"[{sample}] VirulenceFinder failed or produced no results table; see its logs")
    else:
        logger_fn(f"[{sample}] VirulenceFinder disabled by flags.")

    # --- ABRicate VFDB
    if enable_vfdb:
        vfdb_tsv = _run_abricate(sample, fasta, out_root, db="vfdb", threads=threads, logger_fn=logger_fn, min_id=min_id, min_cov=min_cov)
        if vfdb_tsv:
            rows.extend(_parse_abricate(vfdb_tsv, source_label="VFDB", min_id=min_id, min_cov=min_cov, vf_map=vf_map, logger_fn=logger_fn))
        else:
            raise RuntimeError(f"[{sample}] ABRicate/VFDB failed or produced no results table; see its logs")
    else:
        logger_fn(f"[{sample}] ABRicate(VFDB) disabled by flags.")

    # --- ABRicate ecoli_vf
    if enable_ecoli_vf:
        ecoli_tsv = _run_abricate(sample, fasta, out_root, db="ecoli_vf", threads=threads, logger_fn=logger_fn, min_id=min_id, min_cov=min_cov)
        if ecoli_tsv:
            rows.extend(_parse_abricate(ecoli_tsv, source_label="ecoli_vf", min_id=min_id, min_cov=min_cov, vf_map=vf_map, logger_fn=logger_fn))
        else:
            raise RuntimeError(f"[{sample}] ABRicate/ecoli_vf failed or produced no results table; see its logs")
    else:
        logger_fn(f"[{sample}] ABRicate(ecoli_vf) disabled by flags.")

    # Normalize & de-duplicate across sources (prefer VirulenceFinder)
    rows = _dedupe_rows(rows, logger_fn)
    logger_fn(f"[{sample}] Virulence merged rows: {len(rows)}")
    return rows


# ------------------------- map loading & selection ------------------------- #

def _select_vf_map_path(organism: str, extended_vf: bool) -> Path:
    """
    Pick the map file according to organism and extended flag.
    """
    try:
        from importlib.resources import files as resource_files
    except Exception:
        from importlib_resources import files as resource_files  # type: ignore

    data_pkg = resource_files("bactipipe.data")
    if organism.casefold() == "escherichia":
        fname = "all_vf_category_map.tsv" if extended_vf else "dx_vf_category_map.tsv"
    else:
        fname = "all_vf_category_map.tsv"
    return data_pkg.joinpath(fname)


def _read_vf_map(path: Path) -> Dict[str, Dict[str, str]]:
    """
    Read 4-column virulence map:
      gene    display_name    functional_category    notes
    Returns dict keyed by raw 'gene' -> {display_name, category, note}
    """
    out: Dict[str, Dict[str, str]] = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            g = (rec.get("gene") or "").strip()
            disp = (rec.get("display_name") or "").strip()
            cat = (rec.get("functional_category") or "").strip()
            note = (rec.get("notes") or "").strip()
            if not g or not disp or not cat:
                continue
            out[g] = {"display_name": disp, "category": cat, "note": note}
    return out


# --------------------------- tool runners (subprocess) --------------------------- #

def _run_virulencefinder(
    sample: str,
    fasta: str,
    out_root: str,
    db_root: str,
    logger_fn,
    *,
    min_id: float = 0.9,
    min_cov: float = 0.6,
    databases: str = "all",
) -> str | None:
    """
    Run DTU VirulenceFinder via module invocation:
        python -m virulencefinder ...

    Writes directly into:
        <out_root>/raw/virulencefinder/<sample>/

    Returns the path to results_tab.tsv on success, else None.
    """
    import os
    import subprocess

    # Validate inputs
    if not os.path.isfile(fasta):
        logger_fn(f"[{sample}] VirulenceFinder SKIP: FASTA not found: {fasta}")
        return None
    if not db_root:
        logger_fn(f"[{sample}] VirulenceFinder SKIP: --db-dir / $BACTIPIPE_DB_DIR not set.")
        return None

    # VirulenceFinder database root
    vdb_root = os.path.join(db_root, "virulencefinder_db")
    if not os.path.isdir(vdb_root):
        logger_fn(f"[{sample}] VirulenceFinder SKIP: DB root not found: {vdb_root}")
        return None

    # Final output dir
    final_dir = os.path.join(out_root, "raw", "virulencefinder", sample)
    os.makedirs(final_dir, exist_ok=True)

    # Temp dir inside final_dir (keeps all VF temp artifacts contained)
    tmp_dir = os.path.join(final_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # Build command (NOTE: -d expects DB NAMES, not a path)
    cmd = env_cmd("virulencefinder") + [
        "-ifa", fasta,
        "-o", final_dir,
        "-tmp", tmp_dir,
        "-p", vdb_root,
        "-d", databases,      # e.g. "all" or "virulence_ecoli"
        "-t", str(min_id),
        "-l", str(min_cov),
        "-x",
    ]

    logger_fn(
        f"[{sample}] VirulenceFinder TRY | DB_ROOT={vdb_root} | DBs={databases} | "
        f"OUT={final_dir} | TMP={tmp_dir}"
    )
    logger_fn(f"[{sample}] VirulenceFinder CMD: {' '.join(cmd)}")

    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        logger_fn(f"[{sample}] VirulenceFinder FAILED to start: {e}")
        return None

    # Save logs regardless of success (useful for debugging)
    try:
        with open(os.path.join(final_dir, "virulencefinder.stdout.txt"), "w") as f:
            f.write(p.stdout or "")
        with open(os.path.join(final_dir, "virulencefinder.stderr.txt"), "w") as f:
            f.write(p.stderr or "")
    except Exception as e:
        logger_fn(f"[{sample}] VirulenceFinder WARN: failed to save stdout/stderr logs: {e}")

    if p.returncode != 0:
        logger_fn(
            f"[{sample}] VirulenceFinder FAILED (exit {p.returncode}).\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
        return None

    res_tsv = os.path.join(final_dir, "results_tab.tsv")
    if not os.path.isfile(res_tsv):
        logger_fn(
            f"[{sample}] VirulenceFinder OK exit but missing results_tab.tsv.\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
        return None

    logger_fn(f"[{sample}] VirulenceFinder OK → {res_tsv}")
    return res_tsv

def _run_abricate(
    sample: str,
    fasta: str,
    out_root: str,
    *,
    db: str,
    threads: int,
    logger_fn,
    min_id: float = 90,
    min_cov: float = 60,
) -> str | None:
    """
    Run ABRicate for a specific DB (vfdb or ecoli_vf). Returns path to TSV on success.
    """
    if not os.path.isfile(fasta):
        logger_fn(f"[{sample}] ABRicate({db}) SKIP: FASTA not found: {fasta}")
        return None

    tmp_dir = os.path.join(out_root, ".tmp", f"abricate_{db}", sample)
    os.makedirs(tmp_dir, exist_ok=True)
    out_tsv = os.path.join(tmp_dir, f"{sample}.{db}.tsv")

    # Build command; ABRicate accepts --db and --threads
    if env_cmd is not None:
        cmd = env_cmd("abricate") + ["--db", db]
    else:
        cmd = ["abricate", "--db", db]
    if threads and threads > 1:
        cmd += ["--threads", str(threads)]
    cmd += ["--minid", str(min_id), "--mincov", str(min_cov)]
    cmd += [fasta]

    logger_fn(f"[{sample}] ABRicate TRY | db={db} | TMP={tmp_dir}")
    logger_fn(f"[{sample}] ABRicate CMD: {' '.join(cmd)}")

    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        logger_fn(f"[{sample}] ABRicate({db}) FAILED to start: {e}")
        return None

    if p.returncode != 0:
        logger_fn(f"[{sample}] ABRicate({db}) FAILED (exit {p.returncode}).\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
        return None

    # ABRicate prints to STDOUT unless --csv/--tsv redirected; write capture to out_tsv
    try:
        with open(out_tsv, "w") as f:
            f.write(p.stdout or "")
    except Exception as e:
        logger_fn(f"[{sample}] ABRicate({db}) FAILED to write TSV: {e}")
        return None

    if not os.path.isfile(out_tsv) or os.path.getsize(out_tsv) == 0:
        logger_fn(f"[{sample}] ABRicate({db}) produced empty TSV.")
        return None

    # move final TSV to raw/abricate_<db>/<sample>.tsv
    final_dir = os.path.join(out_root, "raw", f"abricate_{db}")
    os.makedirs(final_dir, exist_ok=True)
    final_tsv = os.path.join(final_dir, f"{sample}.tsv")
    try:
        shutil.move(out_tsv, final_tsv)
        with open(os.path.join(final_dir, f"{sample}.stdout.txt"), "w") as f:
            f.write(p.stdout or "")
        with open(os.path.join(final_dir, f"{sample}.stderr.txt"), "w") as f:
            f.write(p.stderr or "")
    except Exception as e:
        logger_fn(f"[{sample}] ABRicate({db}) FAILED to move outputs: {e}")
        return None

    logger_fn(f"[{sample}] ABRicate({db}) OK → {final_tsv}")
    return final_tsv


# ------------------------------ parsers ------------------------------ #

def _parse_vfinder_results(
    tsv_path: str,
    min_id: float,
    min_cov: float,
    vf_map: Dict[str, Dict[str, str]],
    logger_fn,
) -> List[Dict]:
    """
    Parse VirulenceFinder extended results_tab.tsv
    Columns (as seen in your example):
      Database, Virulence factor, Identity, Query / Template length,
      Contig, Position in contig, Protein function, Accession number
    """
    rows: List[Dict] = []
    def _cov_pct(qt: str) -> float:
        # "97 / 117" -> 97/117 * 100
        try:
            q, t = qt.replace(",", "").split("/")
            q = int(str(q).strip())
            t = int(str(t).strip())
            if t <= 0:
                return 0.0
            return 100.0 * (q / t)
        except Exception:
            return 0.0

    with open(tsv_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            raw_gene = (rec.get("Virulence factor") or "").strip()
            if not raw_gene:
                continue
            ident = _safe_float(rec.get("Identity"))
            cov = _cov_pct(rec.get("Query / Template length") or "")
            if ident < min_id or cov < min_cov:
                continue
            map_rec = vf_map.get(raw_gene)
            if not map_rec:
                # drop genes not present in the map
                continue
            rows.append({
                "gene": map_rec["display_name"],
                "base_gene": raw_gene,
                "category": map_rec["category"],
                "note": map_rec["note"],
                "source": "VirulenceFinder",
                "identity": ident,
                "coverage": cov,
            })
    logger_fn(f"[VF parse] {Path(tsv_path).name}: kept {len(rows)} rows after thresholds/map filter")
    return rows


def _parse_abricate(
    tsv_path: str,
    source_label: str,
    min_id: float,
    min_cov: float,
    vf_map: Dict[str, Dict[str, str]],
    logger_fn,
) -> List[Dict]:
    """
    Parse ABRicate TSV output captured from STDOUT.
    ABRicate columns typically include:
      SEQUENCE, START, END, GENE, COVERAGE, COVERAGE_MAP, GAPS, %COVERAGE, %IDENTITY, DATABASE, ...
    We rely on 'GENE', '%COVERAGE', '%IDENTITY', 'DATABASE'.
    """
    rows: List[Dict] = []
    with open(tsv_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        # normalize headers
        headers = [h.strip() for h in (r.fieldnames or [])]
        cov_key = "%COVERAGE" if "%COVERAGE" in headers else "COVERAGE"
        id_key = "%IDENTITY" if "%IDENTITY" in headers else "IDENTITY"
        for rec in r:
            raw_gene = (rec.get("GENE") or "").strip()
            if not raw_gene:
                continue
            ident = _safe_float(rec.get(id_key))
            cov = _safe_float(rec.get(cov_key))
            if ident < min_id or cov < min_cov:
                continue
            map_rec = vf_map.get(raw_gene)
            if not map_rec:
                # drop genes not present in the map
                continue
            rows.append({
                "gene": map_rec["display_name"],
                "base_gene": raw_gene,
                "category": map_rec["category"],
                "note": map_rec["note"],
                "source": source_label,
                "identity": ident,
                "coverage": cov,
            })
    logger_fn(f"[ABRicate parse:{source_label}] {Path(tsv_path).name}: kept {len(rows)} rows after thresholds/map filter")
    return rows


# --------------------------- de-duplication --------------------------- #

def _dedupe_rows(rows: List[Dict], logger_fn) -> List[Dict]:
    """
    De-duplicate rows by display name (gene), preferring VirulenceFinder over ABRicate sources.
    Merge 'source' labels if the same gene was found by multiple sources.
    """
    if not rows:
        return []

    by_gene: Dict[str, Dict] = {}
    for r in rows:
        disp = (r.get("gene") or "").strip()
        if not disp:
            continue
        src = r.get("source", "")
        if disp not in by_gene:
            # copy
            by_gene[disp] = dict(r)
            # ensure 'source' is a set for merging
            by_gene[disp]["source"] = {src} if src else set()
        else:
            # merge source
            cur = by_gene[disp]
            if src:
                if isinstance(cur["source"], set):
                    cur["source"].add(src)
                else:
                    cur["source"] = {cur["source"], src}
            # prefer VirulenceFinder's annotation if present
            if src == "VirulenceFinder":
                by_gene[disp] = dict(r)
                by_gene[disp]["source"] = set(cur["source"])  # keep merged sources

    # finalize
    out: List[Dict] = []
    for disp, rec in by_gene.items():
        src_val = rec.get("source")
        if isinstance(src_val, set):
            rec["source"] = ",".join(sorted(src_val))
        out.append(rec)

    out.sort(key=lambda x: (x.get("category", ""), x.get("gene", "").lower()))
    logger_fn(f"[de-dup] {len(rows)} → {len(out)} by display name")
    return out


# ------------------------------- utils ------------------------------- #

def _safe_float(x) -> float:
    try:
        if x is None:
            return 0.0
        s = str(x).strip().replace("%", "")
        return float(s)
    except Exception:
        return 0.0
