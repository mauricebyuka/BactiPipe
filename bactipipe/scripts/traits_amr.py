# bactipipe/scripts/traits_amr.py
from __future__ import annotations
import os, re, csv, subprocess
from typing import List, Dict, Optional

# --- env helpers --------------------------------------------------------------

def _env_prefix(tool: str, *, viramr_env: Optional[str], abricate_env: Optional[str]) -> List[str]:
    """
    Build a direct argv conda-run prefix for a given tool.
    - amrfinder -> viramr_env
    - abricate  -> abricate_env
    """
    if tool == "amrfinder" and viramr_env:
        return ["conda", "run", "-n", viramr_env, "amrfinder"]
    if tool == "abricate" and abricate_env:
        return ["conda", "run", "-n", abricate_env, "abricate"]
    return [tool]

def _run(cmd: List[str], *, env: Optional[dict] = None, logger_fn=None) -> None:
    """
    Run a direct argv command.
    """
    if logger_fn:
        logger_fn(f"[CMD] {' '.join(cmd)}")
    cp = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if logger_fn and cp.stdout:
        logger_fn(f"[STDOUT]\n{cp.stdout.strip()}")
    if cp.returncode != 0:
        if logger_fn and cp.stderr:
            logger_fn(f"[STDERR]\n{cp.stderr.strip()}", level="Warn")
        raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)

# --- AMRFinder parser ---------------------------------------------------------

_AMR_HDR_NORMAL = {
    "Contig id": ["Contig id", "ContigID", "Contig"],
    "Start": ["Start"],
    "Stop": ["Stop", "End"],
    "Strand": ["Strand"],
    "Element symbol": ["Element symbol", "Gene symbol", "Symbol"],
    "Element name": ["Element name", "Gene name", "Name"],
    "Scope": ["Scope"],
    "Type": ["Type"],
    "Subtype": ["Subtype"],
    "Class": ["Class"],
    "Subclass": ["Subclass"],
    "Method": ["Method"],
    "Target length": ["Target length", "TargetLen"],
    "Reference sequence length": ["Reference sequence length", "RefLen"],
    "% Coverage of reference": ["% Coverage of reference", "Coverage", "PercCoverageOfReference"],
    "% Identity to reference": ["% Identity to reference", "Identity", "PercIdentityToReference"],
    "Alignment length": ["Alignment length", "AlignLen"],
    "Closest reference accession": ["Closest reference accession", "RefAccession"],
    "Closest reference name": ["Closest reference name", "RefName"],
    "HMM accession": ["HMM accession", "HMMAccession"],
    "HMM description": ["HMM description", "HMMDescription"],
}

_MUT_RX = re.compile(r"^[A-Za-z0-9]+_[A-Z]\d+[A-Z]$")  # e.g., gyrA_S83F

def _norm_header_map(headers: List[str]) -> Dict[str, str]:
    hset = {h.strip(): h for h in headers}
    res: Dict[str, str] = {}
    for norm, variants in _AMR_HDR_NORMAL.items():
        for v in variants:
            if v in hset:
                res[norm] = hset[v]
                break
    return res

def _get_field(rec: List[str], headers: List[str], hmap: Dict[str, str], norm: str) -> str:
    src = hmap.get(norm)
    if not src:
        return ""
    try:
        return rec[headers.index(src)].strip()
    except Exception:
        return ""

def parse_amrfinder_tsv(tsv_path: str, *, sample: str, min_id: float, min_cov: float) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not os.path.exists(tsv_path):
        return rows
    with open(tsv_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        headers = next(reader, [])
        if not headers:
            return rows
        hmap = _norm_header_map(headers)
        for rec in reader:
            elem_symbol = _get_field(rec, headers, hmap, "Element symbol")
            if not elem_symbol:
                continue
            subtype = _get_field(rec, headers, hmap, "Subtype")
            rtype = _get_field(rec, headers, hmap, "Type")
            amr_class = _get_field(rec, headers, hmap, "Class")
            amr_subclass = _get_field(rec, headers, hmap, "Subclass")
            cov = _get_field(rec, headers, hmap, "% Coverage of reference")
            ident = _get_field(rec, headers, hmap, "% Identity to reference")

            # apply thresholds if present (values like "99.8", not "%")
            try:
                if ident and float(ident) < float(min_id):
                    continue
            except Exception:
                pass
            try:
                if cov and float(cov) < float(min_cov):
                    continue
            except Exception:
                pass

            is_mut = (subtype.upper() == "POINT") or (_MUT_RX.match(elem_symbol) is not None)

            phenotype = ""
            if amr_class and amr_subclass and amr_subclass.lower() != amr_class.lower():
                phenotype = f"{amr_class} / {amr_subclass}"
            elif amr_class:
                phenotype = amr_class

            rows.append({
                "sample": sample,
                "determinant": elem_symbol,                 # keep exact
                "type": "mutation" if is_mut else (rtype or "AMR"),
                "phenotype": phenotype,
                "class": amr_class,
                "subtype": subtype,
                "identity": ident,
                "coverage": cov,
                "length": _get_field(rec, headers, hmap, "Alignment length"),
                "contig": _get_field(rec, headers, hmap, "Contig id"),
                "start": _get_field(rec, headers, hmap, "Start"),
                "end": _get_field(rec, headers, hmap, "Stop"),
                "strand": _get_field(rec, headers, hmap, "Strand"),
                "source": "AMRFinderPlus",
                "database": "",
                "accession": _get_field(rec, headers, hmap, "Closest reference accession"),
                "tool": "amrfinder",
            })
    return rows

# --- ABRicate parser (resfinder / card) --------------------------------------

_ABR_HDR_MIN = {"SEQUENCE","START","END","STRAND","GENE","%COVERAGE","%IDENTITY","DATABASE","ACCESSION","PRODUCT","RESISTANCE"}

def parse_abricate_tsv(tsv_path: str, *, sample: str, database: str, min_id: float, min_cov: float) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not os.path.exists(tsv_path):
        return rows
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or not set(reader.fieldnames).intersection(_ABR_HDR_MIN):
            return rows
        for rec in reader:
            gene = (rec.get("GENE") or "").strip()
            if not gene:
                continue
            ident = (rec.get("%IDENTITY") or rec.get("IDENTITY") or "").strip()
            cov   = (rec.get("%COVERAGE") or rec.get("COVERAGE") or "").strip()
            try:
                if ident and float(ident) < float(min_id):
                    continue
            except Exception:
                pass
            try:
                if cov and float(cov) < float(min_cov):
                    continue
            except Exception:
                pass
            rows.append({
                "sample": sample,
                "determinant": gene,
                "type": "acquired",
                "phenotype": (rec.get("RESISTANCE") or rec.get("PRODUCT") or "").strip(),
                "class": "",
                "subtype": "",
                "identity": ident,
                "coverage": cov,
                "length": "",
                "contig": (rec.get("SEQUENCE") or "").strip(),
                "start": (rec.get("START") or "").strip(),
                "end": (rec.get("END") or "").strip(),
                "strand": (rec.get("STRAND") or "").strip(),
                "source": "ABRicate",
                "database": database,
                "accession": (rec.get("ACCESSION") or "").strip(),
                "tool": "abricate",
            })
    return rows

# --- Runner (matches your _run_one call) -------------------------------------

def run_amr_for_sample(
    *,
    sample: str,
    fasta: str,
    out_root: str,
    min_id: float,
    min_cov: float,
    enable_amrfinder: bool,
    enable_resfinder: bool,
    enable_card: bool,
    threads: int,
    organism: str,
    enable_ncbi: bool = False,
    logger_fn=None,
) -> List[Dict[str, str]]:
    """
    Executes AMRFinderPlus and optional ABRicate (resfinder, card),
    writes raw outputs under <out_root>/raw/..., and returns normalized AMR rows.
    Uses conda envs by name from env vars:
      VIRAMR_ENV   (for amrfinder)   default 'viramr'
      ABRICATE_ENV (for abricate)    default 'abricate'
    """
    os.makedirs(out_root, exist_ok=True)
    raw_dir = os.path.join(out_root, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    # tmp space for tools that require a writable TMPDIR
    tmp_dir = os.path.join(out_root, ".tmp", "amr", sample)
    os.makedirs(tmp_dir, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = env["TEMP"] = env["TMP"] = tmp_dir

    viramr_env = os.getenv("VIRAMR_ENV", "viramr")
    abricate_env = os.getenv("ABRICATE_ENV", "genepid")

    rows: List[Dict[str, str]] = []

    if enable_ncbi:
        out_tsv = os.path.join(raw_dir, "abricate_ncbi", f"{sample}.tsv")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        pending = out_tsv + ".partial"
        prefix = _env_prefix("abricate", viramr_env=viramr_env, abricate_env=abricate_env)
        cmd = prefix + ["--db", "ncbi", "--threads", str(threads),
                        "--minid", str(min_id), "--mincov", str(min_cov), fasta]
        if logger_fn:
            logger_fn(f"[{sample}] AMR fallback: ABRicate/ncbi (no organism argument)")
            logger_fn(f"[CMD] {' '.join(cmd)}")
        with open(pending, "w") as output:
            cp = subprocess.run(cmd, env=env, stdout=output, stderr=subprocess.PIPE, text=True)
        if cp.returncode != 0:
            if logger_fn:
                logger_fn(cp.stderr, level="Warn")
            raise subprocess.CalledProcessError(cp.returncode, cmd, stderr=cp.stderr)
        with open(pending) as output:
            if "GENE" not in output.readline().split():
                raise RuntimeError(f"[{sample}] ABRicate/ncbi did not produce a valid results table")
        os.replace(pending, out_tsv)
        rows.extend(parse_abricate_tsv(out_tsv, sample=sample, database="ncbi", min_id=min_id, min_cov=min_cov))

    # --- AMRFinderPlus ---
    if enable_amrfinder:
        out_tsv = os.path.join(raw_dir, "amrfinder", f"{sample}.tsv")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        if not os.path.exists(out_tsv):
            # use conda env for amrfinder
            prefix = _env_prefix("amrfinder", viramr_env=viramr_env, abricate_env=abricate_env)
            cmd = prefix + ["-n", fasta, "-O", (organism or "bacteria"), "--threads", str(threads), "--output", out_tsv]
            _run(cmd, env=env, logger_fn=logger_fn)
        rows.extend(parse_amrfinder_tsv(out_tsv, sample=sample, min_id=min_id, min_cov=min_cov))

    # --- ABRicate: ResFinder ---
    if enable_resfinder:
        out_tsv = os.path.join(raw_dir, "abricate_resfinder", f"{sample}.tsv")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        if not os.path.exists(out_tsv):
            prefix = _env_prefix("abricate", viramr_env=viramr_env, abricate_env=abricate_env)
            cmd = prefix + ["--db", "resfinder", fasta]
            with open(out_tsv, "w") as fout:
                cp = subprocess.run(cmd, env=env, stdout=fout, stderr=subprocess.PIPE, text=True)
                if cp.returncode != 0:
                    if logger_fn:
                        logger_fn(cp.stderr, level="Warn")
                    raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
        rows.extend(parse_abricate_tsv(out_tsv, sample=sample, database="resfinder", min_id=min_id, min_cov=min_cov))

    # --- ABRicate: CARD ---
    if enable_card:
        out_tsv = os.path.join(raw_dir, "abricate_card", f"{sample}.tsv")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        if not os.path.exists(out_tsv):
            prefix = _env_prefix("abricate", viramr_env=viramr_env, abricate_env=abricate_env)
            cmd = prefix + ["--db", "card", fasta]
            with open(out_tsv, "w") as fout:
                cp = subprocess.run(cmd, env=env, stdout=fout, stderr=subprocess.PIPE, text=True)
                if cp.returncode != 0:
                    if logger_fn:
                        logger_fn(cp.stderr, level="Warn")
                    raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
        rows.extend(parse_abricate_tsv(out_tsv, sample=sample, database="card", min_id=min_id, min_cov=min_cov))

    if logger_fn:
        logger_fn(f"[{sample}] AMR rows parsed: {len(rows)}")
    return rows
