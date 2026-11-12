#!/usr/bin/env python3
import os, sys, argparse, json, subprocess
import multiprocessing as mp
from datetime import datetime
from pathlib import Path

from bactipipe.reproducibility import database_identity, source_identity
from bactipipe.scripts.traits_organisms import (
    VF_DATABASES, amrfinder_organism, display_organism_name, parse_amrfinder_organisms, virulence_databases,
)

# Reuse your existing BactiPipe utilities
from bactipipe.scripts.utils import time_print, logger as file_logger

from bactipipe.scripts import (
    traits_amr,
    traits_vf,
    traits_reconcile,
    traits_report,
    traits_db,
)

VIRAMR_ENV = os.getenv("VIRAMR_ENV", "viramr")
ABRICATE_ENV = os.getenv("ABRICATE_ENV", "genepid")

SUPPORTED_VF_ORGS = set(VF_DATABASES)

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs['width'] = 100
        super().__init__(*args, **kwargs)
    def _get_help_string(self, action):
        if action.default is not None and action.default != argparse.SUPPRESS:
            return super()._get_help_string(action)
        return action.help

def _conda_run(env_name: str, command: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["conda", "run", "-n", env_name, *command],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

def _list_amrfinder_organisms() -> list[str]:
    """
    Return verified AMRFinderPlus organism tokens, parsed only from the
    'Available --organism options:' line, ignoring Software/Database/Running/timing lines.
    """
    cp = _conda_run(VIRAMR_ENV, ["amrfinder", "--list_organisms"])
    if cp.returncode != 0:
        raise RuntimeError("AMRFinderPlus organism query failed: " + (cp.stderr or cp.stdout or "unknown error").strip())
    return parse_amrfinder_organisms(cp.stdout or "", cp.stderr or "")

def _format_organism_list(orgs: list[str]) -> str:
    """One copyable organism per line; tool selection is an internal concern."""
    names = {display_organism_name(o) for o in set(orgs) | SUPPORTED_VF_ORGS}
    return "\n".join(sorted(names, key=str.casefold) + ["Other"])

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="bactipipe detect",
        description="Detect AMR determinants and virulence genes from assembled genomes.",
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    req = p.add_argument_group("Required arguments")
    opt = p.add_argument_group("Optional arguments")

    # Make these optional at parse-time; we'll validate unless --organism-list is used.
    req.add_argument("--analyst",
        help="Name/initials of the analyst (required).")
    req.add_argument("--genomes_dir", "--genomes-dir", dest="genomes_dir",
        help="Directory of *.fasta assemblies (required).")
    req.add_argument("--sample_sheet", "--sample-sheet", dest="sample_sheet",
        help="TSV (no header): sample_id\\tdisplay_name\\tspecimen. "
             "(required).")
    req.add_argument("--run-name",
        help="Run name identifier, e.g., 250101_10_SEQ10_XYZ (required).")

    opt.add_argument("--accession", help="Run accession identifier (e.g., P25XXXXX)", default="Unknown")
    opt.add_argument("--outdir", help="Output directory. Default: ./Detect_<run-name>/<accession>")

    # Database root shared by CGE tools; expects subdirs like virulencefinder_db, resfinder_db, etc.
    opt.add_argument("--db-dir",
                     default=os.getenv("BACTIPIPE_DB_DIR"),
                     help="Root directory that contains virulencefinder_db, resfinder_db, etc.")
    req.add_argument("--organism",
                     help='Organism name, e.g. Salmonella or "Enterococcus faecium". See --organism-list; use a genus name for an unlisted organism.')
    opt.add_argument("--organism-list",
                     action="store_true",
                     help="List organism names and exit.")

    opt.add_argument("--threads", type=int, default=8, help="CPU threads")
    opt.add_argument("--jobs", type=int, default=4, help="Maximum samples processed concurrently")
    opt.add_argument("--sample-order", choices=("display", "sheet"), default="display", help="Report order: natural display-name order (sample ID when unnamed), or input sheet order; applies to matrix columns too.")
    # Escape % in help strings for argparse
    opt.add_argument("--min_identity", "--min-identity", dest="min_identity", type=float, default=90.0, help="Minimum %% identity")
    opt.add_argument("--min_coverage", "--min-coverage", dest="min_coverage", type=float, default=60.0, help="Minimum %% coverage")

    mode = opt.add_mutually_exclusive_group()
    mode.add_argument("--amr-only", action="store_true", help="Run AMR only")
    mode.add_argument("--virulence-only", action="store_true", help="Run Virulence only")

    opt.add_argument("--extended-virulence", action="store_true",
                 help="For Escherichia: use full virulence database; "
                      "default is diagnostic, PCR-style for Pathotype-specific virulence factors.")

    # Tool selectors (AMR) — AMRFinder runs by default; ABRicate toggles add-on if flagged
    opt.add_argument("--amrfinder", action="store_true",
                     help="(Optional) Explicitly enable AMRFinderPlus (on by default unless --virulence-only).")
    opt.add_argument("--abricate-resfinder", action="store_true",
                     help="Enable ABRicate/ResFinder (in addition to AMRFinder).")
    opt.add_argument("--abricate-card", action="store_true",
                     help="Enable ABRicate/CARD (in addition to AMRFinder).")

    # Tool selectors (Virulence) — if neither is given, choose based on organism
    opt.add_argument("--virulencefinder", action="store_true",
                     help="Enable VirulenceFinder (DTU Python tool).")
    opt.add_argument("--abricate-vfdb", action="store_true",
                     help="Enable ABRicate/VFDB (Set A). If provided, will be used regardless of the organism.")
    opt.add_argument("--no-pdf", action="store_true",
                     help="Do not generate PDF reports.")

    opt.add_argument("-h","--help", action="help", help="Show this help message and exit")
    return p.parse_args(argv)

def _read_sample_sheet(sample_sheet: str, genomes_dir: str, *, sample_order='display') -> list[tuple[str, str, str, str]]:
    """
    Returns a list of (isolate_id, sample_display, specimen, fasta_path).
    - sample_sheet columns: sample_id \t display_name \t specimen  (no header)
    - fasta_path = genomes_dir/{isolate}.fasta
    Only rows where FASTA exists are kept. The sample ID is the stable key used
    for assembly and result filenames; the display name is presentation-only.
    """
    rows: list[tuple[str, str, str, str]] = []
    missing: list[str] = []
    with open(sample_sheet, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip('\r\n')
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            isolate = parts[0].strip()
            sample = parts[1].strip()
            specimen = parts[2].strip() if len(parts) > 2 else ""
            fasta = os.path.join(genomes_dir, f"{isolate}.fasta")
            if os.path.isfile(fasta):
                rows.append((isolate, sample, specimen, fasta))
            else:
                missing.append(isolate)
    if missing:
        sys.stderr.write(f"WARNING: FASTA not found for isolates (skipped): {', '.join(missing)}\n")
    if not rows:
        sys.stderr.write("ERROR: No valid rows found in --sample_sheet (no matching FASTA files).\n")
        sys.exit(2)
    if sample_order == 'display':
        from bactipipe.primary_downstream import report_order
        rows = report_order(rows)
    rows = [(sid, display or sid, specimen, fasta) for sid, display, specimen, fasta in rows]
    isolate_ids = [row[0] for row in rows]
    display_names = [row[1] for row in rows]
    if len(isolate_ids) != len(set(isolate_ids)):
        raise ValueError("Duplicate isolate IDs are not allowed in --sample-sheet.")
    if len(display_names) != len(set(display_names)):
        raise ValueError("Duplicate sample display names are not allowed in --sample-sheet.")
    return rows

def _run_one(args_pack):
    isolate_id, sample_display, specimen, fasta, cfg = args_pack
    out_root = cfg["outdir"]
    log = os.path.join(out_root, "logs", f"{isolate_id}.log")
    os.makedirs(os.path.dirname(log), exist_ok=True)

    def logfn(msg, level="Norm", mode="timestamp"):
        file_logger(log, msg, message_type=level, mode=mode)

    logfn(f"Starting traits-detect for isolate={isolate_id} (report name={sample_display})")
    os.makedirs(os.path.join(out_root, "raw"), exist_ok=True)

    amr_rows = []
    vf_rows = []

    # Internally, we name outputs by isolate_id (matches FASTA and raw/ filenames)
    if not cfg["vir_only"]:
        amr_rows = traits_amr.run_amr_for_sample(
            sample=isolate_id,
            fasta=fasta,
            out_root=out_root,
            min_id=cfg["min_id"],
            min_cov=cfg["min_cov"],
            enable_amrfinder=cfg["amrfinder"],
            enable_resfinder=cfg["resfinder"],
            enable_card=cfg["card"],
            enable_ncbi=cfg["ncbi"],
            threads=cfg["threads"],
            organism=cfg["amrfinder_organism"] or cfg["organism"],
            logger_fn=logfn,
        )

    if not cfg["amr_only"]:
        vf_rows = traits_vf.run_vf_for_sample(
            sample=isolate_id,
            fasta=fasta,
            out_root=out_root,
            min_id=cfg["min_id"],
            min_cov=cfg["min_cov"],
            enable_vfinder=cfg["virulencefinder"],
            enable_vfdb=cfg["vfdb"],
            enable_ecoli_vf=cfg["ecoli_vf"],
            threads=cfg["threads"],
            logger_fn=logfn,
            db_root=cfg["db_root"],
            organism=cfg["organism"],
            extended_vf=cfg["extended_vf"], 
        )

    merged = traits_reconcile.merge_for_sample(isolate_id, amr_rows, vf_rows)

    # Per-sample TSVs are named by isolate_id to match raw outputs
    traits_report.write_sample_tsvs(isolate_id, merged, out_root)

    # Return using the display name as the key for reporting
    return isolate_id, sample_display, specimen, merged

def main(argv=None):
    args = parse_args(argv)

    # --- short-circuit: just list organisms and exit ---
    if args.organism_list:
        try:
            orgs = _list_amrfinder_organisms()
        except (OSError, RuntimeError) as exc:
            sys.stderr.write(f"ERROR: {exc}\n")
            sys.exit(2)
        print(_format_organism_list(orgs))
        sys.exit(0)

    # --- validate required args only for real runs ---
    missing = []
    if not args.analyst:      missing.append("--analyst")
    if not args.genomes_dir:  missing.append("--genomes_dir")
    if not args.sample_sheet: missing.append("--sample_sheet")
    if not args.organism:    missing.append("--organism")
    if not args.run_name:     missing.append("--run-name")

    if missing:
        sys.stderr.write("ERROR: missing required arguments: " + ", ".join(missing) + "\n")
        sys.exit(2)

    if not args.outdir:
        args.outdir = os.path.join(f"Detect_{args.run_name}", args.accession)

    os.makedirs(args.outdir, exist_ok=True)

    # --- discover samples from sample_sheet (maps isolate -> sample display) ---
    sample_rows = _read_sample_sheet(args.sample_sheet, args.genomes_dir, sample_order=args.sample_order)
    # sample_rows: list of (isolate_id, sample_display, specimen, fasta_path)

    # --- decide default tool selection per your logic ---
    # AMR defaults: AMRFinderPlus ON (unless virulence-only); ABRicate toggles add-on if flagged
    supported = _list_amrfinder_organisms() if not args.virulence_only else []
    amr_organism = amrfinder_organism(args.organism, supported)
    amr_amrfinder = bool(amr_organism) and not args.virulence_only
    amr_ncbi = not amr_amrfinder and not args.virulence_only
    if args.amrfinder and not amr_organism:
        sys.stderr.write(f"{args.organism}: unsupported AMRFinderPlus organism; using ABRicate/ncbi.\n")

    amr_resfinder = bool(args.abricate_resfinder)
    amr_card = bool(args.abricate_card)

    # Virulence defaults:
    vf_vfinder = bool(args.virulencefinder)
    vf_vfdb = bool(args.abricate_vfdb)
    if not vf_vfinder and not vf_vfdb:
        if virulence_databases(args.organism):
            vf_vfinder = True
        else:
            vf_vfdb = True
    if vf_vfinder and not virulence_databases(args.organism):
        vf_vfinder, vf_vfdb = False, True

    # Respect mode flags
    run_amr = not args.virulence_only
    run_vf  = not args.amr_only

    # Escherichia => run VirulenceFinder + VFDB + ecoli_vf
    use_ecoli_vf = False
    if run_vf and args.organism.casefold() == "escherichia":
        vf_vfinder = True
        vf_vfdb = True
        use_ecoli_vf = True

    # If virulence disabled by mode, force off all virulence tools
    if not run_vf:
        vf_vfinder = False
        vf_vfdb = False
        use_ecoli_vf = False

    # AMR tool flags (unchanged semantics, but compute once)
    use_amrfinder = bool(amr_amrfinder) and run_amr
    use_resfinder = bool(args.abricate_resfinder) and run_amr
    use_card      = bool(args.abricate_card) and run_amr

    print("=== Settings ===")
    print(f"AMR detection: {'ON' if run_amr else 'OFF'}")
    if run_amr:
        print(f"  AMRFinderPlus: {'ON' if use_amrfinder else 'OFF'}")
        print(f"  ABRicate/NCBI: {'ON' if amr_ncbi else 'OFF'}")
        print(f"  ABRicate/ResFinder: {'ON' if use_resfinder else 'OFF'}")
        print(f"  ABRicate/CARD: {'ON' if use_card else 'OFF'}")
    print(f"Virulence detection: {'ON' if run_vf else 'OFF'}")
    if run_vf:
        print(f"  VirulenceFinder: {'ON' if vf_vfinder else 'OFF'}")
        print(f"  ABRicate/VFDB: {'ON' if vf_vfdb else 'OFF'}")
        print(f"  ABRicate/ecoli_vf: {'ON' if use_ecoli_vf else 'OFF'}")
    print("================================")

    # --- versions / db tracking ---
    versions = traits_db.ensure_and_collect_versions(
        want_amr = use_amrfinder,
        want_ncbi = amr_ncbi,
        want_vf  = run_vf,
        want_resfinder = bool(args.abricate_resfinder) and run_amr,
        want_card      = bool(args.abricate_card) and run_amr,
        want_vfinder   = vf_vfinder,
        want_vfdb      = vf_vfdb,
        want_ecoli_vf  = use_ecoli_vf,
        db_root = args.db_dir,
    )
    versions["thresholds"] = {
        "min_identity": float(args.min_identity),
        "min_coverage": float(args.min_coverage),
    }

    with open(os.path.join(args.outdir, "versions.json"), "w") as f:
        json.dump(versions, f, indent=2)

    cfg = dict(
        outdir=os.path.abspath(args.outdir),
        db_root=args.db_dir,
        min_id=float(args.min_identity),
        min_cov=float(args.min_coverage),
        threads=int(args.threads),

        amr_only=bool(args.amr_only),
        vir_only=bool(args.virulence_only),

        # Final tool toggles after default logic
        amrfinder=bool(amr_amrfinder),
        ncbi=bool(amr_ncbi),
        amrfinder_organism=amr_organism,
        resfinder=bool(amr_resfinder),
        card=bool(amr_card),

        virulencefinder=bool(vf_vfinder),
        vfdb=bool(vf_vfdb),
        ecoli_vf=bool(use_ecoli_vf),
        extended_vf=bool(args.extended_virulence),
        organism=args.organism,
        analyst=args.analyst,
        run_name=args.run_name,
        accession=args.accession,
    )

    started_at = datetime.now().astimezone()
    time_print(f"traits-detect started at {started_at.strftime('%Y-%m-%d %H:%M:%S')}", "Header")
    os.makedirs(os.path.join(cfg["outdir"], "raw"), exist_ok=True)
    os.makedirs(os.path.join(cfg["outdir"], "reports"), exist_ok=True)

    # Treat --threads as the total CPU budget, not a per-sample allocation.
    jobs = max(1, min(len(sample_rows), args.jobs, cfg["threads"]))
    cfg["threads"] = max(1, cfg["threads"] // jobs)
    pool = mp.Pool(processes=jobs)
    try:
        results = pool.map(
            _run_one,
            [(iso, samp, specimen, fa, cfg) for (iso, samp, specimen, fa) in sample_rows],
        )
    finally:
        pool.close(); pool.join()

    # --- unified run-level reporting ---
    # Key reports by the 'sample' (display name from column 2)
    all_merged = {sample_display: merged for (_iso, sample_display, _specimen, merged) in results}

    # Unified TSV for convenience
    traits_report.write_run_summary(all_merged, cfg["outdir"])

    # Single PDF for the entire run
    if not args.no_pdf:
        traits_report.render_run_pdf(
            out_root=cfg["outdir"],
            run_name=cfg["run_name"],
            accession=cfg["accession"],
            tech_name=cfg["analyst"],
            versions=versions,
            thresholds=versions.get("thresholds", {}),
            all_merged=all_merged,
            pipeline_title="Virulence and AMR Gene Detection",
            run_amr=(not args.virulence_only),
            run_vf=(not args.amr_only),
        )
    # Also emit a consolidated matrix TSV (no header)
    accession = cfg["accession"]
    traits_report.write_consolidated_tsv(all_merged, cfg["outdir"], f"{accession}_vir_amr_detect.tsv")

    completed_at = datetime.now().astimezone()
    result_payload = {
        "schema_version": 1,
        "output_contract_version": "bactipipe-detect-v1",
        "pipeline": source_identity(),
        "database": database_identity(cfg["db_root"]),
        "command": "detect",
        "status": "completed",
        "run_name": cfg["run_name"],
        "accession": cfg["accession"],
        "organism": cfg["organism"],
        "started_at": started_at.isoformat(),
        "completed_at": completed_at.isoformat(),
        "parameters": {
            "min_identity": cfg["min_id"],
            "min_coverage": cfg["min_cov"],
            "jobs": jobs,
            "threads_per_sample": cfg["threads"],
            "amr": run_amr,
            "virulence": run_vf,
            "amrfinder": cfg["amrfinder"],
            "amrfinder_organism": cfg["amrfinder_organism"],
            "abricate_ncbi": cfg["ncbi"],
            "virulencefinder_databases": virulence_databases(cfg["organism"]) if cfg["virulencefinder"] else None,
            "abricate_resfinder": cfg["resfinder"],
            "abricate_card": cfg["card"],
            "virulencefinder": cfg["virulencefinder"],
            "abricate_vfdb": cfg["vfdb"],
            "extended_virulence": cfg["extended_vf"],
        },
        "samples": [
            {
                "sample_id": isolate_id,
                "display_name": sample_display,
                "specimen": specimen,
                "amr_tsv": f"{isolate_id}.amr.tsv",
                "virulence_tsv": f"{isolate_id}.vf.tsv",
                "summary_tsv": f"{isolate_id}.summary.tsv",
            }
            for isolate_id, sample_display, specimen, _merged in results
        ],
        "outputs": {
            "run_summary_tsv": "run.summary.tsv",
            "consolidated_tsv": f"reports/{accession}_vir_amr_detect.tsv",
            "consolidated_pdf": None if args.no_pdf else f"reports/{accession}_vir_amr_detect.pdf",
            "versions": "versions.json",
        },
    }
    with open(os.path.join(cfg["outdir"], "analysis_result.json"), "w", encoding="utf-8") as handle:
        json.dump(result_payload, handle, indent=2)


    time_print("traits-detect completed.", "Header")

if __name__ == "__main__":
    main()
