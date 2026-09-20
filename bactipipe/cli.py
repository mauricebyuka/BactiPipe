#!/usr/bin/env python3

import sys
import subprocess
import os

from bactipipe.__version__ import __version__
from bactipipe.scripts.core import (
    list_databases,
    update_databases,
    update_packages,
    check_updates,
)

def main():
    # if config.should_check_for_updates():
    #     updater.check_outdated_pip()
    #     updater.check_outdated_conda()
    #     updater.check_database_versions()
    #     config.update_last_checked()
    
    usage = '''\nUsage: bactipipe <command> [options]

    Pipeline Commands:
    qc-illumina        :  Run the Illumina QC pipeline.
    qc-nanopore        :  Run the Nanopore QC pipeline.
    relate             :  Run genome typing and strain relatedness analysis.
    detect             :  Detect virulence and antimicrobial resistance genes in assemblies.
    identify           :  Identify an assembled bacterial genome with ANI/AF evidence.
                          --place: cached marker placement/confirmation.
                          --render: redraw retained tree and ANI figures.

    Package / Environment Management:
    check-updates      :  Check for available updates to environments and databases.
    update-packages    :  Update all Python and Conda packages in the three environments.

    Database Management:
    list-db            :  List all databases required by BactiPipe and their installation status.
    update-db <dbs>    :  Install or update one or more databases.
                            Example: update-db virulencefinder,amrfinder,kmerfinder
                            Use "update-db all" for the legacy managed databases.
                            SpeciesFinder: update-db speciesfinder --help
                            (explicit verified inputs; not included in "all").

    Other:
    -v, --version      :  Show the current version of BactiPipe.
    -h, --help         :  Show this help message.

    Use 'bactipipe <command> -h' to see help for a specific pipeline.\n'''


    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print(usage)
        sys.exit(0 if len(sys.argv) > 1 else 1)

    if len(sys.argv) > 1 and sys.argv[1] in ("-v", "--version"):
        print(f"bactipipe {__version__}")
        sys.exit(0)


    cmd = sys.argv[1].lower()
    forwarded = sys.argv[2:]
    if cmd == "identify":
        modes = [flag for flag in ("--place", "--render") if flag in forwarded]
        if len(modes) > 1:
            sys.stderr.write("Choose only one of identify --place or --render.\n")
            sys.exit(2)
        if modes:
            forwarded = [arg for arg in forwarded if arg != modes[0]]
            cmd = "place-assembly" if modes[0] == "--place" else "render-placement"
    if cmd == "update-db" and forwarded and forwarded[0].casefold() == "speciesfinder":
        forwarded = forwarded[1:]
        if "--verify" in forwarded:
            cmd = "validate-identify-db"
            forwarded = [arg for arg in forwarded if arg != "--verify"]
        else:
            cmd = "prepare-speciesfinder-db"

    if cmd == "check-updates":
        check_updates()
        sys.exit(0)

    if cmd == "update-packages":
        try:
            update_packages()
            sys.exit(0)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"\n[bactipipe:update-packages] ERROR (exit {e.returncode}): {e}\n")
            sys.exit(e.returncode)

    if cmd == "list-db":
        list_databases()
        sys.exit(0)

    if cmd == "update-db":
        if not forwarded or forwarded[0] in ("-h", "--help"):
            print("Usage: bactipipe update-db <dbs>\n"
                  "  Legacy databases: cgmlstfinder, mlst, serotypefinder, virulencefinder, amrfinder, abricate, kmerfinder\n"
                  "  SpeciesFinder preparation: bactipipe update-db speciesfinder --help\n"
                  "  SpeciesFinder validation: bactipipe update-db speciesfinder --verify --help\n"
                  "  'all' updates legacy tool-managed databases, not immutable SpeciesFinder versions.")
            sys.exit(0)
        try:
            update_databases(forwarded)
            sys.exit(0)
        except SystemExit as e:
            # _parse_db_names may raise SystemExit with a message
            if str(e):
                sys.stderr.write(str(e) + "\n")
            sys.exit(e.code or 1)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"\n[bactipipe:update-db] ERROR (exit {e.returncode}): {e}\n")
            sys.exit(e.returncode)


    tool_manage_commands = ["list-db", "update-db", "check-updates", "update-packages"]

    # Map run commands to the corresponding script
    script_dir = os.path.join(os.path.dirname(__file__), "scripts")
    script_map = {
        "qc-illumina": os.path.join(script_dir, "run_illumina_qc.py"),
        "qc-nanopore": os.path.join(script_dir, "run_nanopore_qc.py"),
        "relate": os.path.join(script_dir, "type_genomes.py"),
        "detect": os.path.join(script_dir, "run_traits.py"),
        "identify": os.path.join(script_dir, "identify_species.py"),
        "place-assembly": os.path.join(script_dir, "place_assembly.py"),
        "render-placement": os.path.join(script_dir, "render_placement.py"),
        "build-identify-db": os.path.join(script_dir, "build_species_database.py"),
        "plan-identify-refs": os.path.join(script_dir, "plan_species_references.py"),
        "plan-identify-inventory": os.path.join(script_dir, "plan_species_inventory.py"),
        "fetch-identify-reports": os.path.join(script_dir, "fetch_species_reports.py"),
        "reconcile-identify-inventory": os.path.join(script_dir, "reconcile_species_inventory.py"),
        "download-identify-refs": os.path.join(script_dir, "download_species_references.py"),
        "augment-identify-inventory": os.path.join(script_dir, "augment_species_inventory.py"),
        "write-identify-manifest": os.path.join(script_dir, "write_species_reference_manifest.py"),
        "validate-identify-manifest": os.path.join(script_dir, "validate_species_reference_manifest.py"),
        "validate-identify-db": os.path.join(script_dir, "validate_species_database.py"),
        "prepare-speciesfinder-db": os.path.join(script_dir, "prepare_speciesfinder_database.py"),
    }

    if cmd not in script_map and cmd not in tool_manage_commands:
        print(f"\nError: Invalid command '{cmd}'.")
        print(usage)
        sys.exit(1)

    elif cmd in script_map:
        # Ensure databases are present before running any pipeline
        script_to_run = script_map[cmd]
    # Forward all remaining arguments to the appropriate script
    # Long historical command names remain accepted but are not advertised as
    # normal workflow commands. Existing installations and web jobs keep working.
    run_command = [sys.executable, script_to_run] + forwarded

    try:
        subprocess.run(run_command, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
