#!/usr/bin/env python3

import sys
import shlex
from shutil import which
import subprocess
import os

from bactipipe.__version__ import __version__
from bactipipe.scripts.core import (
    list_databases,
    update_databases,
    update_packages,
    check_updates,
    DB_ROOT,
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
    place-assembly     :  Run cached marker placement separately from QC reporting.
    render-placement   :  Redraw retained tree and ANI figures without reanalysis.
    build-identify-db  :  Build an inactive reviewed species-reference candidate.
    plan-identify-refs :  Propose canonical references from NCBI assembly reports.
    plan-identify-inventory :  Propose veterinary and near-neighbor assembly accessions.
    fetch-identify-reports :  Fetch current NCBI reports for proposed assemblies.
    reconcile-identify-inventory :  Resolve current accessions and canonical fallbacks.
    download-identify-refs :  Download reconciled assemblies and write a reviewed manifest.
    augment-identify-inventory :  Add audited controlled-gap reference candidates.
    write-identify-manifest :  Validate downloads and write the final candidate manifest.
    validate-identify-manifest :  Verify all candidate references, roles, and checksums.
    validate-identify-db :  Verify an immutable production or candidate database.
    prepare-speciesfinder-db :  Prepare the official bacterial SpeciesFinder candidate database.

    Package / Environment Management:
    check-updates      :  Check for available updates to environments and databases.
    update-packages    :  Update all Python and Conda packages in the three environments.

    Database Management:
    list-db            :  List all databases required by BactiPipe and their installation status.
    update-db <dbs>    :  Install or update one or more databases.
                            Example: update-db virulencefinder,amrfinder,kmerfinder
                            Use "update-db all" to install/update every database.

    Other:
    -v, --version      :  Show the current version of BactiPipe.
    -h, --help         :  Show this help message.

    Use 'bactipipe <command> -h' to see help for a specific pipeline.\n'''


    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print(usage)
        sys.exit(1)

    if len(sys.argv) > 1 and sys.argv[1] in ("-v", "--version"):
        print(f"bactipipe {__version__}")
        sys.exit(0)


    cmd = sys.argv[1].lower()

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
        try:
            update_databases(sys.argv[2:])
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
    run_command = ["python3", script_to_run] + sys.argv[2:]

    try:
        subprocess.run(run_command, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
