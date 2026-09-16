#!/usr/bin/env python3
"""Validate an immutable BactiPipe species-identification database."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.species_database import validate_database
from bactipipe.species_workflow import identification_database_kind
from bactipipe.speciesfinder_database import (
    DATABASE_KIND as SPECIESFINDER_DATABASE_KIND,
    validate_speciesfinder_database,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-root", required=True)
    parser.add_argument(
        "--allow-candidate",
        action="store_true",
        help="Permit validation of an inactive candidate database",
    )
    args = parser.parse_args(argv)
    try:
        if identification_database_kind(args.database_root) == SPECIESFINDER_DATABASE_KIND:
            database = validate_speciesfinder_database(
                args.database_root, allow_candidate=args.allow_candidate
            )
            identity = database.identity
        else:
            identity = validate_database(
                args.database_root, allow_candidate=args.allow_candidate
            )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe species database validation failed: {exc}\n")
        return 1
    if identity.get("database_kind") == SPECIESFINDER_DATABASE_KIND:
        detail = (
            f"templates={identity.get('template_count', 0)}\t"
            f"unique_accessions={identity.get('unique_template_accession_count', 0)}"
        )
    else:
        detail = (
            f"references={identity.get('reference_count', 0)}\t"
            f"species={identity.get('species_count', 0)}"
        )
    print(
        f"version={identity.get('version', 'unknown')}\t{detail}\t"
        f"status={identity.get('status', 'unknown')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
