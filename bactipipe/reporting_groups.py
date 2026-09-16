"""Controlled species-complex labels used by identification and mixture review."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from importlib.resources import files as resource_files
from pathlib import Path


RESOURCE_NAME = "species_reporting_groups.tsv"


@dataclass(frozen=True)
class ReportingGroupMember:
    group_id: str
    display_name: str
    species_taxid: str
    species_name: str
    rationale: str = ""


class ReportingGroupCatalog:
    """Validate and expose explicit groups that cannot support species resolution."""

    def __init__(self, members: list[ReportingGroupMember], checksum: str):
        self.members = tuple(members)
        self.checksum = checksum
        self._by_taxid: dict[str, ReportingGroupMember] = {}
        self._display_by_id: dict[str, str] = {}
        for member in self.members:
            if not member.group_id or not member.display_name:
                raise ValueError("Reporting-group IDs and display names cannot be blank")
            if not member.species_taxid.isdigit() or not member.species_name:
                raise ValueError("Reporting-group members require a numeric taxid and name")
            previous = self._by_taxid.get(member.species_taxid)
            if previous is not None:
                raise ValueError(
                    f"Species taxid {member.species_taxid} belongs to more than one reporting group"
                )
            previous_display = self._display_by_id.get(member.group_id)
            if previous_display is not None and previous_display != member.display_name:
                raise ValueError(
                    f"Reporting group {member.group_id!r} has inconsistent display names"
                )
            self._by_taxid[member.species_taxid] = member
            self._display_by_id[member.group_id] = member.display_name

    @classmethod
    def from_path(cls, path: str | Path) -> "ReportingGroupCatalog":
        raw = Path(path).read_bytes()
        return cls(cls._parse(raw), sha256(raw).hexdigest())

    @classmethod
    def from_package(cls) -> "ReportingGroupCatalog":
        resource = resource_files("bactipipe.data").joinpath(RESOURCE_NAME)
        with resource.open("rb") as handle:
            raw = handle.read()
        return cls(cls._parse(raw), sha256(raw).hexdigest())

    @staticmethod
    def _parse(raw: bytes) -> list[ReportingGroupMember]:
        text = raw.decode("utf-8")
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        required = {"group_id", "display_name", "species_taxid", "species_name"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Reporting-group table is missing columns: " + ", ".join(sorted(missing))
            )
        return [
            ReportingGroupMember(
                group_id=str(row.get("group_id") or "").strip(),
                display_name=" ".join(str(row.get("display_name") or "").split()),
                species_taxid=str(row.get("species_taxid") or "").strip(),
                species_name=" ".join(str(row.get("species_name") or "").split()),
                rationale=" ".join(str(row.get("rationale") or "").split()),
            )
            for row in reader
            if any(str(value or "").strip() for value in row.values())
        ]

    def display_name_for_taxid(self, species_taxid: object) -> str:
        member = self._by_taxid.get(str(species_taxid))
        return member.display_name if member else ""

    @property
    def taxid_to_display_name(self) -> dict[str, str]:
        return {
            taxid: member.display_name for taxid, member in self._by_taxid.items()
        }

    @property
    def member_names_by_display_name(self) -> dict[str, tuple[str, ...]]:
        grouped: dict[str, list[str]] = {}
        for member in self.members:
            grouped.setdefault(member.display_name, []).append(member.species_name)
        return {
            display: tuple(sorted(names, key=str.casefold))
            for display, names in grouped.items()
        }
