"""GTDB bac120 placement and compact species-level tree rendering."""

from __future__ import annotations

import copy
import csv
import gzip
import json
import os
import re
import subprocess
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Sequence


SAFE_GENOME_ID = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class PlacementResult:
    status: str
    reason: str
    tree_newick: str = ""
    tree_svg: str = ""
    closest_species: tuple[str, ...] = ()
    marker_evidence: dict = field(default_factory=dict)
    outgroup: dict = field(default_factory=dict)
    interpretation: str = (
        "Marker-based placement is supporting evidence, not a confirmed species identification. "
        "Distances are evolutionary branch lengths, not ANI or percent match. "
        "This is a pruned reference-tree placement, not a newly inferred species tree; "
        "the display does not provide bootstrap support for this sample."
    )


def run_bac120_placement(
    *,
    assembly: str | Path,
    sample_id: str,
    output_directory: str | Path,
    gtdb_metadata_tsv: str | Path,
    cpus: int = 1,
    pplacer_cpus: int = 1,
    maximum_species: int = 10,
    command_prefix: Sequence[str] = ("gtdbtk",),
    use_pplacer_scratch: bool = True,
    use_full_tree: bool = False,
) -> PlacementResult:
    """Place one non-mixed assembly and render its nearest named species."""

    assembly_path = Path(assembly).resolve()
    if not assembly_path.is_file() or assembly_path.stat().st_size == 0:
        raise ValueError(f"Assembly is missing or empty: {assembly_path}")
    metadata_path = Path(gtdb_metadata_tsv).resolve()
    if not metadata_path.is_file():
        raise ValueError(f"GTDB metadata table is missing: {metadata_path}")
    output = Path(output_directory).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"GTDB placement output must be empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    genome_id = SAFE_GENOME_ID.sub("_", sample_id).strip("._-") or "query"
    batchfile = output / "batchfile.tsv"
    batchfile.write_text(f"{assembly_path}\t{genome_id}\n")
    command = [
        *command_prefix,
        "classify_wf",
        "--batchfile",
        str(batchfile),
        "--out_dir",
        str(output / "gtdbtk"),
        "--prefix",
        "placement",
        "--cpus",
        str(max(1, int(cpus))),
        "--pplacer_cpus",
        str(max(1, int(pplacer_cpus))),
        "--place_species",
        "--keep_intermediates",
    ]
    if use_pplacer_scratch:
        scratch = output / "scratch"
        scratch.mkdir()
        command.extend(["--scratch_dir", str(scratch)])
    if use_full_tree:
        command.append("--full_tree")
    try:
        # Stream logs to disk so a long placement remains diagnosable while running.
        with (output / "gtdbtk.stdout.log").open("w") as stdout, \
                (output / "gtdbtk.stderr.log").open("w") as stderr:
            completed = subprocess.run(command, text=True, stdout=stdout, stderr=stderr)
    except OSError as exc:
        result = PlacementResult(
            status="unavailable",
            reason=f"GTDB-Tk could not be started: {exc}",
        )
        _write_result(output, result)
        return result
    if completed.returncode != 0:
        result = PlacementResult(
            status="unavailable",
            reason=(
                f"GTDB-Tk placement failed with exit {completed.returncode}; "
                "see gtdbtk.stderr.log."
            ),
        )
        _write_result(output, result)
        return result

    return render_bac120_placement(
        gtdb_output=output / "gtdbtk", query_name=genome_id,
        metadata_path=metadata_path, output_directory=output,
        maximum_species=maximum_species,
    )


def render_bac120_placement(
    *, gtdb_output: str | Path, query_name: str, metadata_path: str | Path,
    output_directory: str | Path, maximum_species: int = 10, preferred_accessions=(),
) -> PlacementResult:
    """Render one sample from a completed single- or multi-sample GTDB workflow.

    This cheap publication step can be retried without rerunning marker extraction
    or placement. The source directory and its full reference trees are read-only.
    """
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    root = Path(gtdb_output)
    evidence = _marker_evidence(root, query_name)
    source_tree = _find_query_tree(root, query_name)
    if source_tree is None:
        result = PlacementResult(
            status="unavailable",
            reason="GTDB-Tk did not produce a bacterial placement tree containing the sample.",
            marker_evidence=evidence,
        )
        _write_result(output, result)
        return result
    try:
        result = _compact_tree(
            source_tree=source_tree,
            query_name=query_name,
            metadata_path=Path(metadata_path),
            output_directory=output,
            maximum_species=maximum_species,
            preferred_accessions=preferred_accessions,
        )
    except Exception as exc:
        result = PlacementResult(status="unavailable", reason=str(exc))
    evidence["tree_scope"] = "backbone" if "backbone" in source_tree.name else "class_or_full_reference"
    if evidence["tree_scope"] == "backbone" and result.status == "available":
        result = replace(result, reason="Only coarse backbone placement was available; species-neighbor resolution is limited.")
    result = replace(result, marker_evidence=evidence)
    _write_result(output, result)
    return result


def _find_query_tree(root: Path, query_name: str) -> Path | None:
    from Bio import Phylo

    candidates = sorted(root.rglob("*.bac120.classify.tree*.tree"))
    candidates.extend(sorted(root.rglob("*.bac120.classify.tree")))
    candidates.extend(sorted(root.rglob("*backbone.bac120.classify.tree")))
    # Class-level trees take precedence; backbone placement alone is coarse.
    candidates.sort(key=lambda path: "backbone" in path.name)
    for path in dict.fromkeys(candidates):
        tree = Phylo.read(str(path), "newick")
        if any(tip.name == query_name for tip in tree.get_terminals()):
            return path
    return None


def _marker_evidence(root: Path, query_name: str) -> dict:
    evidence = {"marker_set": "bac120", "total_markers": 120}
    for path in sorted(root.rglob("*.bac120.markers_summary.tsv")):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("name") != query_name:
                    continue
                for source, destination in (
                    ("number_unique_genes", "single_copy_markers"),
                    ("number_multiple_genes", "multiple_copy_markers"),
                    ("number_multiple_unique_genes", "multiple_unique_markers"),
                    ("number_missing_genes", "missing_markers"),
                ):
                    try:
                        evidence[destination] = int(row[source])
                    except (KeyError, TypeError, ValueError):
                        pass  # Unavailable measurements must not become zero.
    for path in sorted(root.rglob("*.bac120.summary.tsv")):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("user_genome") == query_name:
                    for key in ("msa_percent", "classification", "classification_method", "warnings"):
                        if row.get(key) not in (None, "", "N/A"):
                            evidence[key] = row[key]
    return evidence


def _distances_from_query(tree, query):
    """All patristic distances in O(nodes), without a root search per leaf."""
    neighbors = {}
    for parent in tree.find_clades(order="level"):
        neighbors.setdefault(parent, [])
        for child in parent.clades:
            length = child.branch_length or 0.0
            neighbors[parent].append((child, length))
            neighbors.setdefault(child, []).append((parent, length))
    distances = {query: 0.0}
    stack = [query]
    while stack:
        node = stack.pop()
        for other, length in neighbors[node]:
            if other not in distances:
                distances[other] = distances[node] + length
                stack.append(other)
    return distances


def _induced_tree(tree, keep):
    """Copy only retained branches, preserving pairwise evolutionary distances."""
    from Bio.Phylo.BaseTree import Tree

    retained = {}
    for node in tree.find_clades(order="postorder"):
        children = [retained[child] for child in node.clades if child in retained]
        if not node.clades and node.name in keep:
            retained[node] = copy.copy(node)
        elif children:
            clone = copy.copy(node)
            clone.clades = children
            if len(children) == 1:
                clone = children[0]
                clone.branch_length = (clone.branch_length or 0) + (node.branch_length or 0)
            retained[node] = clone
    root = retained[tree.root]
    root.branch_length = 0.0
    return Tree(root=root, rooted=tree.rooted)


def _compact_tree(
    *,
    source_tree: Path,
    query_name: str,
    metadata_path: Path,
    output_directory: Path,
    maximum_species: int,
    preferred_accessions=(),
) -> PlacementResult:
    from Bio import Phylo

    tree = Phylo.read(str(source_tree), "newick")
    query = next((terminal for terminal in tree.get_terminals() if terminal.name == query_name), None)
    if query is None:
        raise ValueError("The GTDB placement tree does not contain the sample leaf")
    terminals = [terminal for terminal in tree.get_terminals() if terminal is not query]
    distances = _distances_from_query(tree, query)
    terminals.sort(key=lambda terminal: (distances[terminal], terminal.name or ""))
    accessions = {_normalise_accession(terminal.name or "") for terminal in terminals}
    metadata = _load_species_names(metadata_path, accessions)
    selected = []
    seen_species: set[str] = set()
    for terminal in terminals:
        accession = _normalise_accession(terminal.name or "")
        species = metadata.get(accession)
        if not species or species in seen_species:
            continue
        selected.append((terminal.name, accession, species, distances[terminal]))
        seen_species.add(species)
        if len(selected) >= maximum_species:
            break
    if len(selected) < 2:
        raise ValueError(
            "Too few named neighboring species were available for a useful phylogenetic tree"
        )

    # Retain exact ANI nominees when present in the source placement, not a
    # different assembly sharing their species name. Always keep the nearest
    # marker neighbor, then fill remaining slots with marker neighbors.
    lookup = {_normalise_accession(t.name or ''):t for t in terminals}
    combined = selected[:1]
    for acc in tuple(preferred_accessions)[:5]:
        terminal = lookup.get(acc)
        if terminal is not None and metadata.get(acc) and terminal.name not in {r[0] for r in combined}:
            combined.append((terminal.name, acc, metadata[acc], distances[terminal]))
    for row in selected:
        if row[0] not in {r[0] for r in combined}:
            combined.append(row)
    selected = combined[:maximum_species]

    # Select a reference outside the MRCA of the displayed ingroup in the
    # *full* rooted placement tree. Never relabel the farthest ingroup tip.
    ingroup = selected[:max(2, maximum_species - 1)]
    ancestor = tree.common_ancestor(query_name, *[row[0] for row in ingroup])
    inside = set(ancestor.get_terminals())
    outgroup = {}
    for terminal in terminals:
        accession = _normalise_accession(terminal.name or "")
        species = metadata.get(accession)
        if terminal not in inside and species and species not in {row[2] for row in ingroup}:
            selected = ingroup + [(terminal.name, accession, species, distances[terminal])]
            outgroup = {"reference_accession": accession, "species": species,
                        "method": "Nearest named reference outside the displayed ingroup MRCA in the full GTDB placement tree"}
            break
    keep = {query_name, *[name for name, _accession, _species, _distance in selected]}
    compact = _induced_tree(tree, keep)
    labels = {name: f"{species} [{accession}]" for name, accession, species, _ in selected}
    labels[query_name] = query_name
    for terminal in compact.get_terminals():
        terminal.name = labels.get(terminal.name, terminal.name)
    if outgroup:
        label = f"{outgroup['species']} [{outgroup['reference_accession']}]"
        compact.root_with_outgroup(label)

    newick = output_directory / "closest_species_tree.nwk"
    Phylo.write(compact, str(newick), "newick")
    svg = output_directory / "closest_species_tree.svg"
    _render_svg(compact, svg, query_name=query_name)
    neighbors = output_directory / "closest_species_tree_neighbors.tsv"
    with neighbors.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "species", "reference_accession", "tree_distance"])
        for rank, (_name, accession, species, distance) in enumerate(selected, start=1):
            writer.writerow([rank, species, accession, f"{distance:.8f}"])
    return PlacementResult(
        status="available",
        reason="The sample was placed using the GTDB bac120 marker set.",
        tree_newick=newick.name,
        tree_svg=svg.name,
        closest_species=tuple(species for _name, _accession, species, _distance in selected),
        outgroup=outgroup,
    )


def _load_species_names(path: Path, wanted: set[str]) -> dict[str, str]:
    names: dict[str, str] = {}
    handle_context = (
        gzip.open(path, mode="rt", newline="")
        if path.suffix.casefold() == ".gz"
        else path.open(mode="r", newline="")
    )
    with handle_context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession", "gtdb_taxonomy"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("GTDB metadata is missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            accession = _normalise_accession(row.get("accession") or "")
            if accession not in wanted:
                continue
            ncbi_species = _clean_species_label(row.get("ncbi_species") or "")
            if not ncbi_species:
                ncbi_species = _taxonomy_species(row.get("ncbi_taxonomy") or "")
            gtdb_species = next(
                (
                    _clean_species_label(value)
                    for value in (row.get("gtdb_taxonomy") or "").split(";")
                    if value.startswith("s__") and len(value) > 3
                ),
                "",
            )
            species = ncbi_species or gtdb_species
            # Preserve familiar labels without merging distinct GTDB clusters
            # that happen to share one NCBI species name (e.g. cecorum/A/B/C).
            if ncbi_species and gtdb_species and ncbi_species != gtdb_species:
                species = f"{ncbi_species} (GTDB: {gtdb_species})"
            if species and not species.casefold().endswith(" sp."):
                names[accession] = species
    return names


def _clean_species_label(value: str) -> str:
    value = " ".join(str(value or "").split())
    if value.startswith("s__"):
        value = value[3:]
    # GTDB alphabetic suffixes are taxonomic identifiers, not whitespace.
    return value.strip()


def _taxonomy_species(value: str) -> str:
    return next(
        (
            _clean_species_label(rank)
            for rank in str(value or "").split(";")
            if rank.strip().startswith("s__") and len(rank.strip()) > 3
        ),
        "",
    )


def _normalise_accession(value: str) -> str:
    value = value.strip().strip("'")
    if value.startswith(("RS_", "GB_")):
        value = value[3:]
    return value


def _render_svg(tree, destination: Path, *, query_name: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    query = next(
        (terminal for terminal in tree.get_terminals() if terminal.name == query_name),
        None,
    )
    if query is not None:
        query.color = "#c62828"
        query.width = 3
    terminal_count = len(tree.get_terminals())
    figure = plt.figure(figsize=(10, max(3.5, 0.45 * terminal_count)))
    axes = figure.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False, show_confidence=False)
    axes.set_title("Closest species: bac120 phylogenetic placement")
    axes.set_xlabel("Evolutionary branch length (not percent identity)")
    axes.set_ylabel("")
    figure.tight_layout()
    figure.savefig(destination, format="svg", bbox_inches="tight")
    plt.close(figure)


def _write_result(output: Path, result: PlacementResult) -> None:
    temporary = output / "phylogenetic_placement.json.tmp"
    temporary.write_text(json.dumps(result.__dict__, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, output / "phylogenetic_placement.json")
