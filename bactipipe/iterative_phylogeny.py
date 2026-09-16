"""Checkpointed, assembly-only phylogenetic follow-up, independent of QC.

The reference installation must be immutable. Cache keys bind its release/package
identity, tool version, processing policy and the actual masked query alignment.
Only tree artifacts are reused across changed assemblies, never GTDB ANI calls.
"""
from __future__ import annotations

import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile

from .phylogenetic_placement import _marker_evidence, render_bac120_placement

POLICY = "bac120-followup-v1"
FILES = ("closest_species_tree.svg", "closest_species_tree.nwk", "closest_species_tree_neighbors.tsv")


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def key(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def marker_sequence(root):
    """Return the exact masked bacterial query MSA, not raw protein/gene hashes."""
    from Bio import SeqIO
    for path in sorted(Path(root).rglob("*.bac120.user_msa.fasta*")):
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == "query":
                    return str(record.seq).upper()
    return ""


def reference_identity(root, metadata, version):
    root = Path(root).resolve()
    # Hash the publisher's package identity and the files defining displayed labels
    # and masks. Never repeatedly hash hundreds of GB of immutable reference data.
    package = root.parent / "MD5SUM.txt"
    if not package.is_file():
        raise ValueError("GTDB installation requires its verified package MD5SUM.txt")
    return {
        "root": str(root), "release": (root / "metadata/metadata.txt").read_text(),
        "package_manifest_sha256": sha256(package), "names_sha256": sha256(metadata),
        "masks": {p.name: sha256(p) for p in sorted((root / "masks").glob("*")) if p.is_file()},
        "tool": version.strip(), "policy": POLICY,
        "implementation": {name: sha256(Path(__file__).with_name(name)) for name in ("iterative_phylogeny.py", "phylogenetic_placement.py")},
        "parameters": {"min_perc_aa": 10, "split_tree": True, "maximum_species": 10, "classification_input": "marker_only_proteins"},
    }


def tool_identity(command_prefix, env):
    """Fingerprint actual executables inside the same environment as GTDB-Tk."""
    script = '''import hashlib,json,shutil,importlib.metadata
result={"packages":{}}
for name in ("gtdbtk","dendropy","numpy","pandas","polars","scipy","biopython"):
    try: result["packages"][name]=importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError: result["packages"][name]=None
for name in ("prodigal","hmmalign","hmmsearch","pplacer","guppy","skani"):
    path=shutil.which(name)
    if not path: raise RuntimeError("Missing placement tool: "+name)
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1048576),b""): h.update(block)
    result[name]=h.hexdigest()
print(json.dumps(result,sort_keys=True))'''
    python = str(Path(command_prefix[-1]).with_name("python"))
    output = subprocess.check_output([*command_prefix[:-1], python, "-c", script], env=env, text=True, timeout=90)
    return json.loads(output)


def valid_receipt(directory):
    try:
        payload = json.loads((directory / "receipt.json").read_text())
        if not payload.get("checksums"):
            return None
        required = {*FILES, "phylogenetic_placement.json"} if payload.get("status") == "available" else {"phylogenetic_placement.json"}
        if not required.issubset(payload["checksums"]):
            return None
        for name, digest in payload["checksums"].items():
            if name not in (*FILES, "phylogenetic_placement.json") or sha256(directory / name) != digest:
                return None
        return payload
    except (OSError, ValueError, KeyError):
        return None


def execute(command, log, env, heartbeat, lock_fd, inherited_lock_fd=None):
    """Bound external lifetime; clean the process group on cancellation/soft limit."""
    with Path(log).open("a") as stream:
        proc = subprocess.Popen(
            ["timeout", "--kill-after=30s", "5h", *command], env=env,
            stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
            pass_fds=(lock_fd,) if inherited_lock_fd is None else (lock_fd, inherited_lock_fd),
        )
        try:
            while True:
                heartbeat()
                try:
                    code = proc.wait(timeout=15)
                    break
                except subprocess.TimeoutExpired:
                    continue
            if code:
                raise RuntimeError(f"GTDB stage failed (exit {code}); see placement.log")
        except BaseException:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
                proc.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            except ProcessLookupError:
                pass
            raise


def place_iteratively(*, assembly, cache_root, reference_root, metadata,
                      command_prefix=("conda", "run", "--no-capture-output", "-n", "gtdbtk", "gtdbtk"),
                      threads=8, pplacer_threads=4, heartbeat=lambda: None, inherited_lock_fd=None,
                      preferred_accessions=()):
    """Return a verified compact output directory and reuse provenance.

    Caller supplies a private per-sample cache and an immutable assembly snapshot.
    Caller alone decides whether the result is current enough to publish.
    """
    assembly = Path(assembly).resolve(strict=True)
    if any(c in str(assembly) for c in "\t\r\n"):
        raise ValueError("Assembly path is not batch-file safe")
    cache = Path(cache_root).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "GTDBTK_DATA_PATH": str(Path(reference_root).resolve())}
    version = subprocess.check_output([*command_prefix, "--version"], env=env, text=True, timeout=60)
    identity = reference_identity(reference_root, metadata, version)
    preferred_accessions = tuple(dict.fromkeys(str(acc) for acc in preferred_accessions))[:5]
    identity['parameters']['preferred_accessions'] = list(preferred_accessions)
    identity["executables"] = tool_identity(command_prefix, env)
    identity["parameters"].update(threads=max(1, min(int(threads), 16)),
                                   pplacer_threads=max(1, min(int(pplacer_threads), 4)))
    source_sha = sha256(assembly)
    assembly_key = key([identity, source_sha])
    target = cache / "assemblies" / assembly_key
    with (cache / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        heartbeat()
        # The inherited lock stays held by a live external stage after a parent
        # crash. Once acquired here, orphaned work directories are no longer active.
        for orphan in cache.glob("work-*"):
            if orphan.is_dir() and not orphan.is_symlink() and (orphan / ".bactipipe-owned").is_file():
                shutil.rmtree(orphan)
        receipt = valid_receipt(target)
        if receipt:
            return {**receipt, "directory": str(target), "reuse": "assembly"}
        with tempfile.TemporaryDirectory(prefix="work-", dir=cache) as temporary:
            work = Path(temporary)
            (work / ".bactipipe-owned").write_text(POLICY)
            batch = work / "batch.tsv"
            batch.write_text(f"{assembly}\tquery\n")
            gtdb = work / "gtdbtk"
            log = cache / "placement.log"
            cpus = str(max(1, min(int(threads), 16)))
            def run(arguments):
                execute([*command_prefix, *arguments], log, env, heartbeat, lock.fileno(), inherited_lock_fd)
            run(["identify", "--batchfile", str(batch), "--out_dir", str(gtdb),
                 "--prefix", "placement", "--cpus", cpus])
            run(["align", "--identify_dir", str(gtdb), "--out_dir", str(gtdb),
                 "--prefix", "placement", "--cpus", cpus, "--min_perc_aa", "10"])
            sequence = marker_sequence(gtdb)
            evidence = _marker_evidence(gtdb, "query")
            if sequence:
                evidence["msa_percent"] = round(100 * sum(c not in "-X.?" for c in sequence) / len(sequence), 2)
            marker_key = key([identity, sequence]) if sequence else ""
            marker_cache = cache / "markers" / marker_key
            marker_receipt = valid_receipt(marker_cache) if marker_key else None
            compact = work / "compact"
            compact.mkdir()
            if marker_receipt:
                for name in FILES:
                    shutil.copy2(marker_cache / name, compact / name)
                payload = json.loads((marker_cache / "phylogenetic_placement.json").read_text())
                evidence["tree_scope"] = payload.get("marker_evidence", {}).get("tree_scope", "")
                # Do not reuse whole-assembly classifications/ANI or prior marker counts.
                payload["marker_evidence"] = evidence
                reuse = "markers"
            elif sequence:
                (work / "scratch").mkdir()
                protein = gtdb / "identify/intermediate_results/marker_genes/query/query_protein.faa"
                if not protein.is_file():
                    raise ValueError("Predicted marker-workflow proteins are missing")
                protein_batch = work / "proteins.tsv"
                protein_batch.write_text(f"{protein}\tquery\n")
                # Disable whole-genome ANI only in this supporting placement stage.
                # Otherwise a changed assembly could affect a supposedly marker-only cache.
                run(["classify", "--batchfile", str(protein_batch), "--align_dir", str(gtdb),
                     "--out_dir", str(gtdb), "--prefix", "placement", "--cpus", cpus,
                     "--pplacer_cpus", str(max(1, min(int(pplacer_threads), 4))),
                     "--genes", "--place_species", "--scratch_dir", str(work / "scratch")])
                result = render_bac120_placement(gtdb_output=gtdb, query_name="query",
                    metadata_path=metadata, output_directory=compact,
                    preferred_accessions=preferred_accessions)
                payload = dict(result.__dict__)
                # GTDB classification is supporting, not the BactiPipe species verdict.
                evidence["tree_scope"] = payload.get("marker_evidence", {}).get("tree_scope", "")
                payload["marker_evidence"] = evidence
                reuse = "none"
            else:
                payload = {"status": "unavailable", "reason": "Too few usable bacterial markers for placement.",
                           "marker_evidence": evidence}
                reuse = "none"
            if sha256(assembly) != source_sha:
                raise ValueError("Assembly changed during placement")
            payload["provenance"] = {"assembly_sha256": source_sha, "marker_key": marker_key,
                                     "reference": identity, "reuse": reuse}
            atomic_json(compact / "phylogenetic_placement.json", payload)
            receipt = {"assembly_sha256": source_sha, "marker_key": marker_key,
                       "reference_key": key(identity), "reuse": reuse,
                       "status": payload["status"],
                       "checksums": {p.name: sha256(p) for p in compact.iterdir() if p.name in (*FILES, "phylogenetic_placement.json")}}
            atomic_json(compact / "receipt.json", receipt)
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                shutil.rmtree(target)  # Exact corrupt cache entry, not user input.
            shutil.copytree(compact, target)
            if marker_key and payload["status"] == "available" and not marker_receipt:
                marker_cache.parent.mkdir(parents=True, exist_ok=True)
                if marker_cache.exists():
                    shutil.rmtree(marker_cache)
                shutil.copytree(compact, marker_cache)
            return {**receipt, "directory": str(target)}
