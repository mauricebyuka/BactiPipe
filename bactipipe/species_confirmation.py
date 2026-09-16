"""Whole-assembly ANI rescue before optional, expensive marker placement.

Search the installed immutable GTDB skani sketches; no new genome database or
FASTQ screen is needed. NCBI species is the reporting identity, GTDB cluster
names remain audit evidence. Never infer ANI from a marker-only cache.
"""
from dataclasses import asdict
import csv
import fcntl
import gzip
import json
import math
import os
from pathlib import Path
import re
import subprocess

from .iterative_phylogeny import atomic_json, execute, key, sha256, place_iteratively
from .species_identification import AniHit, IdentificationPolicy, decide_identification

POLICY = "whole-assembly-confirmation-v1"


def accession(value):
    matches = set(re.findall(r"GC[AF]_\d+\.\d+", str(value)))
    if len(matches) != 1:
        raise ValueError("ANI reference does not resolve to one versioned assembly")
    return matches.pop()


def reference_labels(metadata, wanted):
    labels = {}
    opener = gzip.open if str(metadata).endswith('.gz') else open
    with opener(metadata, 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            acc = accession(row['accession'])
            if acc not in wanted:
                continue
            labels[acc] = row
            if len(labels) == len(wanted):
                break
    # The packaged sketch also contains Archaea. Unknown/unavailable taxonomy
    # must remain a competitor, never disappear and permit a bacterial call.
    return labels


def reference_name(row, acc):
    species = next((v[3:] for v in row.get('ncbi_taxonomy', '').split(';') if v.startswith('s__')), '')
    taxid = row.get('ncbi_species_taxid', '')
    valid = bool(re.fullmatch(r'[A-Z][a-z-]+ [a-z][a-z-]+', species) and taxid.isdigit())
    valid = valid and not re.search(r'\b(uncultured|bacterium|archaeon|sp)\b', species, re.I)
    # GTDB representatives remain eligible regardless of isolate/MAG origin.
    # Source category is audit evidence, not a veto of a resolved species name.
    # Preserve clinically meaningful explicitly annotated subspecies. A generic
    # S. equi reference does not confirm either controlled subspecies.
    if species == 'Streptococcus equi':
        name = row.get('ncbi_organism_name', '')
        sub = re.search(r'Streptococcus equi subsp\. (equi|zooepidemicus)\b', name)
        if sub:
            species = sub.group(0)
            taxid = {'equi': '148942', 'zooepidemicus': '40041'}[sub.group(1)]
        else:
            valid = False
    return (species if valid else f'Unresolved reference {acc}',
            taxid if valid else f'unresolved:{acc}', valid)


def evaluate_ani(rows, labels, *, evidence=None, sample='query', expected='Unknown'):
    evidence = evidence or {}
    policy = IdentificationPolicy.from_saved_thresholds(evidence.get('identification_thresholds'))
    hits, audit, callable_accessions = [], [], set()
    for row in rows:
        acc = accession(row['Ref_file'])
        metadata = labels.get(acc, {})
        name, taxid, callable_name = reference_name(metadata, acc)
        values = [float(row[k]) for k in ('ANI', 'Align_fraction_ref', 'Align_fraction_query')]
        if any(not math.isfinite(v) or not 0 <= v <= 100 for v in values):
            raise ValueError('Invalid ANI/aligned fraction')
        hit = AniHit(acc, taxid, name, row['Ref_file'], *values)
        hits.append(hit)
        if callable_name:
            callable_accessions.add(acc)
        audit.append({**asdict(hit), 'callable': callable_name,
            'ncbi_organism_name': metadata.get('ncbi_organism_name', ''),
            'gtdb_taxonomy': metadata.get('gtdb_taxonomy', ''),
            'ncbi_genome_category': metadata.get('ncbi_genome_category', '')})
    # Preserve earlier close competitors, including clinically important child
    # taxids. Searching a broader panel must not erase prior ambiguity.
    for value in evidence.get('closest_species', []):
        hit = AniHit(**{k: value[k] for k in AniHit.__dataclass_fields__ if k in value})
        hits.append(hit)
    decision = decide_identification(sample_id=sample, ani_hits=hits,
        expected_organism=expected, policy=policy).to_dict()
    if decision['status'] == 'confirmed_species':
        best = decision['closest_species'][0]
        if best['accession'] not in callable_accessions:
            decision.update(status='no_confident_identification', match_percent=None,
                reported_organism='No confident identification',
                reason='The strongest ANI reference lacks a callable clinical taxonomy.')
    # Never erase pre-existing isolate-contamination / secondary-signal concerns.
    if evidence.get('status') not in {None, '', 'no_confident_identification', 'ambiguous_species'}:
        decision.update(status='no_confident_identification', match_percent=None,
            reported_organism='No confident identification', reason='Prior identification requires review.')
    decision['method'] = 'skani whole-assembly search against GTDB representatives'
    decision['reference_audit'] = audit
    decision['original_identification'] = evidence
    return decision


def screen_assembly(*, assembly, cache_root, reference_root, metadata, evidence=None,
                    sample='query', expected='Unknown', threads=8,
                    command_prefix=('conda', 'run', '--no-capture-output', '-n', 'gtdbtk', 'skani'),
                    heartbeat=lambda: None, inherited_lock_fd=None):
    assembly = Path(assembly).resolve(strict=True)
    reference = Path(reference_root).resolve(strict=True)
    database = reference / 'skani/database'
    package = reference.parent / 'MD5SUM.txt'
    for name in ('metadata.tsv', 'index.db', 'sketches.db', 'markers.bin'):
        if not (database / name).is_file():
            raise ValueError(f'Missing GTDB ANI database component: {name}')
    version = subprocess.check_output([*command_prefix, '--version'], text=True, timeout=60).strip()
    identity = {'policy': POLICY, 'implementation': sha256(__file__),
        'decision_implementation': sha256(Path(__file__).with_name('species_identification.py')),
        'package_sha256': sha256(package), 'metadata_sha256': sha256(metadata),
        'sketch_metadata_sha256': sha256(database / 'metadata.tsv'), 'skani': version,
        'release': (reference / 'metadata/metadata.txt').read_text(),
        'assembly_sha256': sha256(assembly), 'evidence': evidence or {},
        'sample': sample, 'expected': expected}
    cache = Path(cache_root).resolve() / 'ani'
    cache.mkdir(parents=True, exist_ok=True)
    directory = cache / key(identity)
    with (cache / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        heartbeat()
        receipt_path = directory / 'receipt.json'
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt.get('identity') == identity and all(
                name in {'ani.tsv', 'confirmation.json'} and sha256(directory / name) == digest
                for name, digest in receipt.get('checksums', {}).items()
            ) and set(receipt.get('checksums', {})) == {'ani.tsv', 'confirmation.json'}:
                return {**receipt, 'directory': str(directory), 'reuse': 'assembly'}
        directory.mkdir(parents=True, exist_ok=True)
        output = directory / 'ani.tsv'
        execute([*command_prefix, 'search', '-d', str(database), '-q', str(assembly),
                 # Retain weak comparisons for supporting figures. The independent
                 # New identification uses ANI >=95% and either AF >=50%;
                 # saved historical policies retain their original thresholds.
                 '-o', str(output), '-t', str(max(1, min(16, int(threads)))), '--min-af', '1'],
                cache / 'screen.log', os.environ.copy(), heartbeat, lock.fileno(), inherited_lock_fd)
        with output.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'Ref_file', 'Query_file', 'ANI', 'Align_fraction_ref', 'Align_fraction_query'}.issubset(reader.fieldnames or []):
                raise ValueError('Malformed skani search output')
            rows = list(reader)
        if any(Path(row['Query_file']).resolve() != assembly for row in rows):
            raise ValueError('ANI results do not belong to the current assembly')
        labels = reference_labels(metadata, {accession(row['Ref_file']) for row in rows}) if rows else {}
        payload = evaluate_ani(rows, labels, evidence=evidence, sample=sample, expected=expected)
        if sha256(assembly) != identity['assembly_sha256']:
            raise ValueError('Assembly changed during ANI confirmation')
        payload['provenance'] = identity
        atomic_json(directory / 'confirmation.json', payload)
        receipt = {'identity': identity, 'status': 'identified' if payload['status'] == 'confirmed_species' else 'unresolved',
            'checksums': {name: sha256(directory / name) for name in ('ani.tsv', 'confirmation.json')}}
        atomic_json(receipt_path, receipt)
        return {**receipt, 'directory': str(directory), 'reuse': 'none'}


def investigate_assembly(*, evidence=None, sample='query', expected='Unknown', **kwargs):
    """ANI first for both platforms. Technical screen failure never starts a tree."""
    screen = screen_assembly(evidence=evidence, sample=sample, expected=expected, **kwargs)
    identification = json.loads((Path(screen['directory']) / 'confirmation.json').read_text())
    if screen['status'] == 'identified':
        return {**screen, 'payload': {'status': 'identified', 'identification': identification}}
    from .placement_visuals import ani_rows
    preferred = [row['accession'] for row in ani_rows(identification)[0][:5] if row.get('accession')]
    placement = place_iteratively(**kwargs, preferred_accessions=preferred)
    placement['payload'] = json.loads((Path(placement['directory']) / 'phylogenetic_placement.json').read_text())
    placement['confirmation'] = screen
    placement['payload']['identification'] = identification
    return placement
