"""Offline, exact-name conventional MLST selection; configured workflows win."""
import csv
import hashlib
import json
from pathlib import Path

from bactipipe.scripts.caps import CAPS


def name_key(value):
    return ' '.join(str(value).replace('_', ' ').split()).casefold()


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def resolve_mlst(organism, database):
    """Never infer genus membership, trim subspecies, or choose a nearby species.

    Existing CAPS tokens preserve their historical scheme and failure behavior.
    New species mappings are emitted by the installer from provider metadata.
    """
    database = Path(database)
    configured = CAPS.get(str(organism).strip().lower())
    if configured is not None:
        return dict(status='configured', scheme=configured.get('mlst_scheme'),
                    organism=organism, selection='BactiPipe configured workflow')
    path = database / 'mlst_catalog.json'
    if not path.is_file():
        return dict(status='catalog_unavailable', scheme=None, organism=organism)
    catalog = json.loads(path.read_text())
    if catalog.get('schema_version') != 1:
        raise ValueError('Unsupported MLST catalogue schema')
    matches = [r for r in catalog['species'] if name_key(r['species']) == name_key(organism)]
    if len(matches) > 1:
        raise ValueError('Duplicate normalized species in MLST catalogue')
    if not matches:
        return dict(status='no_matching_scheme', scheme=None, organism=organism)
    row = matches[0]
    selected = row.get('preferred_scheme')
    evidence = dict(organism=organism, species=row['species'], candidates=row['schemes'],
                    catalog_sha256=file_sha256(path), selection=row.get('selection'))
    if not selected:
        return dict(evidence, status='ambiguous', scheme=None)
    entry = next((r for r in catalog['schemes'] if r['scheme'] == selected), None)
    if not entry or entry.get('status') != 'installed':
        return dict(evidence, status='scheme_unavailable', scheme=None, requested_scheme=selected)
    if not entry.get('files'):
        raise ValueError('MLST catalogue has no verified scheme files')
    for relative, expected in entry['files'].items():
        target = database / relative
        if not target.resolve().is_relative_to(database.resolve()):
            raise ValueError('Unsafe MLST catalogue path')
        if not target.is_file() or file_sha256(target) != expected:
            raise ValueError(f'MLST scheme {selected} is missing or has changed; verify the database')
    return dict(evidence, status='selected', scheme=selected, source=entry['source'],
                access_note=entry.get('access_note'))


def species_rows(entries, preferences=()):
    """Preferences contain exact species + authoritative scheme URL, not rank."""
    preferred = {name_key(r['species']): r['source'] for r in preferences}
    grouped = {}
    for entry in entries:
        for species in entry.get('species', []):
            group = grouped.setdefault(name_key(species), dict(species=species, entries={}))
            group['entries'][entry['scheme']] = entry
    rows = []
    for key, group in sorted(grouped.items()):
        candidates = list(group['entries'].values())
        chosen = [r for r in candidates if r['source'] == preferred.get(key)]
        if not chosen and key not in preferred and len(candidates) == 1:
            chosen = candidates
        rows.append(dict(species=group['species'], schemes=sorted(group['entries']),
                         preferred_scheme=chosen[0]['scheme'] if len(chosen) == 1 else None,
                         selection='explicit preference' if key in preferred else
                         'unique scheme' if len(candidates) == 1 else 'review required'))
    return rows


def read_preferences(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    keys = [name_key(r['species']) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate species in MLST preferences')
    return rows
