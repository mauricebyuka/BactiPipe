"""Resumable conventional MLST catalogue builder. Never edits an active database."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from bactipipe.mlst_catalog import file_sha256, species_rows, read_preferences
from bactipipe.scripts.core import MLST_SCHEMES, _download_legacy_mlst_scheme

XML_URL = 'https://pubmlst.org/static/data/dbases.xml'
API_URL = 'https://rest.pubmlst.org/db'
ALLOWED_HOSTS = {'pubmlst.org', 'rest.pubmlst.org', 'bigsdb.pasteur.fr'}


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def checked_url(url):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != 'https' or parsed.hostname not in ALLOWED_HOSTS:
        raise ValueError('MLST resource is not an approved HTTPS provider')
    return url


class Downloads:
    """Checksum-bound, per-build cache. Failed responses never become receipts."""
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def __call__(self, url):
        checked_url(url)
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.root / key
        receipt = self.root / (key + '.json')
        if path.is_file() and receipt.is_file():
            meta = json.loads(receipt.read_text())
            if meta['url'] == url and meta['sha256'] == file_sha256(path):
                return path.read_bytes()
        # Public access only; no credential discovery or privilege escalation.
        request = urllib.request.Request(url, headers={'User-Agent': 'BactiPipe MLST catalogue installer'})
        with urllib.request.urlopen(request, timeout=45) as response:
            checked_url(response.url)
            data = response.read(128 * 1024 * 1024 + 1)
        if len(data) > 128 * 1024 * 1024:
            raise ValueError('MLST resource exceeds the conventional-scheme size limit')
        tmp = path.with_suffix('.tmp')
        tmp.write_bytes(data)
        tmp.replace(path)
        atomic_json(receipt, dict(url=url, sha256=file_sha256(path)))
        return data


def discover(fetch):
    """Union classic curated listings with current API MLST schemes.

    Do not enumerate/download wgMLST/cgMLST allele collections. Ambiguous names
    are retained as catalogue entries, never converted into guessed species.
    """
    candidates = {}
    def add(url, label='', configured=None):
        checked_url(url)
        row = candidates.setdefault(url, dict(source=url, labels=[]))
        if label and label not in row['labels']:
            row['labels'].append(label)
        if configured:
            row['configured'] = configured
    for element in ET.fromstring(fetch(XML_URL)).findall('species'):
        url = element.findtext('./mlst/database/profiles/url', '').strip()
        if url.endswith('/profiles_csv'):
            add(url.removesuffix('/profiles_csv'), (element.text or '').strip())
    # Preserve the exact scheme IDs/URLs used by established BactiPipe workflows.
    for scheme, species, url in MLST_SCHEMES:
        add(url, species, scheme)
    errors = []
    providers = json.loads(fetch(API_URL))
    for provider in providers:
        for database in provider.get('databases', []):
            if not database['name'].endswith('_seqdef'):
                continue
            url = database['href'] + '/schemes'
            try:
                while url:
                    listing = json.loads(fetch(url))
                    for entry in listing.get('schemes', []):
                        description = entry.get('description', '')
                        if re.search(r'\bMLST\b|multilocus sequence typing', description, re.I) and not re.search(r'cgmlst|wgmlst|rmlst|core.genome|whole.genome', description, re.I):
                            add(entry['scheme'])
                    url = listing.get('paging', {}).get('next')
            except Exception as exc:
                errors.append(dict(source=url, error=str(exc)))
    return list(candidates.values()), errors


def describe(candidate, fetch):
    info = json.loads(fetch(candidate['source']))
    loci = info.get('loci') or []
    if not 2 <= len(loci) <= 20 or not info.get('profiles_csv'):
        raise ValueError('Not an eligible conventional MLST scheme (2–20 loci and ST profiles required)')
    if re.search(r'cgmlst|wgmlst|rmlst|core.genome|whole.genome', info.get('description', ''), re.I):
        raise ValueError('Not a conventional MLST scheme')
    species = {r['name'] for r in info.get('ncbi_taxon', []) if r.get('rank') in ('species', 'subspecies')}
    for label in candidate['labels']:
        label = re.sub(r'#\d+$', '', label).strip()
        # Only unambiguous binomials/subspecies; no genus-wide or complex guessing.
        if re.fullmatch(r'[A-Z][a-z]+ [a-z]+(?: subsp\. [a-z]+)?', label):
            species.add(label)
    parts = candidate['source'].split('/')
    default_id = parts[-3].removeprefix('pubmlst_').removesuffix('_seqdef') + '__' + parts[-1]
    # Host-qualified hash prevents identical provider database names colliding.
    scheme = candidate.get('configured') or default_id + '_' + hashlib.sha256(candidate['source'].encode()).hexdigest()[:8]
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', scheme):
        raise ValueError('Unsafe MLST scheme identifier')
    return dict(candidate, scheme=scheme, species=sorted(species), description=info.get('description'),
                access_note=info.get('message'), loci=[urllib.parse.unquote(u.rsplit('/', 1)[-1]) for u in loci])


def validate_scheme(root, entry):
    folder = Path(root) / entry['scheme']
    profile = folder / (entry['scheme'] + '.tsv')
    fasta = folder / (entry['scheme'] + '.fsa')
    with profile.open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        if header[0].casefold() not in ('st', 'sequence_type') or header[1:1+len(entry['loci'])] != entry['loci']:
            raise ValueError('MLST profile header does not match catalogue loci')
        if not next(reader, None):
            raise ValueError('MLST profile table is empty')
    headers = set()
    with fasta.open() as handle:
        for line in handle:
            if line.startswith('>'):
                headers.add(line[1:].split()[0].rsplit('_', 1)[0])
    if not set(entry['loci']).issubset(headers):
        raise ValueError('MLST allele FASTA is missing a required locus')
    return {p.relative_to(root).as_posix(): file_sha256(p) for p in (profile, fasta)}


def build(output, work, seed=None, preferences=None):
    output, work = Path(output).resolve(), Path(work).resolve()
    if output.exists():
        raise ValueError('Output already exists; use a new version directory')
    if output.is_relative_to(work) or work.is_relative_to(output):
        raise ValueError('Output and build workspace must be separate directories')
    if seed and (output.is_relative_to(Path(seed).resolve()) or work.is_relative_to(Path(seed).resolve())
                 or Path(seed).resolve().is_relative_to(work)):
        raise ValueError('Build must not overwrite the seed database')
    work.mkdir(parents=True, exist_ok=True)
    atomic_json(work / 'status.json', dict(status='DISCOVERING', output=str(output)))
    payload = work / 'payload'
    payload.mkdir(exist_ok=True)
    fetch = Downloads(work / 'downloads')
    # Existing converter is reused; its requests are checksum-cached for resume.
    from bactipipe.scripts import core
    original_fetch = core._download_url
    core._download_url = fetch
    try:
        inventory_path = work / 'inventory.json'
        if inventory_path.exists():
            inventory = json.loads(inventory_path.read_text())
        else:
            candidates, errors = discover(fetch)
            inventory = dict(candidates=candidates, discovery_errors=errors)
            atomic_json(inventory_path, inventory)
        atomic_json(work / 'status.json', dict(status='INSTALLING', total=len(inventory['candidates']), output=str(output)))
        entries = []
        seed_sources = {}
        if seed and (Path(seed) / 'bactipipe_mlst_sources.json').is_file():
            seed_sources = {r['source']: r for r in json.loads((Path(seed) / 'bactipipe_mlst_sources.json').read_text())['schemes']}
        for candidate in inventory['candidates']:
            entry = dict(candidate, status='unavailable')
            try:
                entry = describe(candidate, fetch)
                folder = payload / entry['scheme']
                receipt = work / (entry['scheme'] + '.receipt.json')
                verified = False
                if receipt.exists():
                    saved = json.loads(receipt.read_text())
                    try:
                        verified = saved['source'] == entry['source'] and saved['files'] == validate_scheme(payload, entry)
                    except (OSError, ValueError, StopIteration):
                        verified = False
                    if verified:
                        entry = saved
                if not verified:
                    if folder.exists():
                        # Interrupted payload only, never a source or installed DB.
                        folder.rename(work / (entry['scheme'] + '.incomplete-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')))
                    existing = seed_sources.get(entry['source'])
                    if existing and existing['scheme'] == entry['scheme']:
                        shutil.copytree(Path(seed) / entry['scheme'], folder)
                        entry['loci'] = existing['loci']
                        entry['access_note'] = existing.get('access_note')
                        entry['seed_preserved'] = True
                    else:
                        source = _download_legacy_mlst_scheme(str(payload), entry['scheme'],
                            ', '.join(entry['species']) or entry['description'], entry['source'])
                        entry.update(source)
                    entry['files'] = validate_scheme(payload, entry)
                    entry['status'] = 'installed'
                    atomic_json(receipt, entry)
                print(entry['scheme'], 'installed', flush=True)
            except Exception as exc:
                entry['status'] = 'unavailable'
                entry['error'] = str(exc)
                print(candidate['source'], 'unavailable:', str(exc), flush=True)
            entries.append(entry)
            atomic_json(work / 'progress.json', dict(processed=len(entries), total=len(inventory['candidates']), schemes=entries))
        preferences = read_preferences(preferences) if preferences else []
        # Preserve known choices for full species names, too.
        preferences += [dict(species=species, source=url) for _, species, url in MLST_SCHEMES
                        if species not in {r['species'] for r in preferences}]
        catalog = dict(schema_version=1, generated_at=datetime.now(timezone.utc).isoformat(),
                       schemes=entries, species=species_rows([e for e in entries if 'scheme' in e], preferences),
                       discovery_errors=inventory['discovery_errors'])
        atomic_json(payload / 'mlst_catalog.json', catalog)
        with (payload / 'species_to_mlst.tsv').open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['species', 'preferred_scheme', 'available_schemes', 'selection'])
            for row in catalog['species']:
                writer.writerow([row['species'], row['preferred_scheme'] or '', ','.join(row['schemes']), row['selection']])
        with (payload / 'config').open('w') as handle:
            handle.write('# species_db\tspecies_name\tfolder_content(locilist,profilefile)\n')
            for row in entries:
                if row['status'] == 'installed':
                    handle.write('\t'.join([row['scheme'], ', '.join(row['species']) or row['description'], ','.join(row['loci'])]) + '\n')
        atomic_json(payload / 'bactipipe_mlst_sources.json', dict(format='CGE mlst.py legacy database',
            generated_at=catalog['generated_at'], schemes=[e for e in entries if e['status'] == 'installed']))
        if not any(e['status'] == 'installed' for e in entries):
            raise ValueError('No MLST schemes installed')
        # Refuse a completed snapshot that loses any established scheme.
        if seed and any(e.get('configured') and e['status'] != 'installed' for e in entries):
            raise ValueError('Configured scheme unavailable; incomplete build retained for review')
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(payload, output)
        status = dict(status='READY_FOR_REVIEW', output=str(output), installed=sum(e['status']=='installed' for e in entries),
                      unavailable=sum(e['status']!='installed' for e in entries), discovery_errors=len(inventory['discovery_errors']))
        atomic_json(work / 'status.json', status)
        return status
    finally:
        core._download_url = original_fetch


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-schemes', action='store_true', help='Discover conventional MLST schemes (not cgMLST/wgMLST)')
    parser.add_argument('--output', type=Path, required=True, help='New immutable database directory; never the active path')
    parser.add_argument('--work', type=Path, required=True, help='Persistent resumable build workspace')
    parser.add_argument('--seed', type=Path, help='Preserve installed configured schemes from this database')
    parser.add_argument('--preferences', type=Path, help='TSV with species and authoritative source URL columns')
    args = parser.parse_args(argv)
    try:
        print(json.dumps(build(args.output, args.work, args.seed, args.preferences)))
    except Exception as exc:
        args.work.mkdir(parents=True, exist_ok=True)
        atomic_json(args.work / 'status.json', dict(status='FAILED', error=str(exc)))
        raise


if __name__ == '__main__':
    main()
