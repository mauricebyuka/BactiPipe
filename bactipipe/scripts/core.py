import os
import sys
import time
import json
import shlex
import shutil
import tarfile
import tempfile
import urllib.parse
import urllib.request
from shutil import which
import subprocess
import datetime as dt

def env_cmd(tool: str) -> list:
    """
    Return conda-prefixed command for each external tool.

    Environments:
      VIRAMR_ENV   (default: viramr)   -> amrfinder, virulencefinder
      ABRICATE_ENV (default: genepid)  -> abricate

    Notes:
      - VirulenceFinder is invoked as a Python module:
            python -m virulencefinder
    """
    viramr_env = os.environ.get("VIRAMR_ENV", "viramr")
    abricate_env = os.environ.get("ABRICATE_ENV", "genepid")

    t = tool.lower()

    if t.startswith("amrfinder"):
        return ["conda", "run", "-n", viramr_env, "amrfinder"]

    if t.startswith("virulencefinder"):
        # IMPORTANT: module invocation, not a script
        return ["conda", "run", "-n", viramr_env, "python", "-m", "virulencefinder"]

    if t.startswith("abricate"):
        return ["conda", "run", "-n", abricate_env, "abricate"]

    return [tool]


def set_tmp_env(tmp_dir: str) -> dict:
    """Return a copy of os.environ with TMPDIR/TEMP/TMP set to tmp_dir."""
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir
    env["TEMP"]  = tmp_dir
    env["TMP"]   = tmp_dir
    return env

############CLI TOOL MANAGEMENT HELPERS ############
#Package management
ENV_NAMES = ("bactipipe", "genepid", "viramr")

def _run(cmd, **kwargs):
    """Run a shell command with pretty printing and error bubbling."""
    if isinstance(cmd, str):
        printable = cmd
        cmd = shlex.split(cmd)
    else:
        printable = " ".join(cmd)
    print(f"[bactipipe:update] $ {printable}", flush=True)
    subprocess.run(cmd, check=True, **kwargs)

def _which(cmd):
    from shutil import which
    return which(cmd)

def update_packages():
    """Update Python + conda packages in all envs, then pip-upgrade outdated pip packages."""
    mgr = "mamba" if _which("mamba") else "conda"

    # 1) Update all conda-managed packages (includes python) per env
    for env in ENV_NAMES:
        print(f"\n=== Updating conda packages in env: {env} ===", flush=True)
        _run([mgr, "update", "-n", env, "--all", "-y"])

    # 2) Upgrade pip-managed packages that are outdated in each env
    for env in ENV_NAMES:
        print(f"\n=== Upgrading pip packages in env: {env} ===", flush=True)
        # list outdated pip packages (name==version pins are not upgraded unless compatible)
        try:
            # capture the list of outdated in 'pkg==ver' lines
            result = subprocess.run(
                ["conda", "run", "-n", env, "python", "-m", "pip", "list", "--outdated", "--format=freeze"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            # pip may not exist in some envs (e.g., viramr if no pip installs) -> skip
            if "No module named pip" in e.stderr or "No module named pip" in e.stdout:
                print(f"[bactipipe:update] pip not present in env {env}; skipping pip upgrades.")
                continue
            # other errors bubble up
            raise

        lines = [ln.strip() for ln in result.stdout.splitlines() if ln.strip()]
        if not lines:
            print("[bactipipe:update] No outdated pip packages.")
            continue

        # Each line looks like: 'package==current_version'
        pkgs = [ln.split("==", 1)[0] for ln in lines if "==" in ln]
        for pkg in pkgs:
            _run(["conda", "run", "-n", env, "python", "-m", "pip", "install", "-U", pkg])

    print("\n✅ Update complete. Tip: if you maintain lockfiles, regenerate them now.", flush=True)


# DB root: 
DB_ROOT = os.environ.get("BACTIPIPE_DB_DIR", os.path.expanduser("~/src/database/BactiPipe"))

DB_SPECS = {
    "cgmlstfinder": {
        "label": "cgmlstfinder_db",
        "path": lambda: os.path.join(DB_ROOT, "cgmlstfinder_db"),
        "type": "git+INSTALL.py",
    },
    "kmerfinder": {
        "label": "kmerfinder_db",
        "path": lambda: os.path.join(DB_ROOT, "kmerfinder_db"),
        "type": "git+INSTALL.sh",
    },
    "mlst": {
        "label": "mlst_db",
        "path": lambda: os.path.join(DB_ROOT, "mlst_db"),
        "type": "BIGSdb REST",
    },
    "serotypefinder": {
        "label": "serotypefinder_db",
        "path": lambda: os.path.join(DB_ROOT, "serotypefinder_db"),
        "type": "git-only",
    },
    "virulencefinder": {
        "label": "virulencefinder_db",
        "path": lambda: os.path.join(DB_ROOT, "virulencefinder_db"),
        "type": "git-only",
    },
    # NOTE: tool-managed inside envs, no DB_ROOT path
    "amrfinder": {
        "label": "amrfinder",
        "path": lambda: None,
        "type": "tool-managed",
        "env": "viramr",
    },
    "abricate": {
        "label": "abricate",
        "path": lambda: None,
        "type": "tool-managed",
        "env": "genepid",
    },
}

# The original CGE ``mlst_db`` Bitbucket repository was removed in 2026.
# These endpoints expose the same profile/allele data through BIGSdb's
# documented REST API, allowing us to retain the existing CGE mlst.py
# algorithm and its historical BactiPipe scheme names.
MLST_SCHEMES = (
    ("abaumannii", "Acinetobacter baumannii", "https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1"),
    ("cperfringens", "Clostridium perfringens", "https://rest.pubmlst.org/db/pubmlst_cperfringens_seqdef/schemes/1"),
    ("cfetus", "Campylobacter fetus", "https://rest.pubmlst.org/db/pubmlst_campylobacter_nonjejuni_seqdef/schemes/9"),
    ("cjejuni", "Campylobacter jejuni/coli", "https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef/schemes/1"),
    ("ecoli", "Escherichia coli", "https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/1"),
    ("kpneumoniae", "Klebsiella pneumoniae species complex", "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/1"),
    ("lmonocytogenes", "Listeria monocytogenes", "https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2"),
    ("pmultocida", "Pasteurella multocida", "https://rest.pubmlst.org/db/pubmlst_pmultocida_seqdef/schemes/1"),
    ("senterica", "Salmonella enterica", "https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/2"),
    ("sepidermidis", "Staphylococcus epidermidis", "https://rest.pubmlst.org/db/pubmlst_sepidermidis_seqdef/schemes/1"),
    ("saureus", "Staphylococcus aureus", "https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/1"),
    ("spseudintermedius", "Staphylococcus pseudintermedius", "https://rest.pubmlst.org/db/pubmlst_spseudintermedius_seqdef/schemes/1"),
    ("szooepidemicus", "Streptococcus equi subsp. zooepidemicus", "https://rest.pubmlst.org/db/pubmlst_szooepidemicus_seqdef/schemes/1"),
    ("orhinotracheale", "Ornithobacterium rhinotracheale", "https://rest.pubmlst.org/db/pubmlst_orhinotracheale_seqdef/schemes/1"),
)

def _run(cmd, cwd=None):
    """Run a shell command with printing and error propagation."""
    if isinstance(cmd, str):
        printable = cmd
        cmd = shlex.split(cmd)
    else:
        printable = " ".join(cmd)
    print(f"[bactipipe] $ {printable}", flush=True)
    subprocess.run(cmd, check=True, cwd=cwd)
def list_databases():
    """List known databases and whether they appear to be installed."""
    print(f"\nDatabase root: {DB_ROOT}\n")

    rows = []
    abricate_dbs = []   # we'll fill this from `abricate --list`

    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()
        status = "unknown"
        location = ""

        # --- CGE git-based DBs ---
        if db_type.startswith("git") or db_type == "BIGSdb REST":
            if path and os.path.isdir(path) and os.listdir(path):
                status = "installed"
            else:
                status = "missing"
            location = path

        # --- AMRFinder (tool-managed in viramr env) ---
        elif key == "amrfinder":
            env_name = spec.get("env", "viramr")
            try:
                # Check that amrfinder runs; version output may include DB info.
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "amrfinder", "-V"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                status = "installed"
                db_dir_line = next(
                    (ln for ln in result.stdout.splitlines()
                     if "Database directory:" in ln),
                    None,
                )
                if db_dir_line:
                    location = db_dir_line.split(":", 1)[1].strip()
                else:
                    location = f"(inside {env_name} environment)"
            except Exception:
                status = "missing"
                location = f"(inside {env_name} environment)"

        # --- ABRicate (tool-managed in genepid env) ---
        elif key == "abricate":
            env_name = spec.get("env", "genepid")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "abricate", "--list"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
                status = "installed" if len(lines) > 1 else "missing"
                location = f"(inside {env_name} environment)"

                # Parse the table: header + rows like:
                # DATABASE    SEQUENCES   DBTYPE   DATE
                if len(lines) > 1:
                    header = lines[0].split()
                    for ln in lines[1:]:
                        parts = ln.split()
                        if len(parts) < 4:
                            continue
                        dbname = parts[0]
                        sequences = parts[1]
                        dbtype = parts[2]
                        date = " ".join(parts[3:])
                        abricate_dbs.append(
                            {
                                "name": dbname,
                                "sequences": sequences,
                                "dbtype": dbtype,
                                "date": date,
                            }
                        )
            except Exception:
                status = "missing"
                location = f"(inside {env_name} environment)"

        else:
            location = path or ""

        rows.append((label, key, db_type, status, location or "-"))

    # Main summary table
    headers = ("Label", "Key", "Type", "Status", "Location")
    print("{:<18} {:<14} {:<15} {:<10} {}".format(*headers))
    print("-" * 80)
    for label, key, db_type, status, location in rows:
        print("{:<18} {:<14} {:<15} {:<10} {}".format(label, key, db_type, status, location))

    # Extra: ABRicate databases
    if abricate_dbs:
        print("\nABRicate databases (env: genepid):")
        print("  {:<15} {:>10}   {:<6}   {}".format("DATABASE", "SEQUENCES", "TYPE", "DATE"))
        print("  " + "-" * 55)
        for db in abricate_dbs:
            print(
                "  {:<15} {:>10}   {:<6}   {}".format(
                    db["name"],
                    db["sequences"],
                    db["dbtype"],
                    db["date"],
                )
            )
    
def _parse_db_names(args):
    if not args:
        raise SystemExit("ERROR: 'update-db' requires at least one database name or 'all'.")

    raw = " ".join(args)
    tokens = []
    for part in raw.replace(",", " ").split():
        tokens.append(part.strip())

    if not tokens:
        raise SystemExit("ERROR: No valid database names provided.")

    if any(t.lower() == "all" for t in tokens):
        return list(DB_SPECS.keys())

    # Allow either the "key" or "<label>"/"<label>_db" forms
    normalized = []
    for t in tokens:
        t_low = t.lower().rstrip("/")
        t_low = t_low.replace("_db", "")  # virulencefinder_db → virulencefinder
        if t_low in DB_SPECS:
            normalized.append(t_low)
        else:
            # try matching by label
            for key, spec in DB_SPECS.items():
                if t_low == spec["label"].lower().replace("_db", ""):
                    normalized.append(key)
                    break
            else:
                raise SystemExit(f"ERROR: Unknown database name: {t}")

    # preserve order but drop duplicates
    seen = set()
    result = []
    for n in normalized:
        if n not in seen:
            seen.add(n)
            result.append(n)
    return result

def _ensure_git():
    if which("git") is None:
        raise SystemExit("ERROR: 'git' not found on PATH. Please install git and retry.")

def _download_and_extract_tar_gz(url: str, dest_dir: str, label: str) -> None:
    """
    Download a tar.gz database archive and extract it into dest_dir.

    Strategy:
      - Download to a temporary file
      - Extract to a temporary directory
      - If dest_dir exists and is non-empty, move it to a timestamped backup
      - Move extracted content into dest_dir

    Handles archives that either:
      A) contain a single top-level directory, or
      B) contain files directly at top-level.
    """
    os.makedirs(DB_ROOT, exist_ok=True)

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    backup_dir = f"{dest_dir}.bak-{timestamp}"

    with tempfile.TemporaryDirectory(prefix="bactipipe_db_") as tmpd:
        tgz_path = os.path.join(tmpd, f"{label}.tar.gz")

        print(f"[bactipipe] Downloading {label} from: {url}")
        urllib.request.urlretrieve(url, tgz_path)

        extract_dir = os.path.join(tmpd, "extract")
        os.makedirs(extract_dir, exist_ok=True)

        print(f"[bactipipe] Extracting {label} archive…")
        with tarfile.open(tgz_path, "r:gz") as tf:
            tf.extractall(path=extract_dir)

        # Determine extracted payload location
        entries = [e for e in os.listdir(extract_dir) if e not in (".", "..")]
        if not entries:
            raise RuntimeError(f"{label}: archive extracted but no files were found.")

        # If the tarball contains exactly one top-level directory, treat that as payload
        payload_path = extract_dir
        if len(entries) == 1:
            single = os.path.join(extract_dir, entries[0])
            if os.path.isdir(single):
                payload_path = single

        # Backup existing DB if present and non-empty
        if os.path.isdir(dest_dir) and os.listdir(dest_dir):
            print(f"[bactipipe] Backing up existing {label} to: {backup_dir}")
            shutil.move(dest_dir, backup_dir)
        else:
            # Ensure clean destination
            if os.path.exists(dest_dir):
                shutil.rmtree(dest_dir, ignore_errors=True)

        # Install new DB into dest_dir
        print(f"[bactipipe] Installing {label} into: {dest_dir}")
        shutil.move(payload_path, dest_dir)

        print(f"[bactipipe] {label} install complete.")

def _update_cgmlstfinder():
    db_dir = DB_SPECS["cgmlstfinder"]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)

    url = "https://cge.food.dtu.dk/services/cgMLSTFinder/etc/cgmlst_db.tar.gz"
    _download_and_extract_tar_gz(url=url, dest_dir=db_dir, label="cgmlstfinder_db")


def _update_kmerfinder():
    db_dir = DB_SPECS["kmerfinder"]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)

    url = "https://cge.food.dtu.dk/services/KmerFinder/etc/kmerfinder_db.tar.gz"
    _download_and_extract_tar_gz(url=url, dest_dir=db_dir, label="kmerfinder_db")


def _download_url(url: str, *, attempts: int = 3) -> bytes:
    """Download a BIGSdb resource with bounded retries and a useful agent."""
    headers = {"User-Agent": "BactiPipe database updater"}
    hostname = urllib.parse.urlparse(url).hostname
    authorization_variable = {
        "rest.pubmlst.org": "BACTIPIPE_PUBMLST_AUTHORIZATION",
        "bigsdb.pasteur.fr": "BACTIPIPE_PASTEUR_AUTHORIZATION",
    }.get(hostname)
    authorization = (
        os.environ.get(authorization_variable, "").strip()
        if authorization_variable
        else ""
    )
    if authorization:
        headers["Authorization"] = authorization
    request = urllib.request.Request(
        url,
        headers=headers,
    )
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception as exc:
            last_error = exc
            if attempt < attempts:
                print(
                    f"[bactipipe] Download attempt {attempt} failed for {url}; retrying…",
                    flush=True,
                )
                time.sleep(attempt)
    raise RuntimeError(f"Unable to download {url}: {last_error}") from last_error


def _download_legacy_mlst_scheme(
    destination: str,
    scheme_name: str,
    organism_name: str,
    scheme_url: str,
) -> dict:
    """Create one CGE mlst.py-compatible scheme from a BIGSdb scheme."""
    try:
        scheme_info = json.loads(_download_url(scheme_url).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Invalid BIGSdb scheme response for {scheme_name}") from exc

    loci_urls = scheme_info.get("loci") or []
    profiles_url = scheme_info.get("profiles_csv")
    if not loci_urls or not profiles_url:
        raise RuntimeError(
            f"BIGSdb scheme {scheme_name} has no loci or profile download URL"
        )

    loci = [
        urllib.parse.unquote(url.rstrip("/").rsplit("/", 1)[-1])
        for url in loci_urls
    ]
    profiles = _download_url(profiles_url).decode("utf-8-sig")
    profile_lines = profiles.splitlines()
    if not profile_lines:
        raise RuntimeError(f"BIGSdb scheme {scheme_name} returned no profiles")
    header = profile_lines[0].split("\t")
    if not header or header[0].strip().lower() not in {"st", "sequence_type"}:
        raise RuntimeError(f"BIGSdb scheme {scheme_name} has no ST profile column")
    if header[1 : len(loci) + 1] != loci:
        raise RuntimeError(
            f"BIGSdb scheme {scheme_name} profile loci do not match its scheme loci"
        )

    scheme_dir = os.path.join(destination, scheme_name)
    os.makedirs(scheme_dir)
    profile_path = os.path.join(scheme_dir, f"{scheme_name}.tsv")
    with open(profile_path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(profile_lines) + "\n")

    fasta_path = os.path.join(scheme_dir, f"{scheme_name}.fsa")
    with open(fasta_path, "wb") as handle:
        for locus, locus_url in zip(loci, loci_urls):
            allele_data = _download_url(f"{locus_url}/alleles_fasta")
            if not allele_data.lstrip().startswith(b">"):
                raise RuntimeError(
                    f"BIGSdb locus {locus} returned no FASTA sequences"
                )
            handle.write(allele_data.rstrip() + b"\n")

    return {
        "scheme": scheme_name,
        "organism": organism_name,
        "source": scheme_url,
        "description": scheme_info.get("description"),
        "last_updated": scheme_info.get("last_updated"),
        "records": scheme_info.get("records") or scheme_info.get("profile_count"),
        "loci": loci,
        "access_note": scheme_info.get("message"),
    }


def _update_mlst():
    """Rebuild the legacy CGE MLST layout from authoritative BIGSdb APIs."""
    db_dir = DB_SPECS["mlst"]["path"]()
    if os.path.isfile(os.path.join(db_dir, 'mlst_catalog.json')):
        raise RuntimeError('This database has an expanded MLST catalogue. Use update-db mlst '
                           '--all-schemes --seed ACTIVE_DB --output NEW_VERSION --work BUILD_DIR; '
                           'the configured-only updater must not remove installed schemes.')
    os.makedirs(DB_ROOT, exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    backup_dir = f"{db_dir}.bak-{timestamp}"

    with tempfile.TemporaryDirectory(prefix="bactipipe_mlst_", dir=DB_ROOT) as tmpd:
        payload = os.path.join(tmpd, "mlst_db")
        os.makedirs(payload)
        config_rows = [
            "# Database configuration file generated from BIGSdb REST APIs",
            "# species_db\tspecies_name\tfolder_content(locilist,profilefile)",
        ]
        source_rows = []

        for scheme_name, organism_name, scheme_url in MLST_SCHEMES:
            print(
                f"[bactipipe] Downloading MLST scheme: {scheme_name}",
                flush=True,
            )
            source = _download_legacy_mlst_scheme(
                payload, scheme_name, organism_name, scheme_url
            )
            config_rows.append(
                f"{scheme_name}\t{organism_name}\t{','.join(source['loci'])}"
            )
            source_rows.append(source)

        with open(os.path.join(payload, "config"), "w", encoding="utf-8") as handle:
            handle.write("\n".join(config_rows) + "\n")
        with open(
            os.path.join(payload, "bactipipe_mlst_sources.json"),
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                {
                    "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "format": "CGE mlst.py legacy database",
                    "schemes": source_rows,
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        if os.path.isdir(db_dir) and os.listdir(db_dir):
            print(f"[bactipipe] Backing up existing mlst_db to: {backup_dir}")
            shutil.move(db_dir, backup_dir)
        elif os.path.exists(db_dir):
            shutil.rmtree(db_dir)

        try:
            shutil.move(payload, db_dir)
        except Exception:
            if os.path.isdir(backup_dir) and not os.path.exists(db_dir):
                shutil.move(backup_dir, db_dir)
            raise

    print(f"[bactipipe] MLST database install complete: {db_dir}")

def _update_simple_git_db(key, url):
    _ensure_git()
    db_dir = DB_SPECS[key]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)
    if os.path.isdir(os.path.join(db_dir, ".git")):
        print(f"[bactipipe] Updating {DB_SPECS[key]['label']} in {db_dir}")
        _run(["git", "pull", "--ff-only"], cwd=db_dir)
    else:
        print(f"[bactipipe] Cloning {DB_SPECS[key]['label']} into {db_dir}")
        _run(["git", "clone", url, db_dir])

def update_databases(args):
    dbs = _parse_db_names(args)
    print(f"\nDatabase root: {DB_ROOT}")
    os.makedirs(DB_ROOT, exist_ok=True)

    for db in dbs:
        print(f"\n=== Updating database: {db} ===")
        if db == "cgmlstfinder":
            _update_cgmlstfinder()
        elif db == "kmerfinder":
            _update_kmerfinder()
        elif db == "mlst":
            _update_mlst()
        elif db == "serotypefinder":
            _update_simple_git_db("serotypefinder", "https://bitbucket.org/genomicepidemiology/serotypefinder_db.git")
        elif db == "virulencefinder":
            _update_simple_git_db("virulencefinder", "https://bitbucket.org/genomicepidemiology/virulencefinder_db.git")
        elif db == "amrfinder":
            print("[bactipipe] Updating AMRFinderPlus database (viramr env)…")
            _run(["conda", "run", "-n", "viramr", "amrfinder", "-u"])
        elif db == "abricate":
            print("[bactipipe] Updating ABRicate databases (genepid env)…")
            _run(["conda", "run", "-n", "genepid", "abricate", "--setupdb"])
        else:
            print(f"[bactipipe] WARNING: unknown db key {db}, skipping.")

    print("\n✅ Database update complete.")

def check_updates():
    """Check for environment and database updates without installing them."""

    import datetime as dt

    print("\n=== Checking Conda environment updates ===\n")

    # -----------------------------------------
    # Check outdated conda packages per env
    # -----------------------------------------
    for env in ENV_NAMES:
        try:
            result = subprocess.run(
                ["conda", "search", "--outdated", "--name", env],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            lines = [ln for ln in result.stdout.splitlines() if ln.strip()]

            if len(lines) > 1:
                # First line is header → skip it
                outdated = []
                for ln in lines[1:]:
                    parts = ln.split()
                    if len(parts) > 0:
                        outdated.append(parts[0])
                if outdated:
                    pkgs = ", ".join(outdated)
                    print(f"[{env}] {len(outdated)} updates available: {pkgs}")
                else:
                    print(f"[{env}] up-to-date")
            else:
                print(f"[{env}] up-to-date")

        except Exception as e:
            print(f"[{env}] unable to check (error: {e})")


    print("\n=== Checking databases ===\n")

    today = dt.date.today()
    stale_days = 183  # ~6 months

    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()

        # -------------------------
        # CGE git-based DBs
        # -------------------------
        if db_type == "git-only" or db_type.startswith("git"):
            if path and os.path.isdir(os.path.join(path, ".git")):
                try:
                    subprocess.run(
                        ["git", "fetch"],
                        cwd=path,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=False,
                        text=True,
                    )
                    local = subprocess.check_output(
                        ["git", "rev-parse", "HEAD"], cwd=path, text=True
                    ).strip()
                    remote = subprocess.check_output(
                        ["git", "rev-parse", "@{u}"], cwd=path, text=True
                    ).strip()
                    if local != remote:
                        print(f"{label:18} update available")
                    else:
                        print(f"{label:18} up-to-date")
                except Exception:
                    print(f"{label:18} cannot determine (error)")
            else:
                print(f"{label:18} not installed")

        # -------------------------
        # AMRFinderPlus DB (tool-managed in viramr env)
        # SAFE: show installed DB version, do NOT update
        # -------------------------
        elif key == "amrfinder":
            env_name = spec.get("env", "viramr")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "amrfinder", "--database_version"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                db_version = None
                db_dir = None
                for ln in result.stdout.splitlines():
                    if "Database version:" in ln:
                        db_version = ln.split(":", 1)[1].strip()
                    if "Database directory:" in ln:
                        db_dir = ln.split(":", 1)[1].strip().strip("'\"")

                if db_version:
                    msg = f"installed DB version {db_version}"
                    if db_dir:
                        msg += f" at {db_dir}"
                    msg += " (not checked against NCBI; run 'bactipipe update-db amrfinder' to refresh)"
                else:
                    msg = "installed (database version unknown)"

                print(f"{label:18} {msg}")
            except Exception as e:
                print(f"{label:18} unable to check (error: {e})")

        # -------------------------
        # ABRicate DBs (tool-managed in genepid env)
        # classify by age based on DATE column
        # -------------------------
        elif key == "abricate":
            env_name = spec.get("env", "genepid")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "abricate", "--list"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
                if len(lines) <= 1:
                    print(f"{label:18} no databases installed")
                    continue

                stale = 0
                total = 0
                oldest = None
                newest = None

                # Expected format:
                # DATABASE    SEQUENCES   DBTYPE   DATE
                for ln in lines[1:]:
                    parts = ln.split()
                    if len(parts) < 4:
                        continue
                    date_str = parts[3]
                    total += 1
                    try:
                        db_date = dt.datetime.strptime(date_str, "%Y-%b-%d").date()
                    except ValueError:
                        # skip unparsable dates
                        continue

                    if oldest is None or db_date < oldest:
                        oldest = db_date
                    if newest is None or db_date > newest:
                        newest = db_date

                    if (today - db_date).days > stale_days:
                        stale += 1

                if total == 0 or oldest is None:
                    print(f"{label:18} installed (dates unavailable)")
                else:
                    if stale == 0:
                        print(
                            f"{label:18} all {total} DBs downloaded within 6 months "
                            f"(newest: {newest.isoformat()})"
                        )
                    else:
                        print(
                            f"{label:18} {stale}/{total} DBs older than 6 months "
                            f"(oldest: {oldest.isoformat()})"
                        )
            except Exception as e:
                print(f"{label:18} unable to check (error: {e})")

        else:
            # For any future DB types not explicitly handled
            print(f"{label:18} check not implemented")

    print("\nDone.\n")


