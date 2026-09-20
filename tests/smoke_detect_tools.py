"""Opt-in tiny real-tool check, using a database gene, never a diagnostic sample.

Run in the BactiPipe image with BACTIPIPE_DB_DIR mounted read-only. All output
is temporary and removed after verification. No database updates/downloads.
"""
import json
import os
from pathlib import Path
import tempfile
from bactipipe.scripts.run_traits import main


def smoke():
    db = Path(os.environ["BACTIPIPE_DB_DIR"])
    with tempfile.TemporaryDirectory(prefix="detect-audit-") as tmp:
        root = Path(tmp)
        # A named, mapped Enterococcus positive ensures the chosen VF database
        # actually runs; an empty/random FASTA could only test negative output.
        sequence = []
        with (db / "virulencefinder_db/virulence_ent.fsa").open() as handle:
            collecting = False
            for line in handle:
                if line.startswith(">"):
                    if collecting:
                        break
                    collecting = line[1:].split(":")[0] == "agg"
                    continue
                if collecting:
                    sequence.append(line.strip())
        assert sequence, "Enterococcus agg test sequence missing"
        (root / "S1.fasta").write_text(">synthetic_test\n" + "".join(sequence) + "\n")
        (root / "sheet.tsv").write_text("S1\tSynthetic control\ttest\n")
        main(["--analyst", "test", "--genomes-dir", tmp, "--sample-sheet", str(root / "sheet.tsv"),
              "--run-name", "synthetic", "--organism", "Enterococcus", "--db-dir", str(db),
              "--outdir", str(root / "out"), "--threads", "2", "--jobs", "1"])
        result = json.loads((root / "out/analysis_result.json").read_text())
        assert result["status"] == "completed"
        assert result["parameters"]["abricate_ncbi"] is True
        assert result["parameters"]["amrfinder"] is False
        assert result["parameters"]["virulencefinder_databases"] == "virulence_ent"
        assert "agg" in (root / "out/S1.vf.tsv").read_text()
        assert (root / "out/reports/Unknown_vir_amr_detect.pdf").stat().st_size > 100
        print("DETECT_REAL_TOOL_SMOKE_SUCCESS (synthetic input, TSV/PDF verified)")


if __name__ == "__main__": smoke()
