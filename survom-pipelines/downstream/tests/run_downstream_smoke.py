#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[3]
DOWNSTREAM = ROOT / "survom-pipelines" / "downstream"
FIXTURES = DOWNSTREAM / "tests" / "fixtures"
CLI = DOWNSTREAM / "bin" / "downstream_cli.py"


def run(*args: str) -> None:
    proc = subprocess.run([sys.executable, str(CLI), *args], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(args)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def main() -> int:
    with TemporaryDirectory(prefix="survom_downstream_smoke_") as tmp:
        root = Path(tmp)
        merge = root / "01_merge"
        validate = root / "02_validate"
        filt = root / "03_filter"
        norm = root / "04_norm"
        pca = root / "05_pca"
        umap = root / "06_umap"
        plsda = root / "08_plsda"
        run("merge-featurecounts", "--counts", str(FIXTURES / "featurecounts_raw.txt"), "--outdir", str(merge))
        run(
            "validate",
            "--counts",
            str(merge / "merged_raw_counts.csv"),
            "--metadata",
            str(FIXTURES / "sample_metadata.csv"),
            "--contrasts",
            str(FIXTURES / "contrasts.csv"),
            "--outdir",
            str(validate),
        )
        run("filter", "--counts", str(validate / "validated_counts.csv"), "--metadata", str(validate / "validated_metadata.csv"), "--outdir", str(filt))
        run("normalize", "--counts", str(filt / "filtered_counts.csv"), "--metadata", str(validate / "validated_metadata.csv"), "--outdir", str(norm))
        run("pca", "--expression", str(norm / "vst_expression.csv"), "--metadata", str(validate / "validated_metadata.csv"), "--outdir", str(pca))
        run("umap", "--expression", str(norm / "vst_expression.csv"), "--metadata", str(validate / "validated_metadata.csv"), "--outdir", str(umap))
        run("plsda", "--expression", str(norm / "vst_expression.csv"), "--metadata", str(validate / "validated_metadata.csv"), "--outdir", str(plsda))
        required = [
            merge / "merged_raw_counts.csv",
            validate / "validation_report.json",
            filt / "filtered_counts.csv",
            norm / "vst_expression.csv",
            pca / "pca_scores.csv",
            umap / "umap_coordinates.csv",
            plsda / "plsda_summary.json",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise RuntimeError(f"Missing downstream smoke outputs: {missing}")
        print(f"[downstream smoke ok] {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
