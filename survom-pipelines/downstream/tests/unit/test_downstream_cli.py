from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DOWNSTREAM = ROOT / "survom-pipelines" / "downstream"
FIXTURES = DOWNSTREAM / "tests" / "fixtures"
CLI = DOWNSTREAM / "bin" / "downstream_cli.py"


def run_cli(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(CLI), *args],
        cwd=cwd or DOWNSTREAM,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def test_merge_validate_filter_pca_umap_plsda(tmp_path):
    merge_dir = tmp_path / "01_merge"
    proc = run_cli("merge-featurecounts", "--counts", str(FIXTURES / "featurecounts_raw.txt"), "--outdir", str(merge_dir))
    assert proc.returncode == 0, proc.stderr
    assert (merge_dir / "merged_raw_counts.csv").exists()
    assert (merge_dir / "metadata.json").exists()
    assert (merge_dir / "checksums.sha256").exists()

    validate_dir = tmp_path / "02_validate"
    proc = run_cli(
        "validate",
        "--counts",
        str(merge_dir / "merged_raw_counts.csv"),
        "--metadata",
        str(FIXTURES / "sample_metadata.csv"),
        "--contrasts",
        str(FIXTURES / "contrasts.csv"),
        "--outdir",
        str(validate_dir),
    )
    assert proc.returncode == 0, proc.stderr
    report = json.loads((validate_dir / "validation_report.json").read_text())
    assert report["genes"] == 6
    assert report["samples"] == 6

    filter_dir = tmp_path / "03_filter"
    proc = run_cli(
        "filter",
        "--counts",
        str(validate_dir / "validated_counts.csv"),
        "--metadata",
        str(validate_dir / "validated_metadata.csv"),
        "--outdir",
        str(filter_dir),
    )
    assert proc.returncode == 0, proc.stderr
    filtered = pd.read_csv(filter_dir / "filtered_counts.csv")
    assert "gene_id" in filtered.columns

    norm_dir = tmp_path / "04_norm"
    proc = run_cli(
        "normalize",
        "--counts",
        str(filter_dir / "filtered_counts.csv"),
        "--metadata",
        str(validate_dir / "validated_metadata.csv"),
        "--outdir",
        str(norm_dir),
    )
    assert proc.returncode == 0, proc.stderr

    pca_a = tmp_path / "05_pca_a"
    pca_b = tmp_path / "05_pca_b"
    for outdir in (pca_a, pca_b):
        proc = run_cli(
            "pca",
            "--expression",
            str(norm_dir / "vst_expression.csv"),
            "--metadata",
            str(validate_dir / "validated_metadata.csv"),
            "--outdir",
            str(outdir),
        )
        assert proc.returncode == 0, proc.stderr
    assert (pca_a / "pca_scores.csv").read_text() == (pca_b / "pca_scores.csv").read_text()

    umap_dir = tmp_path / "06_umap"
    proc = run_cli(
        "umap",
        "--expression",
        str(norm_dir / "vst_expression.csv"),
        "--metadata",
        str(validate_dir / "validated_metadata.csv"),
        "--outdir",
        str(umap_dir),
    )
    assert proc.returncode == 0, proc.stderr
    assert (umap_dir / "umap_coordinates.csv").exists()

    plsda_dir = tmp_path / "08_plsda"
    proc = run_cli(
        "plsda",
        "--expression",
        str(norm_dir / "vst_expression.csv"),
        "--metadata",
        str(validate_dir / "validated_metadata.csv"),
        "--outdir",
        str(plsda_dir),
    )
    assert proc.returncode == 0, proc.stderr
    assert (plsda_dir / "plsda_summary.json").exists()


def test_invalid_counts_fail_clearly(tmp_path):
    bad = tmp_path / "bad_counts.csv"
    bad.write_text("gene_id,S1\nGENE1,-1\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("sample_id,condition\nS1,control\n")
    proc = run_cli("validate", "--counts", str(bad), "--metadata", str(metadata), "--outdir", str(tmp_path / "out"))
    assert proc.returncode != 0
    assert "negative counts" in proc.stderr


def test_enrichment_gsea_and_state(tmp_path):
    significant = tmp_path / "significant_genes.csv"
    significant.write_text("gene_id,log2FoldChange\nGENE001,2\nGENE003,-2\n")
    ranked = tmp_path / "ranked_genes.csv"
    ranked.write_text("gene_id,stat\nGENE001,5\nGENE002,4\nGENE003,-4\nGENE004,-3\n")
    universe = tmp_path / "universe.csv"
    universe.write_text("gene_id\nGENE001\nGENE002\nGENE003\nGENE004\nGENE005\n")
    go_dir = tmp_path / "10_go"
    proc = run_cli(
        "go-enrichment",
        "--significant-genes",
        str(significant),
        "--ranked-genes",
        str(ranked),
        "--universe",
        str(universe),
        "--gmt",
        str(FIXTURES / "mini_sets.gmt"),
        "--outdir",
        str(go_dir),
    )
    assert proc.returncode == 0, proc.stderr
    assert (go_dir / "go_overrepresentation_up.csv").exists()

    gsea_dir = tmp_path / "12_gsea"
    proc = run_cli("gsea", "--ranked-genes", str(ranked), "--gmt", str(FIXTURES / "mini_sets.gmt"), "--outdir", str(gsea_dir))
    assert proc.returncode == 0, proc.stderr
    assert (gsea_dir / "gsea_results.csv").exists()

    proc = run_cli("write-state", "--run-id", "demo", "--runs-dir", str(tmp_path / "runs"), "--status", "running")
    assert proc.returncode == 0, proc.stderr
    proc = run_cli("status", "--run-id", "demo", "--runs-dir", str(tmp_path / "runs"))
    assert proc.returncode == 0
    assert '"status": "running"' in proc.stdout


def test_enrichment_maps_ensembl_inputs_to_symbol_gmt(tmp_path):
    significant = tmp_path / "significant_genes.csv"
    significant.write_text("gene_id,log2FoldChange\nENSG001.5,2\nENSG003,-2\n")
    ranked = tmp_path / "ranked_genes.csv"
    ranked.write_text("gene_id,stat\nENSG001.5,5\nENSG002,4\nENSG003,-4\n")
    universe = tmp_path / "universe.csv"
    universe.write_text("gene_id\nENSG001.5\nENSG002\nENSG003\nENSG004\n")
    mapping = tmp_path / "mapping.csv"
    mapping.write_text(
        "gene_id,gene_symbol,ensembl_id,entrez_id\n"
        "ENSG001,GENEA,ENSG001,1\n"
        "ENSG002,GENEB,ENSG002,2\n"
        "ENSG003,GENEC,ENSG003,3\n"
    )
    gmt = tmp_path / "symbol_sets.gmt"
    gmt.write_text("SYMBOL_UP\tsymbol demo\tGENEA\tGENEB\nSYMBOL_DOWN\tsymbol demo\tGENEC\n")
    outdir = tmp_path / "mapped_go"
    proc = run_cli(
        "go-enrichment",
        "--significant-genes",
        str(significant),
        "--ranked-genes",
        str(ranked),
        "--universe",
        str(universe),
        "--gmt",
        str(gmt),
        "--gene-mapping",
        str(mapping),
        "--input-id-type",
        "ensembl_id",
        "--gene-set-id-type",
        "gene_symbol",
        "--outdir",
        str(outdir),
    )
    assert proc.returncode == 0, proc.stderr
    assert "SYMBOL_UP" in (outdir / "go_overrepresentation_up.csv").read_text()
    assert "SYMBOL_DOWN" in (outdir / "go_overrepresentation_down.csv").read_text()
    assert "padj" in (outdir / "go_overrepresentation_up.csv").read_text().splitlines()[0]
