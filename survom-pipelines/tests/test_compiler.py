from pathlib import Path

from builder.compiler import NextflowCompiler
from builder.registry import load_run_request, load_steps


ROOT = Path(__file__).resolve().parents[1]


def compiler():
    return NextflowCompiler(load_steps(ROOT / "registry" / "steps.yaml"), ROOT / "workflows" / "templates")


def test_compile_qc_only(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_qc_only.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "include { FASTQC }" in main
    assert "FASTQC(ch_reads)" in main
    assert "FASTP(ch_reads)" not in main


def test_compile_qc_trim(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_qc_trim.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    params = (out / "params.yaml").read_text()
    assert "include { FASTQC }" in main
    assert "include { FASTP }" in main
    assert "MAKE_PREPROCESS_MANIFEST(FASTP.out.reads)" in main
    assert "MERGE_PREPROCESS_MANIFEST(MAKE_PREPROCESS_MANIFEST.out.fragments.collect())" in main
    assert "quality_threshold: 20" in params


def test_compile_trim_only(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_trim_only.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "include { FASTP }" in main
    assert "include { FASTQC }" not in main
    assert "FASTP(ch_reads)" in main


def test_compile_step_namespaced_custom_trimming_params(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_trim_only.run.yaml")
    req.params = {
        "rnaseq_04_adapter_quality_trimming": {
            "trimming_tool": "fastp",
            "quality_threshold": 25,
            "minimum_read_length": 30,
            "trim_poly_g": "true",
            "trim_front_r1": 5,
        }
    }
    out = compiler().compile(req, tmp_path / "custom_trim")
    params = (out / "params.yaml").read_text()
    assert "quality_threshold: 25" in params
    assert "minimum_read_length: 30" in params
    assert 'trim_poly_g: "true"' in params
    assert "trim_front_r1: 5" in params
