from pathlib import Path

from builder.compiler import NextflowCompiler
from builder.registry import load_run_request, load_steps


ROOT = Path(__file__).resolve().parents[1]


def test_generated_config_has_local_docker_singularity_and_aws_profiles(tmp_path):
    steps = load_steps(ROOT / "registry" / "steps.yaml")
    req = load_run_request(ROOT / "examples" / "rnaseq_qc_trim.run.yaml")
    out = NextflowCompiler(steps, ROOT / "workflows" / "templates").compile(req, tmp_path / req.run_id)
    config = (out / "nextflow.config").read_text()
    assert "local {" in config
    assert "docker {" in config
    assert "singularity {" in config
    assert "aws {" in config
    assert "aws_queue" in config
    assert "aws_region" in config
    assert "aws_workdir" in config
