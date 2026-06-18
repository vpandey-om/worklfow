from pathlib import Path

from builder.planner import WorkflowPlanner
from builder.registry import load_steps


ROOT = Path(__file__).resolve().parents[1]


def planner():
    return WorkflowPlanner(load_steps(ROOT / "registry" / "steps.yaml"))


def test_planner_recommends_trimming_after_raw_qc():
    suggestions = planner().suggest_next_steps(["rnaseq_03_raw_read_qc"], "rnaseq")
    assert suggestions["recommended"][0].step_id == "rnaseq_04_adapter_quality_trimming"


def test_trimming_can_run_individually():
    result = planner().validate_selected_steps(["rnaseq_04_adapter_quality_trimming"], "rnaseq")
    assert result.valid


def test_unknown_step_is_invalid():
    result = planner().validate_selected_steps(["not_a_step"], "rnaseq")
    assert not result.valid


def test_resolves_step_namespaced_trimming_params():
    resolved = planner().resolve_effective_params(
        ["rnaseq_04_adapter_quality_trimming"],
        {
            "rnaseq_04_adapter_quality_trimming": {
                "quality_threshold": 25,
                "minimum_read_length": 30,
                "trim_poly_g": "true",
            }
        },
    )
    assert resolved["quality_threshold"] == 25
    assert resolved["minimum_read_length"] == 30
    assert resolved["trim_poly_g"] == "true"


def test_cutadapt_requires_adapters():
    resolved = planner().resolve_effective_params(
        ["rnaseq_04_adapter_quality_trimming"],
        {"trimming_tool": "cutadapt"},
    )
    result = planner().validate_params(["rnaseq_04_adapter_quality_trimming"], resolved)
    assert not result.valid
    assert any(error["param"] == "adapter_sequence_r1" for error in result.errors)
