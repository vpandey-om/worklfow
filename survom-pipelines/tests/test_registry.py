from pathlib import Path

from builder.registry import load_steps


ROOT = Path(__file__).resolve().parents[1]


def test_registry_loads_and_step_ids_are_unique():
    steps = load_steps(ROOT / "registry" / "steps.yaml")
    assert len(steps) == len({s.step_id for s in steps.values()})
    assert "rnaseq_03_raw_read_qc" in steps
    assert "rnaseq_04_adapter_quality_trimming" in steps


def test_common_steps_are_reusable_across_omics():
    steps = load_steps(ROOT / "registry" / "steps.yaml")
    assert len(steps["rnaseq_03_raw_read_qc"].omics) > 1
    assert len(steps["rnaseq_04_adapter_quality_trimming"].omics) > 1


def test_parameter_exposure_levels_exist():
    steps = load_steps(ROOT / "registry" / "steps.yaml")
    exposures = {p.exposure for s in steps.values() for p in s.parameters.values()}
    assert {"user", "expert", "internal"}.issubset(exposures)
