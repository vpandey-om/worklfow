from pathlib import Path

import yaml

from .models import RunRequest, StepSpec


def load_steps(registry_path: str | Path) -> dict[str, StepSpec]:
    path = Path(registry_path)
    data = yaml.safe_load(path.read_text())
    steps = {}
    for step_id, raw in data.get("steps", {}).items():
        step = StepSpec.model_validate(raw)
        if step.step_id != step_id:
            raise ValueError(f"Registry key {step_id} does not match step_id {step.step_id}")
        steps[step_id] = step
    return steps


def load_run_request(request_path: str | Path) -> RunRequest:
    data = yaml.safe_load(Path(request_path).read_text())
    return RunRequest.model_validate(data)


def project_root_from_builder() -> Path:
    return Path(__file__).resolve().parents[1]

