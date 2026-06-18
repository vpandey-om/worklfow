from typing import Any, Literal

from pydantic import BaseModel, Field

Exposure = Literal["user", "expert", "internal"]


class ParameterSpec(BaseModel):
    type: str
    default: Any = None
    choices: list[Any] | None = None
    min: float | None = None
    max: float | None = None
    exposure: Exposure
    ui_label: str | None = None


class OutputSpec(BaseModel):
    name: str
    type: str
    pattern: str | None = None
    ui_visible: bool = True


class InputSpec(BaseModel):
    name: str
    type: str
    formats: list[str] = Field(default_factory=list)
    required: bool = True
    required_when: str | None = None


class StepSpec(BaseModel):
    step_id: str
    step_name: str
    category: str
    omics: list[str]
    purpose: str
    module_path: str
    process_name: str
    default_tool: str
    fallback_tool: str | None = None
    can_run_individually: bool = True
    can_be_skipped: bool = False
    skip_conditions: list[str] = Field(default_factory=list)
    depends_on: list[str] = Field(default_factory=list)
    optional_depends_on: list[str] = Field(default_factory=list)
    compatible_previous_outputs: list[str] = Field(default_factory=list)
    required_inputs: list[InputSpec] = Field(default_factory=list)
    outputs: list[OutputSpec] = Field(default_factory=list)
    parameters: dict[str, ParameterSpec] = Field(default_factory=dict)
    next_steps: dict[str, list[str]] = Field(default_factory=dict)
    resources: dict[str, Any] = Field(default_factory=dict)
    containers: dict[str, str] = Field(default_factory=dict)
    ui: dict[str, Any] = Field(default_factory=dict)
    chat_guidance: dict[str, str] = Field(default_factory=dict)


class StepSuggestion(BaseModel):
    step_id: str
    label: str
    reason: str
    recommended: bool = False


class RunRequest(BaseModel):
    run_id: str
    omics: str
    execution_profile: str
    input: str
    outdir: str
    selected_steps: list[str]
    params: dict[str, Any] = Field(default_factory=dict)
    user_id: str | None = None
    dataset_id: str | None = None
    run_label: str | None = None
    run_description: str | None = None


class ValidationResult(BaseModel):
    valid: bool
    errors: list[dict[str, Any]] = Field(default_factory=list)
    warnings: list[dict[str, Any]] = Field(default_factory=list)
    suggested_fixes: list[str] = Field(default_factory=list)


class JobRequest(BaseModel):
    job_id: str
    run_id: str
    compiled_dir: str
    command: str
    results_dir: str
    execution_profile: str = "local"
    logs_dir: str | None = None
