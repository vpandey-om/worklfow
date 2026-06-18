import argparse
import json
from pathlib import Path

from .compiler import NextflowCompiler
from .outputs import collect_run_outputs
from .planner import WorkflowPlanner
from .registry import load_run_request, load_steps, project_root_from_builder
from .runner import run_job, run_nextflow


def load_context(root: Path):
    steps = load_steps(root / "registry" / "steps.yaml")
    compiler = NextflowCompiler(steps, root / "workflows" / "templates")
    return steps, compiler


def compile_request(request_path: str, output_dir: str | None = None) -> Path:
    root = project_root_from_builder()
    _, compiler = load_context(root)
    request = load_run_request(request_path)
    out = output_dir or root / "workflows" / "generated" / request.run_id
    return compiler.compile(request, out)


def create_app():
    try:
        from fastapi import FastAPI
    except ImportError as exc:
        raise RuntimeError("Install fastapi and uvicorn to serve the SurvOm API.") from exc

    root = project_root_from_builder()
    steps, compiler = load_context(root)
    planner = WorkflowPlanner(steps)
    app = FastAPI(title="SurvOm Workflow Builder API")

    @app.get("/api/omics")
    def get_omics():
        return sorted({omics for step in steps.values() for omics in step.omics})

    @app.get("/api/steps")
    def get_steps(omics: str | None = None):
        selected = steps.values() if omics is None else planner.get_steps_for_omics(omics)
        return [step.model_dump() for step in selected]

    @app.get("/api/steps/{step_id}")
    def get_step(step_id: str):
        return steps[step_id].model_dump()

    @app.post("/api/workflows/plan")
    def plan_workflow(payload: dict):
        omics = payload["omics"]
        selected_steps = payload.get("selected_steps", [])
        validation = planner.validate_selected_steps(selected_steps, omics)
        suggestions = planner.suggest_next_steps(selected_steps, omics)
        return {
            "valid": validation.valid,
            "selected_steps": selected_steps,
            "errors": validation.errors,
            "warnings": validation.warnings,
            "suggested_fixes": validation.suggested_fixes,
            "recommended_next_steps": [s.model_dump() for s in suggestions["recommended"]],
            "allowed_next_steps": [s.step_id for s in suggestions["allowed"]],
        }

    @app.post("/api/workflows/compile")
    def compile_workflow(payload: dict):
        from .models import RunRequest

        request = RunRequest.model_validate(payload)
        out = compiler.compile(request, root / "workflows" / "generated" / request.run_id)
        return {"compiled_dir": str(out)}

    @app.post("/api/workflows/run")
    def run_workflow(payload: dict):
        proc = run_nextflow(payload["compiled_dir"], payload.get("profile", "local,docker"))
        return {"returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}

    @app.get("/api/runs/{run_id}")
    def get_run(run_id: str):
        manifest = root / "workflows" / "generated" / run_id / "workflow_manifest.json"
        if manifest.exists():
            return json.loads(manifest.read_text())
        return {"run_id": run_id, "status": "not_found"}

    @app.get("/api/runs/{run_id}/outputs")
    def get_run_outputs(run_id: str):
        manifest = root / "workflows" / "generated" / run_id / "workflow_manifest.json"
        if not manifest.exists():
            return {"run_id": run_id, "outputs": []}
        data = json.loads(manifest.read_text())
        outdir = Path(data["params"]["outdir"])
        return {"run_id": run_id, "outputs": [str(p) for p in outdir.rglob("*") if p.is_file()] if outdir.exists() else []}

    return app


def main():
    parser = argparse.ArgumentParser(prog="survom", description="SurvOm workflow builder")
    sub = parser.add_subparsers(dest="command", required=True)

    compile_cmd = sub.add_parser("compile")
    compile_cmd.add_argument("--request", required=True)
    compile_cmd.add_argument("--output")

    build_run = sub.add_parser("build-run")
    build_run.add_argument("--request", required=True)
    build_run.add_argument("--output")

    run_cmd = sub.add_parser("run")
    run_cmd.add_argument("--compiled-dir", required=True)
    run_cmd.add_argument("--profile", default="local,docker")

    collect_cmd = sub.add_parser("collect-outputs")
    collect_cmd.add_argument("--compiled-dir", required=True)
    collect_cmd.add_argument("--status", default="unknown")

    job_cmd = sub.add_parser("run-job")
    job_cmd.add_argument("--job", required=True)

    args = parser.parse_args()
    if args.command in {"compile", "build-run"}:
        out = compile_request(args.request, args.output)
        print(json.dumps({"compiled_dir": str(out)}, indent=2))
    elif args.command == "run":
        proc = run_nextflow(args.compiled_dir, args.profile)
        print(proc.stdout)
        if proc.stderr:
            print(proc.stderr)
        raise SystemExit(proc.returncode)
    elif args.command == "collect-outputs":
        out = collect_run_outputs(args.compiled_dir, args.status)
        print(json.dumps({"execution_record": str(out)}, indent=2))
    elif args.command == "run-job":
        out = run_job(args.job)
        print(json.dumps({"job_record": str(out)}, indent=2))


if __name__ == "__main__":
    main()
