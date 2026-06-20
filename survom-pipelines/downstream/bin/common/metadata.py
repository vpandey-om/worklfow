from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0"


def _run_text(command: list[str], cwd: Path | None = None) -> str | None:
    try:
        proc = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def git_info(root: Path) -> dict[str, Any]:
    return {
        "commit": _run_text(["git", "rev-parse", "HEAD"], root),
        "branch": _run_text(["git", "branch", "--show-current"], root),
        "dirty": bool(_run_text(["git", "status", "--porcelain"], root)),
    }


def base_metadata(step_name: str, params: dict[str, Any], inputs: dict[str, str], repo_root: Path) -> dict[str, Any]:
    version_path = repo_root / "survom-pipelines" / "downstream" / "VERSION"
    return {
        "schema_version": SCHEMA_VERSION,
        "step_name": step_name,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "workflow_version": version_path.read_text(encoding="utf-8").strip() if version_path.exists() else "unknown",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "git": git_info(repo_root),
        "inputs": inputs,
        "parameters": params,
    }


def write_metadata(outdir: Path, step_name: str, params: dict[str, Any], inputs: dict[str, str], repo_root: Path) -> Path:
    target = Path(outdir) / "metadata.json"
    target.write_text(json.dumps(base_metadata(step_name, params, inputs, repo_root), indent=2) + "\n", encoding="utf-8")
    return target
