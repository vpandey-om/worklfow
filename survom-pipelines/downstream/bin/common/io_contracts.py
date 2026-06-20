from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Callable

from .checksums import write_checksums
from .metadata import write_metadata


class StudentError(RuntimeError):
    """Error message intended to be shown directly to students."""


def repo_root_from_file(path: str | Path) -> Path:
    current = Path(path).resolve()
    for parent in [current, *current.parents]:
        if (parent / "survom-pipelines").exists() and (parent / "demo_dash_app").exists():
            return parent
    return Path(__file__).resolve().parents[5]


def write_command(outdir: Path, argv: list[str]) -> None:
    (Path(outdir) / "command.txt").write_text(" ".join(argv) + "\n", encoding="utf-8")


def run_atomic_step(
    *,
    step_name: str,
    outdir: Path,
    argv: list[str],
    params: dict,
    inputs: dict[str, str],
    work: Callable[[Path], None],
    repo_root: Path,
) -> int:
    outdir = Path(outdir).resolve()
    outdir.parent.mkdir(parents=True, exist_ok=True)
    in_place = outdir == Path.cwd().resolve()
    tmp_parent = outdir if in_place else outdir.parent
    tmp = Path(tempfile.mkdtemp(prefix=f".{outdir.name or 'out'}.", dir=tmp_parent))
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    exit_code = 0
    try:
        write_command(tmp, argv)
        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            work(tmp)
        write_metadata(tmp, step_name, params, inputs, repo_root)
        write_checksums(tmp)
        if in_place:
            for path in tmp.iterdir():
                destination = outdir / path.name
                if destination.exists():
                    if destination.is_dir():
                        shutil.rmtree(destination)
                    else:
                        destination.unlink()
                os.replace(path, destination)
            shutil.rmtree(tmp, ignore_errors=True)
        elif outdir.exists():
            shutil.rmtree(outdir)
            os.replace(tmp, outdir)
        else:
            os.replace(tmp, outdir)
    except Exception as exc:
        exit_code = 1
        print(str(exc), file=stderr_buffer)
        failed = outdir.parent / f"{outdir.name}.failed"
        if failed.exists():
            shutil.rmtree(failed)
        os.replace(tmp, failed)
        outdir = failed
    finally:
        (outdir / "stdout.log").write_text(stdout_buffer.getvalue(), encoding="utf-8")
        (outdir / "stderr.log").write_text(stderr_buffer.getvalue(), encoding="utf-8")
    return exit_code


def require_file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists() or not resolved.is_file():
        raise StudentError(f"{label} does not exist or is not a file: {path}")
    return resolved
