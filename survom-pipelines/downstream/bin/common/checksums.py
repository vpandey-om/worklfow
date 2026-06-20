from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_checksums(outdir: Path) -> Path:
    outdir = Path(outdir)
    target = outdir / "checksums.sha256"
    lines = []
    for path in sorted(outdir.rglob("*")):
        if not path.is_file() or path.name == "checksums.sha256":
            continue
        lines.append(f"{sha256_file(path)}  {path.relative_to(outdir)}")
    target.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return target
