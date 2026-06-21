from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from werkzeug.utils import secure_filename


ALLOWED_FASTQ_EXTENSIONS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
ALLOWED_METADATA_EXTENSIONS = (".csv", ".tsv", ".xlsx")
ALLOWED_REFERENCE_EXTENSIONS = (
    ".fa",
    ".fasta",
    ".fa.gz",
    ".fasta.gz",
    ".gtf",
    ".gtf.gz",
    ".gff",
    ".gff3",
    ".gff.gz",
    ".gff3.gz",
    ".tsv",
    ".csv",
    ".json",
    ".zip",
    ".tar",
    ".tar.gz",
    ".tgz",
)
ALLOWED_VENDOR_EXTENSIONS = (".raw", ".mzml", ".mzml.gz", ".zip", ".tar", ".tar.gz", ".tgz")
ALLOWED_UPLOAD_TYPES = {"fastq", "metadata", "reference", "vendor", "other"}


def safe_session_id(session_id: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", session_id or "")
    if not cleaned or cleaned in {".", ".."}:
        raise ValueError("Invalid session id")
    return cleaned[:120]


def safe_slug(value: str, default: str = "unknown") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", value or "")
    if not cleaned or cleaned in {".", ".."}:
        cleaned = default
    return cleaned[:80]


def safe_filename(filename: str) -> str:
    name = secure_filename(filename or "")
    if not name:
        raise ValueError("Invalid filename")
    return name


def has_allowed_extension(filename: str, upload_type: str) -> bool:
    if upload_type == "other":
        return True
    lower = filename.lower()
    if upload_type == "fastq":
        allowed = ALLOWED_FASTQ_EXTENSIONS
    elif upload_type == "metadata":
        allowed = ALLOWED_METADATA_EXTENSIONS
    elif upload_type == "reference":
        allowed = ALLOWED_REFERENCE_EXTENSIONS
    elif upload_type == "vendor":
        allowed = ALLOWED_VENDOR_EXTENSIONS
    else:
        return False
    return lower.endswith(allowed)


class UploadService:
    def __init__(self, upload_root: Path, max_upload_bytes: int):
        self.upload_root = Path(upload_root)
        self.max_upload_bytes = max_upload_bytes
        self.upload_root.mkdir(parents=True, exist_ok=True)

    def session_dir(
        self,
        session_id: str,
        tester_id: str | None = None,
        omics_type: str | None = None,
        project_id: str | None = None,
    ) -> Path:
        sid = safe_session_id(session_id)
        if tester_id and omics_type:
            project = safe_slug(project_id or "demo_project")
            path = self.upload_root / safe_slug(omics_type) / safe_slug(tester_id) / project / sid
        else:
            path = self.upload_root / sid
        path.mkdir(parents=True, exist_ok=True)
        return path

    def upload_type_dir(
        self,
        session_id: str,
        upload_type: str,
        tester_id: str | None = None,
        omics_type: str | None = None,
        project_id: str | None = None,
    ) -> Path:
        if upload_type not in ALLOWED_UPLOAD_TYPES:
            raise ValueError("upload_type must be fastq, metadata, reference, vendor, or other")
        path = self.session_dir(session_id, tester_id, omics_type, project_id) / upload_type
        path.mkdir(parents=True, exist_ok=True)
        return path

    def chunk_path(
        self,
        session_id: str,
        upload_type: str,
        identifier: str,
        chunk_number: int,
        tester_id: str | None = None,
        omics_type: str | None = None,
        project_id: str | None = None,
    ) -> Path:
        ident = safe_filename(identifier)
        chunk_dir = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / ".chunks" / ident
        chunk_dir.mkdir(parents=True, exist_ok=True)
        return chunk_dir / f"{chunk_number:08d}.part"

    def save_chunk(self, session_id: str, upload_type: str, args, file_storage, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> dict:
        filename = safe_filename(args.get("resumableFilename", ""))
        if not has_allowed_extension(filename, upload_type):
            raise ValueError(f"File extension is not allowed for {upload_type}: {filename}")

        total_size = int(args.get("resumableTotalSize", "0"))
        if self.max_upload_bytes and total_size > self.max_upload_bytes:
            raise ValueError("File exceeds configured upload limit")

        identifier = args.get("resumableIdentifier", filename)
        chunk_number = int(args.get("resumableChunkNumber", "1"))
        total_chunks = int(args.get("resumableTotalChunks", "1"))
        part = self.chunk_path(session_id, upload_type, identifier, chunk_number, tester_id, omics_type, project_id)
        file_storage.save(part)

        return {
            "filename": filename,
            "chunk_number": chunk_number,
            "total_chunks": total_chunks,
            "complete": self.is_complete(session_id, upload_type, identifier, total_chunks, tester_id, omics_type, project_id),
        }

    def save_file(self, session_id: str, upload_type: str, file_storage, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> dict:
        filename = safe_filename(file_storage.filename or "")
        if not has_allowed_extension(filename, upload_type):
            raise ValueError(f"File extension is not allowed for {upload_type}: {filename}")

        target = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / filename
        file_storage.save(target)
        if self.max_upload_bytes and target.stat().st_size > self.max_upload_bytes:
            target.unlink(missing_ok=True)
            raise ValueError("File exceeds configured upload limit")

        self.write_manifest(session_id, tester_id, omics_type, project_id)
        return {"filename": filename, "path": str(target), "size": target.stat().st_size}

    def save_bytes(self, session_id: str, upload_type: str, filename: str, content: bytes, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> dict:
        safe_name = safe_filename(filename)
        if not has_allowed_extension(safe_name, upload_type):
            raise ValueError(f"File extension is not allowed for {upload_type}: {safe_name}")
        if self.max_upload_bytes and len(content) > self.max_upload_bytes:
            raise ValueError("File exceeds configured upload limit")

        target = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / safe_name
        target.write_bytes(content)
        self.write_manifest(session_id, tester_id, omics_type, project_id)
        return {
            "filename": safe_name,
            "path": str(target),
            "size": target.stat().st_size,
            "uploaded_at": path_mtime_iso(target),
            "category": upload_type,
            "project_id": safe_slug(project_id or "demo_project"),
        }

    def is_complete(self, session_id: str, upload_type: str, identifier: str, total_chunks: int, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> bool:
        chunk_dir = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / ".chunks" / safe_filename(identifier)
        return all((chunk_dir / f"{i:08d}.part").exists() for i in range(1, total_chunks + 1))

    def assemble_file(self, session_id: str, upload_type: str, args, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> dict:
        filename = safe_filename(args.get("resumableFilename", ""))
        if not has_allowed_extension(filename, upload_type):
            raise ValueError(f"File extension is not allowed for {upload_type}: {filename}")

        identifier = args.get("resumableIdentifier", filename)
        total_chunks = int(args.get("resumableTotalChunks", "1"))
        target = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / filename
        chunk_dir = self.upload_type_dir(session_id, upload_type, tester_id, omics_type, project_id) / ".chunks" / safe_filename(identifier)
        with target.open("wb") as out:
            for i in range(1, total_chunks + 1):
                part = chunk_dir / f"{i:08d}.part"
                if not part.exists():
                    raise FileNotFoundError(f"Missing upload chunk {i}")
                with part.open("rb") as handle:
                    shutil.copyfileobj(handle, out)
        shutil.rmtree(chunk_dir, ignore_errors=True)
        self.write_manifest(session_id, tester_id, omics_type, project_id)
        return {"filename": filename, "path": str(target), "size": target.stat().st_size}

    def list_uploads(self, session_id: str, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> dict:
        session = self.session_dir(session_id, tester_id, omics_type, project_id)
        result = {
            "session_id": safe_session_id(session_id),
            "tester_id": safe_slug(tester_id) if tester_id else None,
            "omics_type": safe_slug(omics_type) if omics_type else None,
            "project_id": safe_slug(project_id or "demo_project") if tester_id and omics_type else None,
            "fastq": [],
            "metadata": [],
            "reference": [],
            "vendor": [],
            "other": [],
        }
        for upload_type in ("fastq", "metadata", "reference", "vendor", "other"):
            folder = session / upload_type
            if not folder.exists():
                continue
            for path in sorted(folder.iterdir()):
                if path.is_file():
                    result[upload_type].append({
                        "name": path.name,
                        "path": str(path),
                        "size": path.stat().st_size,
                        "uploaded_at": path_mtime_iso(path),
                        "category": upload_type,
                    })
        return result

    def list_workspace_uploads(self, tester_id: str, omics_type: str, project_id: str | None = None) -> dict:
        tester = safe_slug(tester_id)
        omics = safe_slug(omics_type)
        project = safe_slug(project_id or "demo_project")
        workspace = self.upload_root / omics / tester / project
        result = {
            "tester_id": tester,
            "omics_type": omics,
            "project_id": project,
            "fastq": [],
            "metadata": [],
            "reference": [],
            "vendor": [],
            "other": [],
        }
        if not workspace.exists():
            return result
        for session_dir in sorted(path for path in workspace.iterdir() if path.is_dir()):
            for upload_type in ("fastq", "metadata", "reference", "vendor", "other"):
                folder = session_dir / upload_type
                if not folder.exists():
                    continue
                for path in sorted(folder.iterdir()):
                    if path.is_file():
                        result[upload_type].append({
                            "name": path.name,
                            "path": str(path),
                            "size": path.stat().st_size,
                            "uploaded_at": path_mtime_iso(path),
                            "category": upload_type,
                            "session_id": session_dir.name,
                            "project_id": project,
                        })
        return result

    def write_manifest(self, session_id: str, tester_id: str | None = None, omics_type: str | None = None, project_id: str | None = None) -> Path:
        path = self.session_dir(session_id, tester_id, omics_type, project_id) / "uploads_manifest.json"
        path.write_text(json.dumps(self.list_uploads(session_id, tester_id, omics_type, project_id), indent=2), encoding="utf-8")
        return path


def path_mtime_iso(path: Path) -> str:
    from datetime import datetime, timezone

    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
