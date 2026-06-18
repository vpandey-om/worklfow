from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class JobStore:
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.init_db()

    def connect(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def init_db(self):
        with self.connect() as conn:
            # This is a demo/tester workspace model, not production authentication.
            # For production, replace this with proper login/auth and user authorization.
            conn.execute(
                """
                create table if not exists jobs (
                    job_id text primary key,
                    run_id text not null,
                    session_id text,
                    tester_id text,
                    tester_label text,
                    omics_type text,
                    project_id text,
                    project_label text,
                    workflow_id text,
                    status text not null,
                    pid integer,
                    run_dir text,
                    compiled_dir text,
                    results_dir text,
                    command text,
                    stdout_path text,
                    stderr_path text,
                    created_at text,
                    started_at text,
                    finished_at text,
                    exit_code integer,
                    metadata_json text
                )
                """
            )
            existing = {row["name"] for row in conn.execute("pragma table_info(jobs)").fetchall()}
            for column in ("tester_id", "tester_label", "omics_type", "project_id", "project_label"):
                if column not in existing:
                    conn.execute(f"alter table jobs add column {column} text")
            conn.execute(
                """
                create table if not exists projects (
                    project_id text not null,
                    project_label text not null,
                    tester_id text not null,
                    omics_type text not null,
                    created_at text,
                    updated_at text,
                    primary key (project_id, tester_id, omics_type)
                )
                """
            )
            conn.execute(
                """
                create table if not exists upload_records (
                    id integer primary key autoincrement,
                    project_id text,
                    tester_id text,
                    omics_type text,
                    session_id text,
                    category text,
                    filename text,
                    stored_path text unique,
                    size_bytes integer,
                    created_at text
                )
                """
            )

    def upsert_job(self, record: dict[str, Any]):
        now = datetime.now(timezone.utc).isoformat()
        record.setdefault("created_at", now)
        record.setdefault("metadata_json", "{}")
        columns = [
            "job_id", "run_id", "session_id", "tester_id", "tester_label", "omics_type",
            "project_id", "project_label", "workflow_id", "status", "pid", "run_dir",
            "compiled_dir", "results_dir", "command", "stdout_path", "stderr_path",
            "created_at", "started_at", "finished_at", "exit_code", "metadata_json",
        ]
        values = [record.get(col) for col in columns]
        placeholders = ", ".join(["?"] * len(columns))
        update = ", ".join([f"{col}=excluded.{col}" for col in columns if col != "job_id"])
        with self.connect() as conn:
            conn.execute(
                f"""
                insert into jobs ({", ".join(columns)})
                values ({placeholders})
                on conflict(job_id) do update set {update}
                """,
                values,
            )

    def update_job(self, job_id: str, **updates):
        if not updates:
            return
        assignments = ", ".join([f"{key}=?" for key in updates])
        values = list(updates.values()) + [job_id]
        with self.connect() as conn:
            conn.execute(f"update jobs set {assignments} where job_id=?", values)

    def get_job(self, job_id: str) -> dict[str, Any] | None:
        with self.connect() as conn:
            row = conn.execute("select * from jobs where job_id=?", (job_id,)).fetchone()
        return dict(row) if row else None

    def list_jobs(self, limit: int = 50) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute("select * from jobs order by created_at desc limit ?", (limit,)).fetchall()
        return [dict(row) for row in rows]

    def list_workspace_jobs(
        self,
        session_id: str,
        tester_id: str,
        omics_type: str,
        project_id: str | None = None,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        query = """
            select * from jobs
            where session_id=? and tester_id=? and omics_type=?
        """
        values: list[Any] = [session_id, tester_id, omics_type]
        if project_id:
            query += " and project_id=?"
            values.append(project_id)
        query += " order by created_at desc limit ?"
        values.append(limit)
        with self.connect() as conn:
            rows = conn.execute(query, values).fetchall()
        return [dict(row) for row in rows]

    def upsert_project(self, tester_id: str, omics_type: str, project_id: str, project_label: str):
        now = datetime.now(timezone.utc).isoformat()
        with self.connect() as conn:
            conn.execute(
                """
                insert into projects (project_id, project_label, tester_id, omics_type, created_at, updated_at)
                values (?, ?, ?, ?, ?, ?)
                on conflict(project_id, tester_id, omics_type) do update set
                    project_label=excluded.project_label,
                    updated_at=excluded.updated_at
                """,
                (project_id, project_label, tester_id, omics_type, now, now),
            )

    def list_projects(self, tester_id: str, omics_type: str) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                select * from projects
                where tester_id=? and omics_type=?
                order by updated_at desc, project_label
                """,
                (tester_id, omics_type),
            ).fetchall()
        return [dict(row) for row in rows]

    def upsert_upload_record(self, record: dict[str, Any]):
        record.setdefault("created_at", datetime.now(timezone.utc).isoformat())
        columns = [
            "project_id", "tester_id", "omics_type", "session_id", "category",
            "filename", "stored_path", "size_bytes", "created_at",
        ]
        values = [record.get(col) for col in columns]
        with self.connect() as conn:
            conn.execute(
                f"""
                insert into upload_records ({", ".join(columns)})
                values ({", ".join(["?"] * len(columns))})
                on conflict(stored_path) do update set
                    project_id=excluded.project_id,
                    tester_id=excluded.tester_id,
                    omics_type=excluded.omics_type,
                    session_id=excluded.session_id,
                    category=excluded.category,
                    filename=excluded.filename,
                    size_bytes=excluded.size_bytes
                """,
                values,
            )

    def list_upload_records(self, tester_id: str, omics_type: str, project_id: str) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                select * from upload_records
                where tester_id=? and omics_type=? and project_id=?
                order by created_at desc, filename
                """,
                (tester_id, omics_type, project_id),
            ).fetchall()
        return [dict(row) for row in rows]

    def refresh_from_job_record(self, job_id: str, compiled_dir: Path):
        record_path = Path(compiled_dir) / "job_record.json"
        if not record_path.exists():
            return
        record = json.loads(record_path.read_text())
        self.update_job(
            job_id,
            status=record.get("status"),
            started_at=record.get("started_at"),
            finished_at=record.get("finished_at"),
            exit_code=record.get("exit_code"),
            stdout_path=record.get("stdout_path"),
            stderr_path=record.get("stderr_path"),
        )
