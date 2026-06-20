from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from .io_contracts import StudentError


FEATURECOUNTS_TECHNICAL_COLUMNS = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}


def read_table(path: Path) -> pd.DataFrame:
    first_line = Path(path).read_text(errors="replace", encoding="utf-8").splitlines()[0]
    sep = "\t" if "\t" in first_line or path.suffix.lower() == ".tsv" or path.name.endswith(".tsv") else ","
    return pd.read_csv(path, sep=sep)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def normalize_gene_column(df: pd.DataFrame) -> pd.DataFrame:
    if "gene_id" in df.columns:
        return df
    if "Geneid" in df.columns:
        return df.rename(columns={"Geneid": "gene_id"})
    first = df.columns[0]
    return df.rename(columns={first: "gene_id"})


def validate_counts_frame(df: pd.DataFrame) -> list[str]:
    errors = []
    if "gene_id" not in df.columns:
        errors.append("Counts matrix must contain a gene_id column.")
        return errors
    if df["gene_id"].isna().any() or (df["gene_id"].astype(str).str.strip() == "").any():
        errors.append("Counts matrix contains blank gene_id values.")
    duplicates = df["gene_id"][df["gene_id"].duplicated()].astype(str).unique().tolist()
    if duplicates:
        errors.append(f"Duplicate gene_id values found: {', '.join(duplicates[:10])}")
    sample_columns = [col for col in df.columns if col != "gene_id"]
    if not sample_columns:
        errors.append("Counts matrix must contain at least one sample column.")
    if len(sample_columns) != len(set(sample_columns)):
        errors.append("Sample columns in counts matrix must be unique.")
    for col in sample_columns:
        values = pd.to_numeric(df[col], errors="coerce")
        if values.isna().any():
            errors.append(f"Counts column '{col}' contains missing or non-numeric values.")
        elif not np.isfinite(values).all():
            errors.append(f"Counts column '{col}' contains non-finite values.")
        elif (values < 0).any():
            errors.append(f"Counts column '{col}' contains negative counts.")
        elif not np.equal(values, np.floor(values)).all():
            errors.append(f"Counts column '{col}' contains non-integer counts.")
    return errors


def clean_counts(df: pd.DataFrame) -> pd.DataFrame:
    df = normalize_gene_column(df)
    keep = ["gene_id"] + [col for col in df.columns if col not in FEATURECOUNTS_TECHNICAL_COLUMNS and col != "gene_id"]
    cleaned = df.loc[:, keep].copy()
    for col in cleaned.columns:
        if col != "gene_id":
            cleaned[col] = pd.to_numeric(cleaned[col], errors="raise").astype(int)
    return cleaned


def validate_metadata(metadata: pd.DataFrame, sample_ids: list[str], group_column: str) -> list[str]:
    errors = []
    if "sample_id" not in metadata.columns:
        errors.append("Metadata must contain a sample_id column.")
        return errors
    if metadata["sample_id"].duplicated().any():
        errors.append("Metadata sample_id values must be unique.")
    if group_column not in metadata.columns:
        errors.append(f"Metadata must contain grouping column '{group_column}'.")
    metadata_samples = set(metadata["sample_id"].astype(str))
    count_samples = set(sample_ids)
    missing_meta = sorted(count_samples - metadata_samples)
    extra_meta = sorted(metadata_samples - count_samples)
    if missing_meta:
        errors.append(f"Metadata is missing samples present in counts: {', '.join(missing_meta[:10])}")
    if extra_meta:
        errors.append(f"Metadata contains samples absent from counts: {', '.join(extra_meta[:10])}")
    return errors


def read_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise StudentError(f"GMT line {line_no} must contain name, description, and at least one gene.")
            gene_sets[parts[0]] = set(gene for gene in parts[2:] if gene)
    return gene_sets
