#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from common.io_contracts import StudentError, repo_root_from_file, require_file, run_atomic_step
from common.validation import clean_counts, read_gmt, read_table, validate_counts_frame, validate_metadata, write_csv


REPO_ROOT = repo_root_from_file(__file__)


def _sample_name_from_count_column(column: str) -> str:
    name = Path(str(column)).name
    for suffix in (".bam", ".sam", ".Aligned.sortedByCoord.out.bam"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name.replace(".featureCounts", "")


def cmd_merge(args: argparse.Namespace) -> int:
    input_paths = [Path(path) for path in (args.counts or [])]
    if args.input_dir:
        folder = Path(args.input_dir)
        input_paths.extend(sorted(folder.glob("*.txt")))
        input_paths.extend(sorted(folder.glob("*.tsv")))
        input_paths.extend(sorted(folder.glob("*.csv")))
    input_paths = [path for path in input_paths if path.exists() and path.is_file()]
    if not input_paths:
        raise StudentError("No featureCounts files or count matrix files were found.")

    def work(outdir: Path) -> None:
        merged = None
        mapping = []
        summaries = []
        for path in input_paths:
            df = read_table(path)
            cleaned = clean_counts(df)
            errors = validate_counts_frame(cleaned)
            if errors:
                raise StudentError(f"{path} is not a valid count file: {'; '.join(errors)}")
            sample_cols = [col for col in cleaned.columns if col != "gene_id"]
            renamed = cleaned[["gene_id", *sample_cols]].copy()
            rename_map = {col: _sample_name_from_count_column(col) for col in sample_cols}
            renamed = renamed.rename(columns=rename_map)
            for old, new in rename_map.items():
                mapping.append({"source_file": str(path), "source_column": old, "sample_id": new})
            if merged is None:
                merged = renamed
            else:
                if not merged["gene_id"].equals(renamed["gene_id"]):
                    raise StudentError(f"Gene IDs in {path} do not match previous count files.")
                overlap = set(merged.columns) & set(renamed.columns) - {"gene_id"}
                if overlap:
                    raise StudentError(f"Duplicate sample IDs after merging: {', '.join(sorted(overlap))}")
                merged = merged.merge(renamed, on="gene_id", how="inner", validate="one_to_one")
            summaries.append({"file": str(path), "genes": int(cleaned.shape[0]), "samples": len(sample_cols)})
        assert merged is not None
        write_csv(merged, outdir / "merged_raw_counts.csv")
        pd.DataFrame(mapping).to_csv(outdir / "sample_column_mapping.csv", index=False)
        (outdir / "merge_summary.json").write_text(
            json.dumps({"input_files": summaries, "genes": int(merged.shape[0]), "samples": int(merged.shape[1] - 1)}, indent=2) + "\n",
            encoding="utf-8",
        )

    return run_atomic_step(
        step_name="MERGE_FEATURECOUNTS",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={},
        inputs={"input_dir": args.input_dir or "", "counts": ",".join(str(path) for path in input_paths)},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_validate(args: argparse.Namespace) -> int:
    counts = require_file(args.counts, "Counts matrix")
    metadata = require_file(args.metadata, "Sample metadata")
    contrasts = Path(args.contrasts).resolve() if args.contrasts else None
    gene_mapping = Path(args.gene_mapping).resolve() if args.gene_mapping else None

    def work(outdir: Path) -> None:
        count_df = clean_counts(read_table(counts))
        errors = validate_counts_frame(count_df)
        meta_df = read_table(metadata)
        sample_ids = [col for col in count_df.columns if col != "gene_id"]
        errors.extend(validate_metadata(meta_df, sample_ids, args.group_column))
        warnings = []
        total_counts = {sample: int(count_df[sample].sum()) for sample in sample_ids}
        zero_samples = [sample for sample, total in total_counts.items() if total == 0]
        if zero_samples:
            warnings.append(f"Zero-library samples: {', '.join(zero_samples)}")
        contrast_rows = []
        if contrasts and contrasts.exists():
            contrast_df = read_table(contrasts)
            required = {"contrast_id", "numerator", "denominator", "design_column"}
            missing = required - set(contrast_df.columns)
            if missing:
                errors.append(f"Contrasts file is missing columns: {', '.join(sorted(missing))}")
            else:
                for _, row in contrast_df.iterrows():
                    column = row["design_column"]
                    if column not in meta_df.columns:
                        errors.append(f"Contrast {row['contrast_id']} uses missing design column '{column}'.")
                        continue
                    levels = set(meta_df[column].astype(str))
                    for level in (str(row["numerator"]), str(row["denominator"])):
                        if level not in levels:
                            errors.append(f"Contrast {row['contrast_id']} references absent level '{level}' in {column}.")
                    contrast_rows.append(row.to_dict())
        else:
            levels = sorted(meta_df[args.group_column].astype(str).unique()) if args.group_column in meta_df.columns else []
            if len(levels) >= 2:
                contrast_rows.append({
                    "contrast_id": f"{levels[1]}_vs_{levels[0]}",
                    "numerator": levels[1],
                    "denominator": levels[0],
                    "design_column": args.group_column,
                })
            else:
                warnings.append("No contrasts provided and fewer than two groups were found.")
        if errors:
            raise StudentError("; ".join(errors))
        count_df = count_df[["gene_id", *sample_ids]]
        meta_df = meta_df.set_index("sample_id").loc[sample_ids].reset_index()
        write_csv(count_df, outdir / "validated_counts.csv")
        write_csv(meta_df, outdir / "validated_metadata.csv")
        pd.DataFrame(contrast_rows).to_csv(outdir / "validated_contrasts.csv", index=False)
        pd.DataFrame({"sample_id": sample_ids, "count_matrix_column": sample_ids}).to_csv(outdir / "sample_id_mapping.csv", index=False)
        if gene_mapping and gene_mapping.exists():
            gm = read_table(gene_mapping)
            gm.to_csv(outdir / "gene_mapping_validation.csv", index=False)
        else:
            pd.DataFrame({"gene_id": count_df["gene_id"], "mapping_status": "not_provided"}).to_csv(
                outdir / "gene_mapping_validation.csv", index=False
            )
        report = {
            "valid": True,
            "genes": int(count_df.shape[0]),
            "samples": int(len(sample_ids)),
            "total_counts_per_sample": total_counts,
            "zero_library_samples": zero_samples,
            "missing_values": int(count_df.isna().sum().sum() + meta_df.isna().sum().sum()),
            "duplicate_genes": [],
            "duplicate_samples": [],
            "condition_counts": meta_df[args.group_column].astype(str).value_counts().to_dict(),
            "warnings": warnings,
            "errors": [],
        }
        (outdir / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return run_atomic_step(
        step_name="VALIDATE_DOWNSTREAM_INPUTS",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"group_column": args.group_column},
        inputs={"counts": str(counts), "metadata": str(metadata), "contrasts": str(contrasts or ""), "gene_mapping": str(gene_mapping or "")},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_filter(args: argparse.Namespace) -> int:
    counts = require_file(args.counts, "Validated counts")
    metadata = require_file(args.metadata, "Validated metadata")

    def work(outdir: Path) -> None:
        count_df = read_table(counts)
        meta_df = read_table(metadata)
        sample_ids = [col for col in count_df.columns if col != "gene_id"]
        group_sizes = meta_df[args.group_column].astype(str).value_counts()
        min_group_n = int(group_sizes.min()) if not group_sizes.empty else 1
        min_samples = int(args.min_samples or min_group_n)
        count_values = count_df[sample_ids].astype(float)
        lib_sizes = count_values.sum(axis=0)
        cpm = count_values.div(lib_sizes.replace(0, np.nan), axis=1) * 1_000_000
        keep_mask = (cpm >= float(args.cpm)).sum(axis=1) >= min_samples
        filtered = count_df.loc[keep_mask].copy()
        removed = count_df.loc[~keep_mask, ["gene_id"]].copy()
        removed["reason"] = f"CPM < {args.cpm} in fewer than {min_samples} samples"
        write_csv(filtered, outdir / "filtered_counts.csv")
        write_csv(removed, outdir / "removed_genes.csv")
        pd.DataFrame([{
            "input_genes": int(count_df.shape[0]),
            "retained_genes": int(filtered.shape[0]),
            "removed_genes": int(removed.shape[0]),
            "cpm_threshold": float(args.cpm),
            "min_samples": min_samples,
        }]).to_csv(outdir / "gene_filter_summary.csv", index=False)
        (outdir / "filter_parameters.json").write_text(
            json.dumps({"cpm": float(args.cpm), "min_samples": min_samples, "group_column": args.group_column}, indent=2) + "\n",
            encoding="utf-8",
        )

    return run_atomic_step(
        step_name="FILTER_LOW_EXPRESSION",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"cpm": args.cpm, "min_samples": args.min_samples, "group_column": args.group_column},
        inputs={"counts": str(counts), "metadata": str(metadata)},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_normalize(args: argparse.Namespace) -> int:
    counts = require_file(args.counts, "Filtered counts")
    metadata = require_file(args.metadata, "Validated metadata")

    def work(outdir: Path) -> None:
        count_df = read_table(counts)
        sample_ids = [col for col in count_df.columns if col != "gene_id"]
        values = count_df[sample_ids].astype(float)
        lib_sizes = values.sum(axis=0)
        median_lib = float(np.median(lib_sizes[lib_sizes > 0])) if (lib_sizes > 0).any() else 1.0
        size_factors = (lib_sizes / median_lib).replace(0, np.nan).fillna(1.0)
        normalized = values.div(size_factors, axis=1)
        vst = np.log2(normalized + 1.0)
        normalized_df = pd.concat([count_df[["gene_id"]], normalized.round(6)], axis=1)
        vst_df = pd.concat([count_df[["gene_id"]], vst.round(6)], axis=1)
        write_csv(normalized_df, outdir / "normalized_counts.csv")
        write_csv(vst_df, outdir / "vst_expression.csv")
        pd.DataFrame({"sample_id": sample_ids, "size_factor": [float(size_factors[s]) for s in sample_ids]}).to_csv(
            outdir / "size_factors.csv", index=False
        )
        (outdir / "transformation_summary.json").write_text(
            json.dumps({
                "method": "median_library_size_log2_plus_one",
                "note": "Python fallback for CLI/test mode. Use DESeq2 VST for production DE workflows.",
                "metadata": str(metadata),
            }, indent=2) + "\n",
            encoding="utf-8",
        )

    return run_atomic_step(
        step_name="NORMALIZE_AND_TRANSFORM",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"method": args.method},
        inputs={"counts": str(counts), "metadata": str(metadata)},
        work=work,
        repo_root=REPO_ROOT,
    )


def _expression_matrix(expression: Path) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    expr = read_table(expression)
    sample_ids = [col for col in expr.columns if col != "gene_id"]
    x = expr[sample_ids].astype(float).T.to_numpy()
    return expr, sample_ids, x


def cmd_pca(args: argparse.Namespace) -> int:
    expression = require_file(args.expression, "VST expression")
    metadata = require_file(args.metadata, "Validated metadata")

    def work(outdir: Path) -> None:
        expr, sample_ids, x = _expression_matrix(expression)
        meta = read_table(metadata).set_index("sample_id").loc[sample_ids].reset_index()
        n_components = min(int(args.components), x.shape[0], x.shape[1])
        try:
            from sklearn.decomposition import PCA

            model = PCA(n_components=n_components, random_state=int(args.seed))
            scores = model.fit_transform(x)
            loadings = model.components_.T
            explained = model.explained_variance_ratio_
            method = "sklearn_pca"
        except Exception:
            centered = x - x.mean(axis=0)
            u, s, vt = np.linalg.svd(centered, full_matrices=False)
            scores = u[:, :n_components] * s[:n_components]
            loadings = vt[:n_components, :].T
            variance = (s ** 2) / max(1, x.shape[0] - 1)
            explained = variance[:n_components] / variance.sum() if variance.sum() else np.zeros(n_components)
            method = "numpy_svd_fallback"
        score_df = meta.copy()
        for idx in range(n_components):
            score_df[f"PC{idx + 1}"] = scores[:, idx]
        loading_df = pd.DataFrame(loadings, columns=[f"PC{i + 1}" for i in range(n_components)])
        loading_df.insert(0, "gene_id", expr["gene_id"])
        variance_df = pd.DataFrame({
            "component": [f"PC{i + 1}" for i in range(n_components)],
            "explained_variance_ratio": explained,
        })
        write_csv(score_df, outdir / "pca_scores.csv")
        write_csv(loading_df, outdir / "pca_loadings.csv")
        write_csv(variance_df, outdir / "pca_explained_variance.csv")
        (outdir / "pca_parameters.json").write_text(json.dumps({"components": n_components, "seed": int(args.seed), "method": method}, indent=2) + "\n")

    return run_atomic_step(
        step_name="PCA",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"components": args.components, "seed": args.seed},
        inputs={"expression": str(expression), "metadata": str(metadata)},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_umap(args: argparse.Namespace) -> int:
    expression = require_file(args.expression, "VST expression")
    metadata = require_file(args.metadata, "Validated metadata")

    def work(outdir: Path) -> None:
        expr, sample_ids, x = _expression_matrix(expression)
        meta = read_table(metadata).set_index("sample_id").loc[sample_ids].reset_index()
        warnings = []
        if x.shape[0] < 4:
            warnings.append("Too few samples for stable UMAP; using first two centered expression dimensions.")
            coords = x[:, :2] if x.shape[1] >= 2 else np.column_stack([x[:, 0], np.zeros(x.shape[0])])
            coords = coords - coords.mean(axis=0)
        else:
            try:
                import umap

                model = umap.UMAP(
                    n_neighbors=min(int(args.n_neighbors), max(2, x.shape[0] - 1)),
                    min_dist=float(args.min_dist),
                    metric=args.metric,
                    random_state=int(args.seed),
                )
                coords = model.fit_transform(x)
            except Exception as exc:
                warnings.append(f"umap-learn unavailable or failed ({exc}); using deterministic centered-expression fallback coordinates.")
                centered = x - x.mean(axis=0)
                if centered.shape[1] >= 2:
                    coords = centered[:, :2]
                else:
                    coords = np.column_stack([centered[:, 0], np.zeros(centered.shape[0])])
        out = meta.copy()
        out["UMAP1"] = coords[:, 0]
        out["UMAP2"] = coords[:, 1]
        write_csv(out, outdir / "umap_coordinates.csv")
        (outdir / "umap_parameters.json").write_text(
            json.dumps({
                "n_neighbors": int(args.n_neighbors),
                "min_dist": float(args.min_dist),
                "metric": args.metric,
                "seed": int(args.seed),
                "warnings": warnings,
            }, indent=2) + "\n",
            encoding="utf-8",
        )

    return run_atomic_step(
        step_name="UMAP",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"n_neighbors": args.n_neighbors, "min_dist": args.min_dist, "metric": args.metric, "seed": args.seed},
        inputs={"expression": str(expression), "metadata": str(metadata)},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_plsda(args: argparse.Namespace) -> int:
    expression = require_file(args.expression, "VST expression")
    metadata = require_file(args.metadata, "Validated metadata")

    def work(outdir: Path) -> None:
        expr, sample_ids, x = _expression_matrix(expression)
        meta = read_table(metadata).set_index("sample_id").loc[sample_ids].reset_index()
        if args.group_column not in meta.columns:
            raise StudentError(f"Metadata is missing group column '{args.group_column}'.")
        y_labels = meta[args.group_column].astype(str)
        if y_labels.nunique() < 2:
            raise StudentError("PLS-DA requires at least two groups.")
        try:
            from sklearn.cross_decomposition import PLSRegression
            from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
            from sklearn.model_selection import StratifiedKFold, cross_val_predict
            from sklearn.preprocessing import LabelEncoder

            encoder = LabelEncoder()
            y = encoder.fit_transform(y_labels)
            components = min(int(args.components), x.shape[0] - 1, x.shape[1])
            model = PLSRegression(n_components=max(1, components))
            scores = model.fit_transform(x, y)[0]
            loadings_values = model.x_loadings_
            vip = np.sqrt((model.x_weights_ ** 2).sum(axis=1))
            folds = min(3, min(pd.Series(y).value_counts()))
            metrics = {"method": "sklearn_plsda", "warning": "PLS-DA is exploratory and can overfit, especially with small sample sizes."}
            if folds >= 2:
                cv = StratifiedKFold(n_splits=int(folds), shuffle=True, random_state=int(args.seed))
                pred = cross_val_predict(model, x, y, cv=cv)
                pred_class = np.rint(pred).astype(int).clip(0, len(encoder.classes_) - 1).ravel()
                metrics.update({
                    "accuracy": float(accuracy_score(y, pred_class)),
                    "balanced_accuracy": float(balanced_accuracy_score(y, pred_class)),
                    "f1_macro": float(f1_score(y, pred_class, average="macro")),
                    "folds": int(folds),
                })
            else:
                metrics["cross_validation_warning"] = "Not enough samples per group for stratified CV."
        except Exception as exc:
            classes = sorted(y_labels.unique())
            group_a = x[y_labels == classes[0]]
            group_b = x[y_labels == classes[1]]
            direction = group_b.mean(axis=0) - group_a.mean(axis=0)
            denom = np.linalg.norm(direction) or 1.0
            unit = direction / denom
            scores = (x @ unit).reshape(-1, 1)
            loadings_values = unit.reshape(-1, 1)
            vip = np.abs(unit)
            metrics = {
                "method": "numpy_group_projection_fallback",
                "warning": "scikit-learn unavailable; wrote deterministic exploratory group-projection fallback, not a full PLS-DA model.",
                "fallback_reason": str(exc),
            }
        score_df = meta.copy()
        for idx in range(scores.shape[1]):
            score_df[f"LV{idx + 1}"] = scores[:, idx]
        loadings = pd.DataFrame(loadings_values, columns=[f"LV{i + 1}" for i in range(loadings_values.shape[1])])
        loadings.insert(0, "gene_id", expr["gene_id"])
        pd.DataFrame({"gene_id": expr["gene_id"], "vip_score": vip}).to_csv(outdir / "plsda_vip_scores.csv", index=False)
        write_csv(score_df, outdir / "plsda_scores.csv")
        write_csv(loadings, outdir / "plsda_loadings.csv")
        pd.DataFrame([metrics]).to_csv(outdir / "plsda_cross_validation.csv", index=False)
        pd.DataFrame([{"status": "not_run", "reason": "Permutation test placeholder for first implementation."}]).to_csv(
            outdir / "plsda_permutation_test.csv", index=False
        )
        (outdir / "plsda_summary.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    return run_atomic_step(
        step_name="PLSDA",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"components": args.components, "group_column": args.group_column, "seed": args.seed},
        inputs={"expression": str(expression), "metadata": str(metadata)},
        work=work,
        repo_root=REPO_ROOT,
    )


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        original_rank = n - rank + 1
        value = pvalues[idx] * n / max(1, original_rank)
        running = min(running, value)
        adjusted[idx] = min(1.0, running)
    return adjusted


def cmd_differential_expression(args: argparse.Namespace) -> int:
    counts = require_file(args.counts, "Filtered counts")
    metadata = require_file(args.metadata, "Validated metadata")
    contrasts = require_file(args.contrasts, "Validated contrasts")
    gene_mapping = Path(args.gene_mapping).resolve() if args.gene_mapping else None

    def work(outdir: Path) -> None:
        count_df = read_table(counts)
        meta = read_table(metadata)
        contrast_df = read_table(contrasts)
        gene_symbols = {}
        if gene_mapping and gene_mapping.exists():
            gm = read_table(gene_mapping)
            if {"gene_id", "gene_symbol"}.issubset(gm.columns):
                gene_symbols = dict(zip(gm["gene_id"].astype(str), gm["gene_symbol"].astype(str)))
        sample_cols = [col for col in count_df.columns if col != "gene_id"]
        values = count_df[sample_cols].astype(float)
        log_values = np.log2(values + 1.0)
        meta_indexed = meta.set_index("sample_id")
        for _, contrast in contrast_df.iterrows():
            contrast_id = str(contrast["contrast_id"])
            numerator = str(contrast["numerator"])
            denominator = str(contrast["denominator"])
            design_column = str(contrast["design_column"])
            if design_column not in meta_indexed.columns:
                raise StudentError(f"Contrast {contrast_id} uses missing design column '{design_column}'.")
            numerator_samples = [sample for sample in sample_cols if str(meta_indexed.loc[sample, design_column]) == numerator]
            denominator_samples = [sample for sample in sample_cols if str(meta_indexed.loc[sample, design_column]) == denominator]
            if len(numerator_samples) < 2 or len(denominator_samples) < 2:
                raise StudentError(f"Contrast {contrast_id} requires at least two samples per group for fallback DE.")
            num = log_values[numerator_samples]
            den = log_values[denominator_samples]
            num_mean = num.mean(axis=1)
            den_mean = den.mean(axis=1)
            num_var = num.var(axis=1, ddof=1)
            den_var = den.var(axis=1, ddof=1)
            denom = np.sqrt((num_var / len(numerator_samples)) + (den_var / len(denominator_samples)))
            stat = np.divide(num_mean - den_mean, denom, out=np.zeros_like(num_mean), where=denom > 0)
            # Normal approximation keeps this fallback dependency-free. DESeq2 remains the production target.
            pvalue = np.array([math.erfc(abs(float(value)) / math.sqrt(2.0)) for value in stat])
            padj = _bh_adjust(pvalue)
            numerator_mean = values[numerator_samples].mean(axis=1)
            denominator_mean = values[denominator_samples].mean(axis=1)
            log2fc = np.log2((numerator_mean + 1.0) / (denominator_mean + 1.0))
            base_mean = values.mean(axis=1)
            result = pd.DataFrame({
                "gene_id": count_df["gene_id"].astype(str),
                "gene_symbol": [gene_symbols.get(str(gene), "") for gene in count_df["gene_id"].astype(str)],
                "baseMean": base_mean,
                "log2FoldChange": log2fc,
                "lfcSE": np.nan,
                "stat": stat,
                "pvalue": pvalue,
                "padj": padj,
            })
            result["significant"] = (result["padj"] <= float(args.padj)) & (result["log2FoldChange"].abs() >= float(args.log2fc))
            contrast_dir = outdir / "results" / contrast_id
            contrast_dir.mkdir(parents=True, exist_ok=True)
            result.sort_values(["padj", "pvalue", "gene_id"]).to_csv(contrast_dir / "deseq2_results.csv", index=False)
            result.loc[result["significant"]].to_csv(contrast_dir / "significant_genes.csv", index=False)
            result.sort_values("stat", ascending=False).to_csv(contrast_dir / "ranked_genes.csv", index=False)
            (contrast_dir / "deseq2_summary.json").write_text(
                json.dumps({
                    "contrast_id": contrast_id,
                    "method": "numpy_welch_normal_approx_fallback",
                    "note": "Testing/teaching fallback with DESeq2-compatible output columns; use DESeq2 for production inference.",
                    "numerator": numerator,
                    "denominator": denominator,
                    "design_column": design_column,
                    "genes": int(result.shape[0]),
                    "significant_genes": int(result["significant"].sum()),
                    "padj_threshold": float(args.padj),
                    "abs_log2fc_threshold": float(args.log2fc),
                }, indent=2) + "\n",
                encoding="utf-8",
            )

    return run_atomic_step(
        step_name="DIFFERENTIAL_EXPRESSION",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"padj": args.padj, "log2fc": args.log2fc, "method": "numpy_welch_normal_approx_fallback"},
        inputs={"counts": str(counts), "metadata": str(metadata), "contrasts": str(contrasts), "gene_mapping": str(gene_mapping or "")},
        work=work,
        repo_root=REPO_ROOT,
    )


def _hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    # P[X >= k] with only stdlib math, suitable for small local test sets.
    denom = math.comb(population, draws)
    total = 0
    for i in range(k, min(successes, draws) + 1):
        total += math.comb(successes, i) * math.comb(population - successes, draws - i)
    return float(total / denom) if denom else 1.0


def _bh_adjust_list(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda item: item[1], reverse=True)
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank_from_high, (index, pvalue) in enumerate(indexed, start=1):
        rank = total - rank_from_high + 1
        running = min(running, float(pvalue) * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _strip_ensembl_version(value: str) -> str:
    text = str(value).strip()
    if text.startswith("ENS") and "." in text:
        return text.split(".", 1)[0]
    return text


def _mapping_column(id_type: str) -> str:
    normalized = str(id_type or "").lower().replace("-", "_")
    columns = {
        "gene_id": "gene_id",
        "ensembl": "ensembl_id",
        "ensembl_id": "ensembl_id",
        "symbol": "gene_symbol",
        "gene_symbol": "gene_symbol",
        "entrez": "entrez_id",
        "entrez_id": "entrez_id",
    }
    return columns.get(normalized, normalized)


def _load_id_mapping(path: Path | None, input_id_type: str, gene_set_id_type: str) -> dict[str, str]:
    if not path:
        return {}
    mapping = read_table(path)
    source_col = _mapping_column(input_id_type)
    target_col = _mapping_column(gene_set_id_type)
    missing = [col for col in (source_col, target_col) if col not in mapping.columns]
    if missing:
        raise StudentError(
            "Gene mapping file is missing required column(s): "
            + ", ".join(missing)
            + ". Expected columns like gene_id,gene_symbol,ensembl_id,entrez_id."
        )
    result = {}
    for _, row in mapping.iterrows():
        source = _strip_ensembl_version(str(row[source_col]))
        target = str(row[target_col]).strip()
        if source and source.lower() != "nan" and target and target.lower() != "nan":
            result[source] = target
    return result


def _convert_gene_ids(genes: set[str], mapping: dict[str, str]) -> tuple[set[str], list[dict[str, str]]]:
    if not mapping:
        return {_strip_ensembl_version(gene) for gene in genes}, []
    converted = set()
    unmapped = []
    for gene in genes:
        normalized = _strip_ensembl_version(gene)
        mapped = mapping.get(normalized)
        if mapped:
            converted.add(mapped)
        else:
            unmapped.append({"gene_id": gene, "normalized_gene_id": normalized, "reason": "not found in mapping file"})
    return converted, unmapped


def _convert_ranked_index(ranked: pd.Series, mapping: dict[str, str]) -> tuple[pd.Series, list[dict[str, str]]]:
    if not mapping:
        ranked = ranked.copy()
        ranked.index = [_strip_ensembl_version(gene) for gene in ranked.index]
        return ranked, []
    rows = []
    unmapped = []
    for gene, value in ranked.items():
        normalized = _strip_ensembl_version(str(gene))
        mapped = mapping.get(normalized)
        if mapped:
            rows.append((mapped, value))
        else:
            unmapped.append({"gene_id": str(gene), "normalized_gene_id": normalized, "reason": "not found in mapping file"})
    if not rows:
        return pd.Series(dtype=float), unmapped
    return pd.Series({gene: value for gene, value in rows}, dtype=float), unmapped


def _run_overrepresentation(genes: set[str], universe: set[str], gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    query = genes & universe
    for name, members in gene_sets.items():
        tested_members = members & universe
        overlap = query & tested_members
        if not overlap:
            continue
        pvalue = _hypergeom_sf(len(overlap), len(universe), len(tested_members), len(query))
        rows.append({
            "term": name,
            "overlap_count": len(overlap),
            "query_size": len(query),
            "term_size": len(tested_members),
            "universe_size": len(universe),
            "pvalue": pvalue,
            "genes": ";".join(sorted(overlap)),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["padj"] = _bh_adjust_list(df["pvalue"].astype(float).tolist())
        df = df.sort_values(["pvalue", "term"])
    return df


def cmd_enrichment(args: argparse.Namespace, step_name: str, prefix: str) -> int:
    significant = require_file(args.significant_genes, "Significant genes")
    ranked = require_file(args.ranked_genes, "Ranked genes")
    universe_path = require_file(args.universe, "Tested-gene universe")
    gmt = require_file(args.gmt, "GMT gene set file")
    gene_mapping = Path(args.gene_mapping).resolve() if args.gene_mapping else None

    def work(outdir: Path) -> None:
        sig = read_table(significant)
        universe_df = read_table(universe_path)
        gene_sets = read_gmt(gmt)
        gene_col = "gene_id" if "gene_id" in sig.columns else sig.columns[0]
        mapping = _load_id_mapping(gene_mapping, args.input_id_type, args.gene_set_id_type)
        universe, unmapped_universe = _convert_gene_ids(set(universe_df["gene_id"].astype(str)), mapping)
        lfc = sig["log2FoldChange"] if "log2FoldChange" in sig.columns else pd.Series([1] * sig.shape[0])
        up, unmapped_up = _convert_gene_ids(set(sig.loc[lfc >= 0, gene_col].astype(str)), mapping)
        down, unmapped_down = _convert_gene_ids(set(sig.loc[lfc < 0, gene_col].astype(str)), mapping)
        _run_overrepresentation(up, universe, gene_sets).to_csv(outdir / f"{prefix}_overrepresentation_up.csv", index=False)
        _run_overrepresentation(down, universe, gene_sets).to_csv(outdir / f"{prefix}_overrepresentation_down.csv", index=False)
        pd.DataFrame(unmapped_up + unmapped_down + unmapped_universe).to_csv(outdir / "unmapped_genes.csv", index=False)
        pd.DataFrame([{
            "input_gene_count": int(sig.shape[0]),
            "mapped_up": len(up),
            "mapped_down": len(down),
            "mapped_universe": len(universe),
            "gene_sets": len(gene_sets),
        }]).to_csv(outdir / "gene_id_mapping_summary.csv", index=False)
        (outdir / "gene_set_database_metadata.json").write_text(
            json.dumps({
                "gmt": str(gmt),
                "source": "local",
                "organism": args.organism,
                "input_id_type": args.input_id_type,
                "gene_set_id_type": args.gene_set_id_type,
                "gene_mapping": str(gene_mapping or ""),
                "checksum_recorded_in": "checksums.sha256",
            }, indent=2) + "\n",
            encoding="utf-8",
        )
        (outdir / f"{prefix}_enrichment_summary.json").write_text(
            json.dumps({"significant_genes": int(sig.shape[0]), "gene_sets": len(gene_sets), "ranked_genes": str(ranked), "organism": args.organism}, indent=2) + "\n",
            encoding="utf-8",
        )

    return run_atomic_step(
        step_name=step_name,
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"organism": args.organism, "input_id_type": args.input_id_type, "gene_set_id_type": args.gene_set_id_type},
        inputs={"significant_genes": str(significant), "ranked_genes": str(ranked), "universe": str(universe_path), "gmt": str(gmt), "gene_mapping": str(gene_mapping or "")},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_gsea(args: argparse.Namespace) -> int:
    ranked = require_file(args.ranked_genes, "Ranked genes")
    gmt = require_file(args.gmt, "GMT gene set file")
    gene_mapping = Path(args.gene_mapping).resolve() if args.gene_mapping else None

    def work(outdir: Path) -> None:
        ranked_df = read_table(ranked)
        if "gene_id" not in ranked_df.columns or args.rank_field not in ranked_df.columns:
            raise StudentError(f"Ranked genes must contain gene_id and {args.rank_field}.")
        gene_sets = read_gmt(gmt)
        mapping = _load_id_mapping(gene_mapping, args.input_id_type, args.gene_set_id_type)
        ranks, unmapped = _convert_ranked_index(ranked_df.set_index("gene_id")[args.rank_field].astype(float), mapping)
        rows = []
        leading = []
        for name, members in gene_sets.items():
            hits = [gene for gene in ranks.index if gene in members]
            if not hits:
                continue
            score = float(ranks.loc[hits].mean())
            rows.append({"term": name, "score": score, "genes_in_set": len(hits), "rank_field": args.rank_field})
            for gene in hits[:20]:
                leading.append({"term": name, "gene_id": gene})
        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            result_df = result_df.sort_values("score", ascending=False)
        result_df.to_csv(outdir / "gsea_results.csv", index=False)
        pd.DataFrame(leading).to_csv(outdir / "gsea_leading_edge_genes.csv", index=False)
        pd.DataFrame(unmapped).to_csv(outdir / "unmapped_genes.csv", index=False)
        (outdir / "gsea_parameters.json").write_text(json.dumps({"rank_field": args.rank_field, "seed": int(args.seed), "organism": args.organism, "input_id_type": args.input_id_type, "gene_set_id_type": args.gene_set_id_type}, indent=2) + "\n")
        (outdir / "gene_set_database_metadata.json").write_text(json.dumps({"gmt": str(gmt), "source": "local", "gene_mapping": str(gene_mapping or "")}, indent=2) + "\n")

    return run_atomic_step(
        step_name="RANKED_GSEA",
        outdir=Path(args.outdir),
        argv=sys.argv,
        params={"rank_field": args.rank_field, "seed": args.seed, "organism": args.organism, "input_id_type": args.input_id_type, "gene_set_id_type": args.gene_set_id_type},
        inputs={"ranked_genes": str(ranked), "gmt": str(gmt), "gene_mapping": str(gene_mapping or "")},
        work=work,
        repo_root=REPO_ROOT,
    )


def cmd_state(args: argparse.Namespace) -> int:
    run_dir = Path(args.runs_dir) / args.run_id
    state_dir = run_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)
    state = {
        "schema_version": "1.0",
        "run_id": args.run_id,
        "status": args.status,
        "current_step": args.current_step,
        "completed_steps": args.completed_steps or [],
        "failed_step": args.failed_step,
        "error_message": args.error_message,
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    tmp = state_dir / "run_state.json.tmp"
    tmp.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    tmp.replace(state_dir / "run_state.json")
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    state = Path(args.runs_dir) / args.run_id / "state" / "run_state.json"
    if not state.exists():
        raise StudentError(f"No run_state.json found for run '{args.run_id}' in {args.runs_dir}.")
    print(state.read_text(encoding="utf-8"))
    return 0


def not_implemented_r(args: argparse.Namespace, name: str) -> int:
    raise StudentError(f"{name} requires the R implementation and dependencies. This entry point is reserved; no plots are generated.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="SurvOm downstream bulk RNA-seq CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("merge-featurecounts")
    p.add_argument("--input-dir")
    p.add_argument("--counts", action="append")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("validate")
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--contrasts")
    p.add_argument("--gene-mapping")
    p.add_argument("--group-column", default="condition")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("filter")
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group-column", default="condition")
    p.add_argument("--cpm", default=1.0, type=float)
    p.add_argument("--min-samples", type=int)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_filter)

    p = sub.add_parser("normalize")
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--method", default="vst")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_normalize)

    p = sub.add_parser("pca")
    p.add_argument("--expression", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--components", default=5, type=int)
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_pca)

    p = sub.add_parser("umap")
    p.add_argument("--expression", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--n-neighbors", default=15, type=int)
    p.add_argument("--min-dist", default=0.1, type=float)
    p.add_argument("--metric", default="euclidean")
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_umap)

    p = sub.add_parser("plsda")
    p.add_argument("--expression", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group-column", default="condition")
    p.add_argument("--components", default=2, type=int)
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_plsda)

    p = sub.add_parser("differential-expression")
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--contrasts", required=True)
    p.add_argument("--gene-mapping")
    p.add_argument("--padj", default=0.05, type=float)
    p.add_argument("--log2fc", default=1.0, type=float)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_differential_expression)

    p = sub.add_parser("go-enrichment")
    p.add_argument("--significant-genes", required=True)
    p.add_argument("--ranked-genes", required=True)
    p.add_argument("--universe", required=True)
    p.add_argument("--gmt", required=True)
    p.add_argument("--gene-mapping")
    p.add_argument("--organism", default="human")
    p.add_argument("--input-id-type", default="ensembl_id")
    p.add_argument("--gene-set-id-type", default="ensembl_id")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=lambda args: cmd_enrichment(args, "GO_ENRICHMENT", "go"))

    p = sub.add_parser("pathway-enrichment")
    p.add_argument("--significant-genes", required=True)
    p.add_argument("--ranked-genes", required=True)
    p.add_argument("--universe", required=True)
    p.add_argument("--gmt", required=True)
    p.add_argument("--gene-mapping")
    p.add_argument("--organism", default="human")
    p.add_argument("--input-id-type", default="ensembl_id")
    p.add_argument("--gene-set-id-type", default="ensembl_id")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=lambda args: cmd_enrichment(args, "PATHWAY_ENRICHMENT", "pathway"))

    p = sub.add_parser("gsea")
    p.add_argument("--ranked-genes", required=True)
    p.add_argument("--gmt", required=True)
    p.add_argument("--gene-mapping")
    p.add_argument("--organism", default="human")
    p.add_argument("--input-id-type", default="ensembl_id")
    p.add_argument("--gene-set-id-type", default="ensembl_id")
    p.add_argument("--rank-field", default="stat")
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_gsea)

    p = sub.add_parser("deseq2")
    p.set_defaults(func=lambda args: not_implemented_r(args, "DIFFERENTIAL_EXPRESSION"))

    p = sub.add_parser("oplsda")
    p.set_defaults(func=lambda args: not_implemented_r(args, "OPLSDA"))

    p = sub.add_parser("write-state")
    p.add_argument("--run-id", required=True)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--status", required=True)
    p.add_argument("--current-step")
    p.add_argument("--completed-steps", action="append")
    p.add_argument("--failed-step")
    p.add_argument("--error-message")
    p.set_defaults(func=cmd_state)

    p = sub.add_parser("status")
    p.add_argument("--run-id", required=True)
    p.add_argument("--runs-dir", required=True)
    p.set_defaults(func=cmd_status)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except StudentError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
