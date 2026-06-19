# SurvOm Dash + Nextflow Workflow Testing

This tutorial shows how to test RNA-seq atomic steps and supported chained workflows without opening the Dash app.

## 1. Activate the environment

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh
```

If that activation script is not available in your shell, use the known conda Python directly:

```bash
/home/vikash/miniconda3/envs/survom-preprocess/bin/python --version
```

## 2. Fast compile-only smoke test

This checks that each workflow can generate:

- `main.nf`
- `params.yaml`
- `workflow_plan.json`
- atomic workflow files

It does not run Docker or Nextflow processes.

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
/home/vikash/miniconda3/envs/survom-preprocess/bin/python scripts/smoke_workflows.py
```

Expected result:

```text
[compile ok] reference_only
[compile ok] raw_qc_only
[compile ok] trim_only
...
```

## 3. Test one workflow compile

```bash
/home/vikash/miniconda3/envs/survom-preprocess/bin/python scripts/smoke_workflows.py --case salmon_count_matrix
```

Useful cases:

```text
reference_only
raw_qc_only
trim_only
post_trim_qc_only
strandedness_only
salmon_quant_only
star_align_only
hisat2_align_only
salmon_count_matrix
star_count_matrix
qc_trim_strandedness
```

## 4. Run a real Nextflow smoke test

This executes the generated workflow with Docker/local profile:

```bash
/home/vikash/miniconda3/envs/survom-preprocess/bin/python scripts/smoke_workflows.py \
  --case salmon_count_matrix \
  --run \
  --profile local,docker \
  --keep /tmp/survom_salmon_smoke
```

Compiled workflows are written under the normal Nextflow generated-workflow area so relative module includes work correctly:

```text
/data/shared/vikash/mult-omics/survom-pipelines/workflows/generated/salmon_count_matrix/
```

Inputs and result outputs are kept under:

```text
/tmp/survom_salmon_smoke/
```

Check the execution record:

```bash
cat /data/shared/vikash/mult-omics/survom-pipelines/workflows/generated/salmon_count_matrix/execution_record.json
```

## 5. Manual compile from a run request

Dash creates `run_request.yaml` inside each run folder. You can recompile it manually:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
/home/vikash/miniconda3/envs/survom-preprocess/bin/python -m builder.api compile \
  --request /path/to/run_request.yaml \
  --output workflows/generated/debug_run
```

Then run:

```bash
cd workflows/generated/debug_run
nextflow run main.nf -params-file params.yaml -profile local,docker -work-dir work -resume
```

Collect outputs for Dash-like status:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
/home/vikash/miniconda3/envs/survom-preprocess/bin/python -m builder.api collect-outputs \
  --compiled-dir workflows/generated/debug_run \
  --status succeeded
```

## 6. How to read results

Workflow failure:

- Nextflow exits non-zero.
- `execution_record.json` has a failed step.
- Check `.command.err`, `.command.log`, and `.command.sh` links.

Workflow success with zero counts:

- Nextflow exits zero.
- Count matrix exists.
- Logs may show `0 mapped` or `0 assigned`.
- This means the tiny demo reference does not biologically match the reads well. It is not a pipeline crash.

## 7. Dash app sanity checks

Run this after editing callbacks or layout:

```bash
cd /data/shared/vikash/mult-omics/demo_dash_app
/home/vikash/miniconda3/envs/survom-preprocess/bin/python - <<'PY'
import app
ids=set()
def walk(x):
    if isinstance(x,(list,tuple)):
        for y in x: walk(y)
        return
    if hasattr(x,'id') and x.id is not None:
        ids.add(x.id)
    if hasattr(x,'children'):
        walk(x.children)
walk(app.app.layout)
missing=[]
for callback in app.app.callback_map.values():
    deps=[]
    deps.extend(callback.get('inputs',[]))
    deps.extend(callback.get('state',[]))
    out=callback.get('output',[])
    if not isinstance(out, list): out=[out]
    deps.extend(out)
    for dep in deps:
        cid = dep.get('id') if isinstance(dep,dict) else getattr(dep,'component_id',None)
        if isinstance(cid,str) and cid not in ids:
            missing.append(cid)
print('callbacks', len(app.app.callback_map), 'ids', len(ids), 'missing', sorted(set(missing)))
PY
```

Expected:

```text
missing []
```

## 8. Current organization notes

`app.py` still owns callbacks and high-level Dash layout. Display rendering for job status/logs now lives in:

```text
demo_dash_app/job_view.py
```

Good next modules to extract later:

- `dash_callbacks.py` for callback registration
- `workflow_ui.py` for parameter panels
- `validation_ui.py` for Dash-facing validation messages

Do those one at a time and run the callback sanity check after each move.
