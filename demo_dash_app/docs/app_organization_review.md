# Dash App Organization Review

Current status: `app.py` still owns layout, callbacks, route handlers, and run orchestration. It works, but it is too large for safe long-term development.

## What Was Improved

- Job/status rendering moved to `demo_dash_app/job_view.py`.
- The job screen now hides idle pending steps and shows concise step status.
- The raw Nextflow log is no longer shown as a large default block.
- `builder.outputs` now handles successful resumed workflows better:
  - completed outputs win over stale failed attempts
  - steps with outputs are not left pending after a successful run
  - task process names are parsed from the top of `.command.run`

## Remaining Risk Areas

- `app.py` still mixes UI layout, callback registration, validation, and Flask API routes.
- Parameter panels are still large and should be moved out next.
- Callback functions still depend on many globals.
- Workflow IDs are still partly hardcoded in route-selection logic.
- Dash UI and Nextflow compiler tests are separate; they should share more fixtures later.

## Safe Refactor Order

1. Keep `app.py` as the callback registry for now.
2. Move parameter panel builders into `workflow_ui.py`.
3. Move Dash-facing validation helpers into `validation_ui.py`.
4. Move Flask upload/job API routes into `api_routes.py`.
5. Only after those moves, consider splitting callbacks into `callbacks/*.py`.

Run these checks after each move:

```bash
cd /data/shared/vikash/mult-omics
/home/vikash/miniconda3/envs/survom-preprocess/bin/python -m py_compile \
  demo_dash_app/app.py \
  demo_dash_app/job_view.py \
  demo_dash_app/services/workflow_service.py \
  survom-pipelines/builder/outputs.py
```

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

## Command-Line Workflow Test

Use:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
/home/vikash/miniconda3/envs/survom-preprocess/bin/python scripts/smoke_workflows.py
```

For real Nextflow smoke tests:

```bash
/home/vikash/miniconda3/envs/survom-preprocess/bin/python scripts/smoke_workflows.py \
  --case salmon_count_matrix \
  --case star_count_matrix \
  --case qc_trim_strandedness \
  --run \
  --profile local,docker \
  --keep /tmp/survom_chain_smoke
```
