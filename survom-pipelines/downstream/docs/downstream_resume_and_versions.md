# Downstream Resume And Versions

The downstream workflow is versioned in:

```text
survom-pipelines/downstream/VERSION
```

Each module records:

```text
workflow_version
python version
platform
git commit
git branch
dirty working tree flag
inputs
parameters
checksums.sha256
```

Nextflow resume:

```bash
nextflow -C nextflow_downstream.config run main_downstream.nf ... -resume
```

User-facing state:

```bash
python bin/downstream_cli.py write-state \
  --run-id downstream_demo_001 \
  --runs-dir /tmp/survom_downstream_runs \
  --status running \
  --current-step 05_pca

python bin/downstream_cli.py status \
  --run-id downstream_demo_001 \
  --runs-dir /tmp/survom_downstream_runs
```

The first implementation writes atomic metadata/checksums for each step and provides state/status commands. Incompatible resume enforcement is intentionally conservative and should be expanded when pinned Python/R/container lockfiles are finalized.
