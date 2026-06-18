# SurvOm Codex Agent Prompts

Use these prompts with Codex when you add steps gradually. They are designed to protect existing working modules.

---

## 1. General Codex System Prompt

```text
You are Codex working on SurvOm, a modular multiomics workflow builder using Nextflow DSL2.

Your job is to make safe, incremental changes.

Always preserve existing working steps, process names, registry contracts, output names, tests, and examples.

The architecture is:

registry/steps.yaml -> builder/planner.py -> builder/compiler.py -> generated Nextflow workflow -> reusable DSL2 modules

Do not create a monolithic pipeline.
Do not hard-code RNA-seq-only logic into COMMON modules.
Do not delete or rename working code unless explicitly asked and unless you add a migration note.

Before changing code, inspect the relevant registry entries, modules, planner/compiler code, and tests.

After changing code, show:
- files changed
- what changed
- what did not change
- compatibility impact
- how to run the step individually
- how to run it in a workflow chain
- tests added or updated
```

---

## 2. Prompt for Adding a New Step

```text
Add a new modular workflow step to SurvOm.

Step details:
- step_id: <STEP_ID>
- step_name: <STEP_NAME>
- category: <COMMON | RNA_SEQ_SPECIFIC | GENOMICS_COMMON | STATISTICS_COMMON | REPORTING_COMMON>
- omics: <omics list>
- default_tool: <tool>
- fallback_tool: <fallback or none>
- required inputs: <inputs>
- outputs: <outputs>
- recommended next steps: <next steps>

Constraints:
- Do not modify existing step behavior.
- Do not rename existing step IDs, process names, params, or outputs.
- Add registry entry first.
- Reuse existing modules when possible.
- If a module already exists under modules/common, do not duplicate it.
- Add tests for registry validation, planner next-step suggestion, and compiler generation.
- Add an example run request YAML.

Before coding, inspect:
- registry/steps.yaml
- modules/common
- modules/rnaseq
- builder/planner.py
- builder/compiler.py
- tests

After coding, show how to run:
- this step alone
- this step after its previous step
- a workflow chain including this step
```

---

## 3. Prompt for Modifying a Working Step Safely

```text
Modify an existing SurvOm workflow step safely.

Target step:
<STEP_ID>

Requested change:
<CHANGE>

Constraints:
- Preserve step_id.
- Preserve process_name unless explicitly required.
- Preserve existing output names and patterns unless this is a versioned migration.
- Preserve existing params; new params must have safe defaults.
- Existing examples must still compile.
- Existing tests must still pass.
- Add a migration note if any public contract changes.

Before coding, answer:
1. Is this additive or breaking?
2. Which registry fields are affected?
3. Which Nextflow emits are affected?
4. Which tests prove compatibility?

Then implement the smallest safe change.
```

---

## 4. Prompt for Creating a Reusable COMMON Module

```text
Create a reusable COMMON Nextflow DSL2 module for SurvOm.

Module:
<TOOL_OR_PROCESS>

Requirements:
- Place it under modules/common/<module_name>/main.nf.
- Keep it omics-agnostic.
- Accept typed tuple inputs, usually tuple val(meta), path(input_files).
- Emit stable named outputs.
- Emit versions.yml.
- Do not hard-code RNA-seq-specific metadata.
- Do not hard-code local or S3 paths.
- Configure containers/resources in nextflow.config.
- Add or update registry entries that use this module.
- Add compiler tests proving it can be included in generated workflows.

Show example usage for:
- one-step run
- chained workflow run
```

---

## 5. Prompt for Planner Changes

```text
Update the SurvOm workflow planner.

Requested behavior:
<BEHAVIOR>

Constraints:
- Planner must read from registry, not hard-code one RNA-seq order.
- Existing step suggestions must continue to work.
- Individual step execution must remain supported if required inputs are supplied.
- Invalid chains must return clear errors, not crashes.
- Add tests for recommended_next_steps and allowed_next_steps.

Do not change compiler behavior unless planner output format requires it.
```

---

## 6. Prompt for Compiler Changes

```text
Update the SurvOm workflow compiler.

Requested behavior:
<BEHAVIOR>

Constraints:
- Compiler must generate Nextflow DSL2 from selected registry steps.
- Do not duplicate module code into generated workflows.
- Include only modules required by selected steps.
- Preserve existing generated workflow behavior.
- Generated workflows must support local and AWS profiles through config, not separate code.
- Add snapshot tests or text assertions for generated main.nf and params.yaml.

Test at least:
- one-step workflow
- two-step chain
- existing example workflows
```

---

## 7. Prompt for Regression Review

```text
Review this SurvOm change for safety.

Check:
- Are any existing step IDs renamed?
- Are any process names renamed?
- Are any output contracts changed?
- Are any params removed or changed without defaults?
- Are COMMON modules still omics-agnostic?
- Do existing generated workflows still compile?
- Are tests updated?
- Is a migration note needed?

Return:
- safe / unsafe
- breaking changes found
- tests that must be added
- recommended fix
```

---

## 8. Prompt for Documentation Update

```text
Update SurvOm documentation for the changed workflow modules.

Constraints:
- Do not document unimplemented parameters.
- Separate user, expert, and internal parameters.
- Explain how to run the step individually.
- Explain how to run it in a chain.
- List outputs that UI/chat can show.
- Keep examples local/AWS compatible.
```
