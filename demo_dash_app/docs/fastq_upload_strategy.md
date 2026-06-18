# FASTQ Upload Strategy

## Recommendation

For production FASTQ uploads, use:

```text
Browser UI: Uppy Dashboard + @uppy/tus + Golden Retriever
Upload server: tusd standalone service
Workflow app: Dash/FastAPI-style backend receives completed-upload records
Storage: local filesystem now, MinIO/S3 later
```

Keep the existing Resumable.js uploader for the current teaching demo until the tus service is installed and tested. The Nextflow builder, job JSON, run directories, and output rendering do not need to change.

## Why tus for FASTQ

FASTQ files are commonly many GB. Uploads may run through campus Wi-Fi, VPN, or browser sessions that are interrupted. Tus is a better fit than a simple chunk endpoint because it has protocol-level resume. The client can ask the server how many bytes were already stored and continue from that offset.

Uppy also gives the user-facing upload experience: queue, progress, pause/resume, retry behavior, and multiple-file handling.

## Current Demo Upload Flow

```text
Browser
  -> Resumable.js chunks
  -> Dash Flask endpoints
  -> demo_dash_app/uploads/<session_id>/<upload_type>/
  -> run button creates run_request.yaml
  -> builder compiles Nextflow
  -> job runner executes run_command.sh
  -> outputs rendered from execution_record.json
```

This works for testing and small classes, but the Dash process is doing upload handling itself.

## Production Upload Flow

```text
Browser
  -> Uppy Dashboard
  -> @uppy/tus
  -> tusd /files endpoint
  -> tusd storage directory or S3/MinIO bucket
  -> tusd finish hook
  -> backend registers completed file
  -> demo_dash_app/uploads/<session_id>/<upload_type>/<original_filename>
  -> run button uses registered files
```

The backend should not hold the upload request open. It should only:

1. Authenticate the user/session.
2. Issue upload metadata or signed information if needed.
3. Receive completion events from tusd.
4. Store file records in the database.
5. Let the workflow runner use completed files only.

## Required Upload Metadata

Every tus upload should carry enough metadata for the backend to register it:

```text
session_id: browser/user session id
upload_type: fastq or metadata
original_filename: original user filename
sample_id: optional, if known at upload time
project_id: optional, useful in production
user_id: optional, useful with authentication
```

The browser should reject obvious wrong files before upload:

```text
FASTQ: .fastq, .fastq.gz, .fq, .fq.gz
Metadata: .csv, .tsv, .xlsx
```

The backend must still validate again after upload.

## Session and Collision Rules

Two users uploading at the same time must never write to the same logical input directory.

Use this layout:

```text
uploads/<session_id>/fastq/
uploads/<session_id>/metadata/
runs/<run_id>/
survom-pipelines/workflows/generated/<run_id>/
```

Use a new `run_id` for every execution, even if the same user runs the same step twice.

## Long Upload and Long Run UI

Upload status and workflow status are different states.

Upload states:

```text
selected
uploading
paused
resuming
uploaded
failed
```

Workflow states:

```text
compiled
queued
running
succeeded
failed
cancelled
```

The UI should poll workflow job status every 30 seconds. Upload progress should come from Uppy events in the browser, not from the workflow poller.

## Suggested Uppy Client Shape

```javascript
const uppy = new Uppy({
  autoProceed: false,
  restrictions: {
    allowedFileTypes: [".fastq", ".fq", ".gz", ".csv", ".tsv", ".xlsx"],
  },
  meta: {
    session_id: sessionId,
    upload_type: "fastq",
  },
});

uppy
  .use(Dashboard, { inline: true, target: "#fastq-upload-target" })
  .use(GoldenRetriever, { serviceWorker: true })
  .use(Tus, {
    endpoint: "/files/",
    retryDelays: [0, 1000, 3000, 5000, 10000],
    limit: 2,
    allowedMetaFields: [
      "session_id",
      "upload_type",
      "original_filename",
      "sample_id",
      "project_id",
      "user_id",
    ],
  });
```

For very large FASTQ uploads, keep concurrency low. Two simultaneous upload streams per browser is a sensible starting point.

## tusd Deployment Shape

Run tusd as its own process:

```text
Dash app: http://127.0.0.1:8050
tusd:     http://127.0.0.1:1080/files/
Nginx:    / uploads UI/API to Dash, /files/ to tusd
```

Nginx should allow large request bodies and long upload timeouts.

## Completed Upload Registration

When tusd finishes an upload, use a finish hook to call the backend and register the file.

The backend registration should create or update a durable record like:

```json
{
  "session_id": "session-or-user-id",
  "upload_type": "fastq",
  "original_filename": "sample_R1.fastq.gz",
  "stored_path": "/path/from/tusd/or/object/key",
  "size": 1234567890,
  "status": "uploaded"
}
```

For local filesystem mode, the hook can move or hard-link the completed tusd file into:

```text
demo_dash_app/uploads/<session_id>/<upload_type>/<original_filename>
```

For MinIO/S3 mode, do not move the file. Store the object key and make the workflow staging layer download or mount it for Nextflow.

## Implementation Order

1. Keep current Resumable.js demo working.
2. Install tusd and run it locally with disk storage.
3. Add Uppy/tus upload UI behind a feature flag.
4. Add a tusd finish hook that registers completed files.
5. Test with one small FASTQ, then interrupted upload resume, then two users at the same time.
6. Add MinIO/S3 storage only after local tusd is reliable.

## Do Not Change

The upload upgrade should not rewrite:

```text
survom-pipelines registry logic
Nextflow modules
builder compiler
run-job JSON contract
execution_record.json output rendering
```

Only the browser upload component and completed-file registration need to change.
