# Run SurvOm Demo App

## Start the App

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh

cd /data/shared/vikash/mult-omics/demo_dash_app
python app.py
```

Keep this terminal open while students use the app.

## Link to Share

If students are on the same network as this server, share:

```text
http://obilab24:8056
```

If that does not open, find the server IP on the server terminal:

```bash
hostname -I
```

Then share:

```text
http://SERVER_IP:8056
```

If you are opening the app on the same computer where it runs:

```text
http://localhost:8056
```

## Test Files

FASTQ files:

```text
/data/shared/vikash/mult-omics/survom-pipelines/testdata/data/ggal/ggal_gut_1.fq
/data/shared/vikash/mult-omics/survom-pipelines/testdata/data/ggal/ggal_gut_2.fq
```

Sample sheet:

```text
/data/shared/vikash/mult-omics/survom-pipelines/testdata/ggal_gut_samplesheet.csv
```

## Student Flow

1. Open the link.
2. Upload one or more FASTQ files.
3. Upload the sample sheet.
4. Choose one step, such as FastQC, or click `Run full workflow`.
5. Wait for status updates. The app checks job status every 30 seconds.
6. Download outputs from the output panel.

## If Browser Shows Old JavaScript Error

Use a hard refresh:

```text
Ctrl + Shift + R
```

The old custom uploader script was removed. The app now uses Dash upload boxes only.
