#!/usr/bin/env bash

# Source this file before compiling/running local preprocessing workflows:
# source scripts/activate_survom_preprocess.sh

source /home/vikash/miniconda3/etc/profile.d/conda.sh
conda activate survom-preprocess

# The conda environment currently provides OpenJDK 25, while Nextflow 24.10
# supports Java 8-22. Keep Nextflow on the system Java 17 and use conda tools
# for FastQC/FastP/Cutadapt.
export JAVA_CMD=/usr/bin/java
unset JAVA_HOME
unset JAVA_LD_LIBRARY_PATH

echo "Activated survom-preprocess"
echo "nextflow: $(command -v nextflow)"
echo "fastqc:   $(command -v fastqc)"
echo "fastp:    $(command -v fastp)"
echo "cutadapt: $(command -v cutadapt)"
echo "java:     ${JAVA_CMD}"
