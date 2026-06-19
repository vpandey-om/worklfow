from pathlib import Path

from builder.compiler import NextflowCompiler
from builder.models import RunRequest
from builder.registry import load_run_request, load_steps


ROOT = Path(__file__).resolve().parents[1]


def compiler():
    return NextflowCompiler(load_steps(ROOT / "registry" / "steps.yaml"), ROOT / "workflows" / "templates")


def make_fastqs(tmp_path):
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    r1.write_text("")
    r2.write_text("")
    return r1, r2


def make_raw_manifest(tmp_path):
    r1, r2 = make_fastqs(tmp_path)
    manifest = tmp_path / "raw_manifest.csv"
    manifest.write_text(
        "sample_id,fastq_1,fastq_2,single_end,strandedness,platform\n"
        f"sample,{r1},{r2},false,unknown,unknown\n"
    )
    return manifest


def make_trim_manifest(tmp_path):
    r1, r2 = make_fastqs(tmp_path)
    manifest = tmp_path / "trim_manifest.tsv"
    manifest.write_text(
        "sample_id\tfastq_1\tfastq_2\tsingle_end\tstrandedness\n"
        f"sample\t{r1}\t{r2}\tfalse\tunknown\n"
    )
    return manifest


def demo_reference_params():
    demo = ROOT / "assets" / "demo_reference"
    return {
        "reference_mode": "demo_reference",
        "selected_route": "salmon",
        "organism": "demo",
        "genome_build": "demo_build",
        "genome_fasta": str(demo / "genome.fa"),
        "gtf": str(demo / "genes.gtf"),
        "transcriptome_fasta": str(demo / "transcriptome.fa"),
        "tx2gene": str(demo / "tx2gene.tsv"),
        "salmon_index": str(demo / "salmon_index"),
        "star_index": str(demo / "star_index"),
        "strandedness_method": "salmon_auto",
    }


def test_compile_qc_only(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_qc_only.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "include { FASTQC }" in main
    assert "FASTQC(ch_reads)" in main
    assert "FASTP(ch_reads)" not in main


def test_compile_qc_trim(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_qc_trim.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    params = (out / "params.yaml").read_text()
    assert "include { FASTQC }" in main
    assert "include { FASTP }" in main
    assert "ch_trimmed_reads = FASTP.out.reads" in main
    assert "MAKE_PREPROCESS_MANIFEST(ch_trimmed_reads)" in main
    assert "MERGE_PREPROCESS_MANIFEST(MAKE_PREPROCESS_MANIFEST.out.fragments.collect())" in main
    assert (out / "samplesheet.resolved.tsv").exists()
    assert (out / "workflow_plan.json").exists()
    assert "quality_threshold: 20" in params


def test_compile_trim_only(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_trim_only.run.yaml")
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "include { FASTP }" in main
    assert "include { FASTQC }" not in main
    assert "FASTP(ch_reads)" in main


def test_compile_step_namespaced_custom_trimming_params(tmp_path):
    req = load_run_request(ROOT / "examples" / "rnaseq_trim_only.run.yaml")
    req.params = {
        "rnaseq_04_adapter_quality_trimming": {
            "trimming_tool": "fastp",
            "quality_threshold": 25,
            "minimum_read_length": 30,
            "trim_poly_g": "true",
            "trim_front_r1": 5,
        }
    }
    out = compiler().compile(req, tmp_path / "custom_trim")
    params = (out / "params.yaml").read_text()
    assert "quality_threshold: 25" in params
    assert "minimum_read_length: 30" in params
    assert 'trim_poly_g: "true"' in params
    assert "trim_front_r1: 5" in params


def test_compile_post_trim_qc_only_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    req = RunRequest(
        run_id="post_trim_qc_only",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=["rnaseq_05_post_trim_quality_control"],
        params={"trim_manifest": str(manifest), "trim_input_mode": "previous_manifest"},
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    resolved = (out / "input_manifest.resolved.tsv").read_text()
    assert "include { POST_TRIM_FASTQC }" in main
    assert "POST_TRIM_FASTQC(ch_trimmed_reads)" in main
    assert "fastq_1" in resolved
    assert "single_end" in resolved
    assert "strandedness" in resolved


def test_compile_strandedness_only_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    params = demo_reference_params()
    params.update({"trim_manifest": str(manifest), "trim_input_mode": "previous_manifest"})
    req = RunRequest(
        run_id="strandedness_only",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=["rnaseq_06a_reference_build_validation", "rnaseq_06_strandedness_inference"],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "REFERENCE_VALIDATE(" in main
    assert "SALMON_STRANDEDNESS_INFERENCE(ch_trimmed_reads, reference_salmon_index)" in main
    assert "FASTP(ch_reads)" not in main


def test_compile_salmon_quant_only_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    params = demo_reference_params()
    params.update({"trim_manifest": str(manifest), "trim_input_mode": "previous_manifest", "approved_library_type": "A"})
    req = RunRequest(
        run_id="salmon_quant_only",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=["rnaseq_06a_reference_build_validation", "rnaseq_07_salmon_quantification"],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "SALMON(ch_trimmed_reads, reference_salmon_index)" in main
    assert "FASTP(ch_reads)" not in main


def test_compile_star_alignment_only_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    params = demo_reference_params()
    params.update({"selected_route": "star", "trim_manifest": str(manifest), "trim_input_mode": "previous_manifest"})
    req = RunRequest(
        run_id="star_align_only",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=["rnaseq_06a_reference_build_validation", "rnaseq_09_star_alignment"],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "STAR_ALIGN_COUNTS(ch_trimmed_reads, reference_star_index)" in main
    assert "FASTP(ch_reads)" not in main


def test_compile_atomic_chain_raw_trim_post_trim_strandedness(tmp_path):
    manifest = make_raw_manifest(tmp_path)
    params = demo_reference_params()
    params.update({
        "trimming_tool": "fastp",
        "quality_threshold": 20,
        "minimum_read_length": 20,
        "trim_poly_g": "auto",
    })
    req = RunRequest(
        run_id="atomic_chain",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=[
            "rnaseq_03_raw_read_qc",
            "rnaseq_04_adapter_quality_trimming",
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06a_reference_build_validation",
            "rnaseq_06_strandedness_inference",
        ],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "FASTQC(ch_reads)" in main
    assert "FASTP(ch_reads)" in main
    assert "POST_TRIM_FASTQC(ch_trimmed_reads)" in main
    assert "REFERENCE_VALIDATE(" in main
    assert "SALMON_STRANDEDNESS_INFERENCE(ch_trimmed_reads, reference_salmon_index)" in main
    assert (out / "atomic_rnaseq_03_raw_read_qc.nf").exists()


def test_compile_chain_salmon_count_matrix_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    params = demo_reference_params()
    params.update({
        "trim_manifest": str(manifest),
        "trim_input_mode": "previous_manifest",
        "approved_library_type": "A",
        "selected_route": "salmon",
        "salmon_index": None,
    })
    req = RunRequest(
        run_id="salmon_count_matrix",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=[
            "rnaseq_06a_reference_build_validation",
            "rnaseq_07_salmon_quantification",
            "rnaseq_08_tximport_gene_summarization",
            "rnaseq_10_count_matrix_qc",
        ],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "BUILD_SALMON_INDEX(" in main
    assert "SALMON(ch_trimmed_reads, reference_salmon_index)" in main
    assert "TXIMPORT(SALMON.out.quant.map { meta, quant_dir -> quant_dir }.collect())" in main
    assert "COUNT_MATRIX_QC(TXIMPORT.out.counts)" in main
    assert "FASTP(ch_reads)" not in main
    assert "STAR_ALIGN_COUNTS" not in main


def test_compile_chain_star_count_matrix_from_trim_manifest(tmp_path):
    manifest = make_trim_manifest(tmp_path)
    params = demo_reference_params()
    params.update({
        "selected_route": "star",
        "trim_manifest": str(manifest),
        "trim_input_mode": "previous_manifest",
        "star_index": None,
    })
    req = RunRequest(
        run_id="star_count_matrix",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=[
            "rnaseq_06a_reference_build_validation",
            "rnaseq_09_star_alignment",
            "rnaseq_09b_bam_sort_index",
            "rnaseq_09c_featurecounts_gene_counting",
            "rnaseq_10_count_matrix_qc",
        ],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "BUILD_STAR_INDEX(" in main
    assert "STAR_ALIGN_COUNTS(ch_trimmed_reads, reference_star_index)" in main
    assert "BAM_SORT_INDEX(STAR_ALIGN_COUNTS.out.bam)" in main
    assert "FEATURECOUNTS(ch_featurecounts_alignments, file(params.gtf, checkIfExists: true))" in main
    assert "COUNT_MATRIX_QC(FEATURECOUNTS.out.counts)" in main
    assert "FASTP(ch_reads)" not in main
    assert "SALMON(ch_trimmed_reads" not in main


def test_compile_chain_qc_trim_strandedness_demo_mode(tmp_path):
    manifest = make_raw_manifest(tmp_path)
    params = demo_reference_params()
    params.update({
        "selected_route": "salmon",
        "salmon_index": None,
        "trimming_tool": "fastp",
        "quality_threshold": 20,
        "minimum_read_length": 20,
        "trim_poly_g": "auto",
    })
    req = RunRequest(
        run_id="qc_trim_strandedness",
        omics="rnaseq",
        execution_profile="local",
        input=str(manifest),
        outdir=str(tmp_path / "results"),
        selected_steps=[
            "rnaseq_03_raw_read_qc",
            "rnaseq_04_adapter_quality_trimming",
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06a_reference_build_validation",
            "rnaseq_06_strandedness_inference",
        ],
        params=params,
    )
    out = compiler().compile(req, tmp_path / req.run_id)
    main = (out / "main.nf").read_text()
    assert "FASTQC(ch_reads)" in main
    assert "FASTP(ch_reads)" in main
    assert "MERGE_PREPROCESS_MANIFEST(MAKE_PREPROCESS_MANIFEST.out.fragments.collect())" in main
    assert "POST_TRIM_FASTQC(ch_trimmed_reads)" in main
    assert "BUILD_SALMON_INDEX(" in main
    assert "SALMON_STRANDEDNESS_INFERENCE(ch_trimmed_reads, reference_salmon_index)" in main
