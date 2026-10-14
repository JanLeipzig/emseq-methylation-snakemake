# Stage 1: trimming -> Bismark mapping -> deduplication -> methylation calls (+ MultiQC)
import os

require_keys("alignment", ["genome_dir", "raw_dir", "raw_r1_name", "raw_r2_name"])

GENOME_DIR = config["genome_dir"]
RAW_R1 = os.path.join(config["raw_dir"], config["raw_r1_name"])   # keeps {sample} wildcard
RAW_R2 = os.path.join(config["raw_dir"], config["raw_r2_name"])

TRIM_DIR = config.get("trim_dir", "results/trimmed")
QC_DIR   = config.get("qc_dir", "results/qc")

_trim = config.get("trimming", {})
TRIM_CLIP_5P = _trim.get("clip_5p", {"R1": 10, "R2": 10})
TRIM_CLIP_3P = _trim.get("clip_3p", {"R1": 10, "R2": 10})
_bm = config.get("bismark", {})
BISMARK_SCORE_MIN = _bm.get("score_min", "L,0,-0.6")
EXTRACT_EXTRA     = _bm.get("extractor_extra", "--scaffolds")
EXTRACT_BUFFER    = _bm.get("extractor_buffer", "10G")

KEEP_TRIMMED_FASTQ = config.get("keep_trimmed_fastq", False)
RUN_MULTIQC        = config.get("run_multiqc", True)

TRIM_THREADS    = config.get("trim_threads", 4)
MAP_THREADS     = config.get("map_threads", 8)
DEDUP_THREADS   = config.get("dedup_threads", 1)
EXTRACT_THREADS = config.get("extract_threads", 2)

MAP_REPORT   = f"{SAMPLE_OUT}/{{sample}}_val_1_bismark_bt2_PE_report.txt"
TRIM_R1      = f"{TRIM_DIR}/{{sample}}_val_1.fq.gz"
TRIM_R2      = f"{TRIM_DIR}/{{sample}}_val_2.fq.gz"
MULTIQC_HTML = f"{QC_DIR}/multiqc_report.html"


def maybe_temp(path):
    return path if KEEP_TRIMMED_FASTQ else temp(path)


# Early, readable warnings instead of cryptic errors halfway through the run
if not os.path.isdir(os.path.join(GENOME_DIR, "Bisulfite_Genome")):
    logger.warning(f"No Bisulfite_Genome folder in {GENOME_DIR} - was bismark_genome_preparation run?")
for _s in SAMPLES:
    for _f in (RAW_R1, RAW_R2):
        if not os.path.exists(_f.format(sample=_s)):
            logger.warning(f"Missing raw file: {_f.format(sample=_s)}")

TARGETS += expand(f"{STEM}.deduplicated.bedGraph.gz", sample=SAMPLES)
if RUN_MULTIQC:
    TARGETS.append(MULTIQC_HTML)


rule trimming:
    input:
        r1=RAW_R1,
        r2=RAW_R2,
    output:
        r1=maybe_temp(TRIM_R1),
        r2=maybe_temp(TRIM_R2),
    params:
        clip_r1=TRIM_CLIP_5P["R1"], clip_r2=TRIM_CLIP_5P["R2"],
        clip3_r1=TRIM_CLIP_3P["R1"], clip3_r2=TRIM_CLIP_3P["R2"],
        outdir=TRIM_DIR,
    threads: TRIM_THREADS
    conda: "../envs/trimming.yaml"
    log: f"{LOG_DIR}/trimming/{{sample}}.log"
    shell:
        """
        trim_galore --paired --fastqc -j {threads} --basename {wildcards.sample} \
            --clip_R1 {params.clip_r1} --clip_R2 {params.clip_r2} \
            --three_prime_clip_R1 {params.clip3_r1} --three_prime_clip_R2 {params.clip3_r2} \
            --output_dir {params.outdir} {input.r1} {input.r2} &> {log}
        """


# Bismark's --parallel N starts N instances, each needing ~4 cores (2 bowtie2 processes
# with 2 threads each, plus helpers). So `threads: 8` means `--parallel 2`.
# The custom resource `mapping_slots` limits how many mapping jobs run simultaneously
# (set via --resources mapping_slots=N or in the profile).
rule bismark_mapping:
    input:
        genome=GENOME_DIR,
        r1=TRIM_R1,
        r2=TRIM_R2,
    output:
        bam=f"{STEM}.bam",
        report=MAP_REPORT,
    params:
        score=BISMARK_SCORE_MIN,
        instances=lambda wildcards, threads: max(1, threads // 4),
        outdir=SAMPLE_OUT,
    threads: MAP_THREADS
    resources: mapping_slots=1
    conda: "../envs/bismark.yaml"
    log: f"{LOG_DIR}/bismark_mapping/{{sample}}.log"
    shell:
        """
        bismark --score_min {params.score} --parallel {params.instances} \
            --genome {input.genome} -1 {input.r1} -2 {input.r2} \
            --output_dir {params.outdir} &> {log}
        """


rule bismark_deduplication:
    input:
        bam=f"{STEM}.bam",
    output:
        bam=f"{STEM}.deduplicated.bam",
        report=f"{STEM}.deduplication_report.txt",
    params:
        outdir=SAMPLE_OUT,
    threads: DEDUP_THREADS
    conda: "../envs/bismark.yaml"
    log: f"{LOG_DIR}/bismark_deduplication/{{sample}}.log"
    shell:
        "deduplicate_bismark -p {input.bam} --output_dir {params.outdir} &> {log}"


rule bismark_methylation_extraction:
    input:
        bam=f"{STEM}.deduplicated.bam",
    output:
        # CpG/CHG/CHH x OT/OB context files
        contexts=expand(
            f"{SAMPLE_OUT}/{{ctx}}_{{strand}}_{{sample}}_val_1_bismark_bt2_pe.deduplicated.txt.gz",
            ctx=["CpG", "CHG", "CHH"], strand=["OT", "OB"], allow_missing=True,
        ),
        splitting=f"{STEM}.deduplicated_splitting_report.txt",
        mbias=f"{STEM}.deduplicated.M-bias.txt",
        bedgraph=f"{STEM}.deduplicated.bedGraph.gz",
        cov=f"{STEM}.deduplicated.bismark.cov.gz",
    params:
        outdir=SAMPLE_OUT,
        buffer=EXTRACT_BUFFER,
        extra=EXTRACT_EXTRA,
    threads: EXTRACT_THREADS
    conda: "../envs/bismark.yaml"
    log: f"{LOG_DIR}/bismark_methylation_extraction/{{sample}}.log"
    shell:
        """
        bismark_methylation_extractor -p {params.extra} --buffer_size {params.buffer} \
            --parallel {threads} --bedGraph --gzip {input.bam} \
            --output_dir {params.outdir} &> {log}
        """


rule multiqc:
    input:
        expand(MAP_REPORT, sample=SAMPLES),
        expand(f"{STEM}.deduplication_report.txt", sample=SAMPLES),
        expand(f"{STEM}.deduplicated_splitting_report.txt", sample=SAMPLES),
        expand(f"{STEM}.deduplicated.M-bias.txt", sample=SAMPLES),
    output:
        MULTIQC_HTML,
    params:
        outdir=QC_DIR,
    conda: "../envs/multiqc.yaml"
    log: f"{LOG_DIR}/multiqc.log"
    shell:
        "multiqc --force --outdir {params.outdir} {input} &> {log}"
