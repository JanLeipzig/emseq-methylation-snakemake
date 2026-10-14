# Stage 2: methylation analysis of the Bismark output
#   per-sample CpG tables -> combined table -> overlap with genes / exons+introns /
#   intergenic regions / flanking regions -> summary report
import os

require_keys("analysis", ["analysis_dir", "exp_name", "genome_fasta",
                          "annotation.genes_gff", "annotation.exons_gff"])

OUT          = config["analysis_dir"]
EXP          = config["exp_name"]
GENES_GFF    = config["annotation"]["genes_gff"]
EXONS_GFF    = config["annotation"]["exons_gff"]
GENOME_FASTA = config["genome_fasta"]

SPIKE_CONTIGS = config.get("spike_contigs", ["Lambda_NEB", "pUC19"])
COV_CUTOFF    = config.get("coverage_cutoff", 5)     # CpGs with coverage > cutoff are kept
FLANK_SIZE    = config.get("flank_size", 1000)
SORT_THREADS  = config.get("sort_threads", 4)

SCRIPT = {name: os.path.join(SCRIPTS_DIR, fname) for name, fname in {
    "combine_cpgs":     "combineCpGs_v0.2.py",
    "gene_methrate":    "genes_methrate.py",
    "exon_intron_list": "make_exon_intron-list_v0.3.py",
    "fasta_length":     "fastalength_v2.py",
    "report":           "combineReport.py",
}.items()}

# output of stage 1 (naming scheme STEM comes from common.smk)
BISMARK_COV = f"{STEM}.deduplicated.bismark.cov.gz"

CPG_ALL      = f"{OUT}/{{sample}}.AllCpG.lst"
CPG_SPIKES   = f"{OUT}/{{sample}}.spikeIns.AllCpG.lst"
CPG_COMBINED = f"{OUT}/Combined.AllCpG.lst"
CPG_FILTERED = f"{OUT}/Combined.AllCpG.short.ge{COV_CUTOFF}.bed"

GENES_VS_CPG  = f"{OUT}/{EXP}_Genes_vs_CpG.lst"
CPG_VS_IGR    = f"{OUT}/CpG_vs_IGR.lst"
GENE_METHRATE = f"{OUT}/{EXP}_Genes.Methrate.ge{COV_CUTOFF}.ge3.lst"

EXON_INTRON_GFF    = f"{OUT}/{EXP}.ExonIntron.sort.gff"
EXON_INTRON_VS_CPG = f"{OUT}/{EXP}_ExonIntron_vs_CpG.lst"

GENOME_SIZE = f"{OUT}/genome.size"
UP_GFF      = f"{OUT}/Upstream.sort.gff"
DOWN_GFF    = f"{OUT}/Downstream.sort.gff"
UP_VS_CPG   = f"{OUT}/Upstream_vs_CpG.lst"
DOWN_VS_CPG = f"{OUT}/Downstream_vs_CpG.lst"

REPORT = f"{OUT}/{EXP}.MethylReport.txt"

TARGETS += [REPORT, GENE_METHRATE]


# Spike-in contigs (lambda, pUC19) go to a separate file; the rest is sorted
rule split_cpg_files:
    input:
        BISMARK_COV,
    output:
        cpg_all=CPG_ALL,
        cpg_spikes=CPG_SPIKES,
    params:
        spike_regex="|".join(SPIKE_CONTIGS),
    threads: SORT_THREADS
    log: f"{LOG_DIR}/analysis/split_cpg_files_{{sample}}.log"
    shell:
        """
        exec 2> {log}
        zcat {input} \
            | awk -v spikes={output.cpg_spikes} -v pat="{params.spike_regex}" \
                  '$1 ~ pat {{print > spikes; next}} {{print}}' \
            | sort --parallel={threads} -k1,1 -k2,2n > {output.cpg_all}
        touch {output.cpg_spikes}
        """


rule combine_cpg_table:
    input:
        script=SCRIPT["combine_cpgs"],
        cpgs=expand(CPG_ALL, sample=SAMPLES),
    output:
        CPG_COMBINED,
    conda: "../envs/python.yaml"
    log: f"{LOG_DIR}/analysis/combine_cpg_table.log"
    shell:
        "python {input.script} {input.cpgs} > {output} 2> {log}"


rule filter_cpgs:
    input:
        CPG_COMBINED,
    output:
        CPG_FILTERED,
    params:
        cutoff=COV_CUTOFF,
    log: f"{LOG_DIR}/analysis/filter_cpgs.log"
    shell:
        """
        exec 2> {log}
        awk -v OFS='\\t' -v c={params.cutoff} '!/^#/ && $3 > c {{print $1,$2,$3,$4}}' {input} \
            | sort -k1,1 -k2,2n > {output}
        """


# --- genes and intergenic regions -------------------------------------------
rule overlap_genes:
    input:
        genes=GENES_GFF,
        cpgs=CPG_FILTERED,
    output:
        genes_vs_cpg=GENES_VS_CPG,
        igr=CPG_VS_IGR,
    conda: "../envs/bedtools.yaml"
    log: f"{LOG_DIR}/analysis/overlap_genes.log"
    shell:
        """
        exec 2> {log}
        bedtools intersect -wao -a {input.genes} -b {input.cpgs} > {output.genes_vs_cpg}
        bedtools intersect -v -a {input.cpgs} -b {input.genes} > {output.igr}
        """


rule gene_methylation:
    input:
        script=SCRIPT["gene_methrate"],
        genes_vs_cpg=GENES_VS_CPG,
    output:
        GENE_METHRATE,
    conda: "../envs/python.yaml"
    log: f"{LOG_DIR}/analysis/gene_methylation.log"
    shell:
        "python {input.script} {input.genes_vs_cpg} > {output} 2> {log}"


# --- exons and introns ----------------------------------------------------------
rule make_exon_intron_list:
    input:
        script=SCRIPT["exon_intron_list"],
        exons=EXONS_GFF,
    output:
        EXON_INTRON_GFF,
    conda: "../envs/python.yaml"
    log: f"{LOG_DIR}/analysis/make_exon_intron_list.log"
    shell:
        """
        exec 2> {log}
        python {input.script} {input.exons} | sort -k1,1 -k4,4n > {output}
        """


rule overlap_exon_intron:
    input:
        features=EXON_INTRON_GFF,
        cpgs=CPG_FILTERED,
    output:
        EXON_INTRON_VS_CPG,
    conda: "../envs/bedtools.yaml"
    log: f"{LOG_DIR}/analysis/overlap_exon_intron.log"
    shell:
        "bedtools intersect -wao -a {input.features} -b {input.cpgs} > {output} 2> {log}"


# --- flanking regions (up-/downstream of genes) ----------------------------------
rule genome_size:
    input:
        script=SCRIPT["fasta_length"],
        genome=GENOME_FASTA,
    output:
        GENOME_SIZE,
    conda: "../envs/python.yaml"
    log: f"{LOG_DIR}/analysis/genome_size.log"
    shell:
        "python {input.script} {input.genome} > {output} 2> {log}"


rule flanking_regions:
    input:
        genes=GENES_GFF,
        genome_size=GENOME_SIZE,
    output:
        up=UP_GFF,
        down=DOWN_GFF,
    params:
        flank=FLANK_SIZE,
    conda: "../envs/bedtools.yaml"
    log: f"{LOG_DIR}/analysis/flanking_regions.log"
    shell:
        """
        exec 2> {log}
        bedtools flank -i {input.genes} -g {input.genome_size} -l {params.flank} -r 0 -s \
            | sort -k1,1 -k4,4n > {output.up}
        bedtools flank -i {input.genes} -g {input.genome_size} -l 0 -r {params.flank} -s \
            | sort -k1,1 -k4,4n > {output.down}
        """


rule overlap_flanking:
    input:
        up=UP_GFF,
        down=DOWN_GFF,
        cpgs=CPG_FILTERED,
    output:
        up=UP_VS_CPG,
        down=DOWN_VS_CPG,
    conda: "../envs/bedtools.yaml"
    log: f"{LOG_DIR}/analysis/overlap_flanking.log"
    shell:
        """
        exec 2> {log}
        bedtools intersect -wao -a {input.up} -b {input.cpgs} > {output.up}
        bedtools intersect -wao -a {input.down} -b {input.cpgs} > {output.down}
        """


# --- summary --------------------------------------------------------------------
rule write_report:
    input:
        script=SCRIPT["report"],
        exon_intron=EXON_INTRON_VS_CPG,
        genes=GENES_VS_CPG,
        genome=CPG_FILTERED,
        igr=CPG_VS_IGR,
        up=UP_VS_CPG,
        down=DOWN_VS_CPG,
    output:
        REPORT,
    conda: "../envs/python.yaml"
    log: f"{LOG_DIR}/analysis/write_report.log"
    shell:
        """
        python {input.script} {input.exon_intron} {input.genes} {input.genome} \
            {input.igr} {input.up} {input.down} > {output} 2> {log}
        """


# Open questions from the original file:
#  - Do exons and introns need separate handling?
#  - "At least 3 CpGs per feature" rule: apply also to the simple statistics?
#  - Compare with the Lde paper
#  - Percent methylation: multiply by 100?
