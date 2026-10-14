# EM-seq methylation workflow (Snakemake)

Snakemake workflow for EM-seq data: Bismark alignment and methylation calling, followed by
an analysis of CpG methylation in genes, exons/introns, intergenic and flanking regions.

> Working title - rename the repository as you like.

## Stages

| Stage | File | What it does |
|---|---|---|
| 1 `alignment` | `workflow/rules/01_alignment_and_calling.smk` | Trim Galore, Bismark (paired-end), deduplication, methylation extraction, MultiQC |
| 2 `analysis` | `workflow/rules/02_methylation_analysis.smk` | CpG tables, combined table, overlaps with genes / exons+introns / IGR / flanks (bedtools), report |

Stage 2 reads the `*.deduplicated.bismark.cov.gz` files of stage 1.

## Requirements

- Snakemake >= 8
- conda (or mamba); all tools are installed automatically from `workflow/envs/*.yaml`
- A Bismark-prepared genome (`bismark_genome_preparation`), a genome FASTA and gene/exon GFFs for stage 2
- The helper scripts in `workflow/scripts/` (see the README there)

## Quick start

```bash
cp config/config_template.yaml config/config_myexp.yaml   # edit paths, samples, parameters
snakemake --profile profiles/default --configfile config/config_myexp.yaml -n   # dry-run
snakemake --profile profiles/default --configfile config/config_myexp.yaml
```

Run from the repository root; Snakemake finds `workflow/Snakefile` itself.
A config file is always required (there is deliberately no default, so settings of different
experiments are never merged).

Useful overrides on the command line:

```bash
--config stages=alignment            # only stage 1 (or: stages=analysis)
--config groups=adult,embryo         # sample groups to process
--cores 16 --resources mapping_slots=3 --config map_threads=4
```

## Configuration

- **`config/config_<experiment>.yaml`** - everything specific to one experiment (data, samples,
  parameters). Documented in `config/config_template.yaml`.
- **`profiles/<name>/config.yaml`** - settings of the machine (cores, conda location, resource limits).
- Output goes to `results/` by default (folders can be changed in the config).

### Cores and memory for mapping

Bismark's `--parallel N` starts N instances of about 4 cores each. The workflow therefore runs
`map_threads // 4` instances per sample. The custom resource `mapping_slots` limits how many
samples are mapped at the same time.

## Layout

```
workflow/
  Snakefile                 # includes the stage files
  rules/                    # common.smk, 01_..., 02_...
  envs/                     # conda environments (one per tool group)
  scripts/                  # helper scripts for stage 2
config/                     # config template + experiment configs
profiles/                   # machine settings
```

## Citation / License

Add a `LICENSE` file and, if you publish, a `CITATION.cff` before making the repository public.
