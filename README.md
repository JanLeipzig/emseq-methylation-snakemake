# EM-seq methylation workflow (Snakemake)

Snakemake workflow for EM-seq data: Bismark alignment and methylation calling, followed by
an analysis of CpG methylation in genes, exons/introns, intergenic and flanking regions.

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
Copyright © 2026 Jan Engelhardt

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details. 


Please cite: Länger et al. 2025, J Exp Zool B 344(7): 454–469