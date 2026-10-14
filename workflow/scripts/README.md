# Helper scripts

Copy the five Python helper scripts used by `02_methylation_analysis.smk` into this folder
(or set `scripts_dir` in your config):

- `combineCpGs_v0.2.py`
- `genes_methrate.py`
- `make_exon_intron-list_v0.3.py`
- `fastalength_v2.py`
- `combineReport.py`

If their names change, adapt the `SCRIPT` table at the top of `workflow/rules/02_methylation_analysis.smk`.
Packages the scripts need go into `workflow/envs/python.yaml`.
