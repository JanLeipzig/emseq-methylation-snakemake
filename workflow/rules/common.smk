# Shared by all stages: config checks, sample list, common folders and file naming.
import os
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger

if not config:
    raise WorkflowError("No config given. Run with: --configfile config/<your_config>.yaml")


def require_keys(stage, keys):
    """Stop with a readable message if required config keys are missing (dotted = nested)."""
    missing = []
    for key in keys:
        node = config
        for part in key.split("."):
            if not isinstance(node, dict) or part not in node:
                missing.append(key)
                break
            node = node[part]
    if missing:
        raise WorkflowError(f"[{stage}] missing required config key(s): {', '.join(missing)}")


# --- which stages to run ---------------------------------------------------------
STAGES = config.get("stages", ["alignment", "analysis"])
if isinstance(STAGES, str):                   # from the command line: --config stages=alignment
    STAGES = STAGES.split(",")
_bad = [s for s in STAGES if s not in ("alignment", "analysis")]
if _bad:
    raise WorkflowError(f"Unknown stage(s) {_bad}; choose from: alignment, analysis")

# --- samples ---------------------------------------------------------------------
require_keys("common", ["sample_groups"])
SAMPLE_GROUPS = config["sample_groups"]
_groups = config.get("groups", list(SAMPLE_GROUPS))
if isinstance(_groups, str):                  # from the command line: --config groups=a,b
    _groups = _groups.split(",")
_unknown = [g for g in _groups if g not in SAMPLE_GROUPS]
if _unknown:
    raise WorkflowError(f"Unknown group(s) {_unknown}; available: {list(SAMPLE_GROUPS)}")
SAMPLES = [str(s) for g in _groups for s in SAMPLE_GROUPS[g]]

# --- folders and naming shared between the stages ----------------------------------
BISMARK_DIR = config.get("bismark_dir", "results/bismark")
LOG_DIR     = config.get("log_dir", "results/logs")

SAMPLE_OUT = f"{BISMARK_DIR}/{{sample}}"
# Bismark names its files after the first trimmed read file; defined ONCE, used by both stages
STEM = f"{SAMPLE_OUT}/{{sample}}_val_1_bismark_bt2_pe"

# Helper scripts (default: workflow/scripts)
SCRIPTS_DIR = config.get("scripts_dir", os.path.join(workflow.basedir, "scripts"))

# final files requested by `rule all`; each stage adds its own
TARGETS = []

wildcard_constraints:
    sample="[^/]+"
