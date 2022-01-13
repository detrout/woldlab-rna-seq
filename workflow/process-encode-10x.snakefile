#
# Copyright 2021 California Institute of Technology
#
# Licensed under a BSD like license that includes a
# non-endorsement clause.
#

DEFAULT_10X_CB_LENGTH = 16

configfile: "config.yaml"

include:
    "rules/star_solo_config.smk"

# Set defaults for this version of the ENCODE scRNA-seq
# pipeline
config.setdefault(
    "alignment_step_run",
    "barbara-wold:starsolo-10x-alignment-step-run"
)
config.setdefault(
    "quantification_step_run",
    "barbara-wold:starsolo-10x-quantification-step-run",
)

rule ALL:
    input:
        list_default_encode_targets()

include:
    "rules/provide_star_index.smk"

include:
    "rules/fastqs.smk"

include:
    "rules/star_co_file.smk"

include:
    "rules/star_solo_10x.smk"

include:
    "rules/archive_solo_counts.smk"

include:
    "rules/prepare-qc.smk"

include:
    "rules/encode-submit.smk"
