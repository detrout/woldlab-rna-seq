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


def list_encode_10x_targets():
    default = [
        "Log.final.out",
        "Aligned.sortedByCoord.out.bam",
        "{}_Unique_filtered.tar.gz".format(get_gene_model()),
        "{}_EM_filtered.tar.gz".format(get_gene_model()),
        "{}_Unique_raw.tar.gz".format(get_gene_model()),
        "{}_EM_raw.tar.gz".format(get_gene_model()),
        "SJ_Unique_raw.tar.gz",
        "umi_per_cell.png",
        "qc_metric_violin.{}_EM_filtered.png".format(get_gene_model()),
        "pct_count_mt.{}_EM_filtered.png".format(get_gene_model()),
        "n_genes_by_counts.{}_EM_filtered.png".format(get_gene_model()),
        "metadata.{}.csv".format(get_submit_host(config)),
    ]
    if automatic_submission:
        default.extend([
            "posted.{}.csv".format(get_submit_host(config)),
            "Log.final.out.{}.qc-upload".format(get_submit_host(config)),
            SOLO_ROOT / get_gene_model() / "Summary.csv.{}.qc-upload".format(
                get_submit_host(config)
            ),
            "pct_count_mt.{}_EM_filtered.png.{}.qc-upload".format(
                get_gene_model(), get_submit_host(config)
            ),
        ])
    return default


rule ALL:
    input:
        list_encode_10x_targets()

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
    "rules/prepare_star_counts_qc.smk"

include:
    "rules/create_md5s.smk"

include:
    "rules/post-star-qc.smk"

include:
    "rules/post-star-solo-qc.smk"

include:
    "rules/encode_submit_10x.smk"
