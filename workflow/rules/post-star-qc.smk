import logging
from encoded_client.encoded import ENCODED, DCCValidator

from woldrnaseq.models import (
    load_star_final_log,
)


def prepare_star_qc_metric(config, metric_of, filename):
    star_log = load_star_final_log(filename)
    star_quality_metric = {
        "Mapping speed, Million of reads per hour": star_log.loc[("", "Mapping speed, Million of reads per hour")],
        "Number of input reads": star_log.loc[("", "Number of input reads")],
        "Average input read length": star_log.loc[("", "Average input read length")],

        "Average mapped length": star_log.loc[("UNIQUE READS", "Average mapped length")],
        "Deletion average length": star_log.loc[("UNIQUE READS", "Deletion average length")],
        "Deletion rate per base": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Deletion rate per base")]),
        "Insertion average length": star_log.loc[("UNIQUE READS", "Insertion average length")],
        "Insertion rate per base": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Insertion rate per base")]),
        "Mismatch rate per base, %": "{:0.2f}%".format(star_log.loc[("UNIQUE READS", "Mismatch rate per base, %")]),
        "Number of splices: AT/AC": star_log.loc[("UNIQUE READS", "Number of splices: AT/AC")],
        "Number of splices: Annotated (sjdb)": star_log.loc[("UNIQUE READS", "Number of splices: Annotated (sjdb)")],
        "Number of splices: GC/AG": star_log.loc[("UNIQUE READS", "Number of splices: GC/AG")],
        "Number of splices: GT/AG": star_log.loc[("UNIQUE READS", "Number of splices: GT/AG")],
        "Number of splices: Non-canonical": star_log.loc[("UNIQUE READS", "Number of splices: Non-canonical")],
        "Number of splices: Total": star_log.loc[("UNIQUE READS", "Number of splices: Total")],
        "Uniquely mapped reads %": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Uniquely mapped reads %")]),
        "Uniquely mapped reads number": star_log.loc[("UNIQUE READS", "Uniquely mapped reads number")],

        "Number of reads mapped to multiple loci": star_log.loc[("MULTI-MAPPING READS", "Number of reads mapped to multiple loci")],
        "% of reads mapped to multiple loci": "{:0.2f}%".format(star_log.loc[("MULTI-MAPPING READS", "% of reads mapped to multiple loci")]),
        "Number of reads mapped to too many loci": star_log.loc[("MULTI-MAPPING READS", "Number of reads mapped to too many loci")],
        "% of reads mapped to too many loci": "{:0.2f}%".format(star_log.loc[("MULTI-MAPPING READS", "% of reads mapped to too many loci")]),

        "% of reads unmapped: too many mismatches": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: too many mismatches")]),
        "% of reads unmapped: other": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: other")]),
        "% of reads unmapped: too short": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: too short")]),

        "% of chimeric reads": "{:0.2f}%".format(star_log.loc[("CHIMERIC READS", "% of chimeric reads")]),
        "Number of chimeric reads": star_log.loc[("CHIMERIC READS", "Number of chimeric reads")],

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["alignment_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }
    return star_quality_metric


rule post_star_qc:
    input:
        log_final = "Log.final.out",
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
    output:
        "Log.final.out.{}.qc-upload".format(get_submit_host())
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "post_star_qc.log"
    run:
        logger = logging.getLogger("post_star_qc")
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.FileHandler(log[0]))
        logger.addHandler(logging.StreamHandler(sys.stderr))
        host = get_submit_host()
        server = ENCODED(host)
        uploaded = pandas.read_csv(input.metadata_posted)
        accession = uploaded[uploaded["file_format"] == "bam"]["accession"].to_list()
        qc = prepare_star_qc_metric(config, accession, input.log_final)
        try:
            validator = DCCValidator(server)
            validator.validate(qc, "star_quality_metric")
            results = server.post_json("/star_quality_metric/", qc)
            with open(output[0], "wt") as outstream:
                outstream.write(str(results))
        except Exception as e:
            logger.error(e)
            raise(e)


