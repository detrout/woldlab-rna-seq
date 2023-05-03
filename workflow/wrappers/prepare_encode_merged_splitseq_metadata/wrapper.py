import datetime
from pathlib import Path
import pandas

from encoded_client.metadata import (
    compute_alignment_alias,
    compute_subpool_matrix_alias,
    to_array_sheet,
)
from encoded_client.hashfile import make_md5sum


def generate_star_solo_merged_metadata(config, records):
    datestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_type_to_file_type = {
        "unfiltered sparse gene count matrix of unique reads": "tar",
        "unfiltered sparse gene count matrix of all reads": "tar",
    }

    alignment_alias = compute_alignment_alias(
        config["alias_prefix"], config["biosample"], datestamp
    )
    rows = []
    for output_type in records:
        filename = records[output_type]
        file_type = output_type_to_file_type[output_type]

        derived_from = []
        for library_id in config["library"]:
            derived_from.append(
                compute_subpool_matrix_alias(config, library_id, output_type, datestamp)
            )

        obj = {
            "uuid": None,
            "accession": None,
            "dataset": config["experiment_accession"],
            "file_format": file_type,
            "output_type": output_type,
            "assembly": config["assembly"],
            "genome_annotation": config["genome_annotation"],
            "derived_from": derived_from,
            "md5sum": make_md5sum(filename),
            "file_size": Path(filename).stat().st_size,
            "submitted_file_name": str(filename),
            "award": config["award"],
            "lab": config["lab"],
        }
        if file_type == "bam":
            obj["step_run"] = config["alignment_step_run"]
            obj["aliases"] = [alignment_alias]
        elif file_type == "tar":
            obj["step_run"] = config["quantification_step_run"]
        else:
            print("Unknown file type {}".format(file_type))
        rows.append(obj)

    return rows


def generate_star_solo_merged_sheet(config, records):
    config = config.copy()

    rows = generate_star_solo_merged_metadata(config, records)

    sheet = pandas.DataFrame(rows)

    for array in ["aliases", "derived_from"]:
        if array in sheet.columns:
            sheet[array] = sheet[array].apply(to_array_sheet)

    sheet = sheet.rename(
        {
            "aliases": "aliases:array",
            "derived_from": "derived_from:array",
            "file_size": "file_size:integer",
        },
        axis=1,
    )

    return sheet


records = {
    "unfiltered sparse gene count matrix of unique reads": snakemake.input.gene_unique_raw,
    "unfiltered sparse gene count matrix of all reads": snakemake.input.gene_multi_raw,
}


metadata = generate_star_solo_merged_sheet(snakemake.config, records)
metadata = pandas.DataFrame(metadata)
metadata.to_csv(snakemake.output[0], index=False)
