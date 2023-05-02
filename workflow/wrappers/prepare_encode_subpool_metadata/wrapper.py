import pandas

from encoded_client.metadata import generate_star_solo_subpool_sheet


records = {
    "alignments": snakemake.input.bam,
    "unfiltered sparse gene count matrix of unique reads": snakemake.input.gene_unique_raw,
    "unfiltered sparse gene count matrix of all reads": snakemake.input.gene_multi_raw,
    "unfiltered sparse splice junction count matrix of unique reads": snakemake.input.sj_unique_raw,
}
library_id = snakemake.wildcards.library_id

metadata = generate_star_solo_subpool_sheet(
    snakemake.config, records, library_id=library_id
)
metadata = pandas.DataFrame(metadata)
metadata.to_csv(snakemake.output[0], index=False)
