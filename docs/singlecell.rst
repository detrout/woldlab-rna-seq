.. _reference.scRNA-seq:

Single Cell RNA-seq 
===================

There are two different single cell RNA-seq pipelines available.  One
is designed for 10x Genomic kits and should work with Chromium version
1 through 3.1 kits, as well as the new combined RNA and multiome
kit. The other is designed to work the Parse Biosciences Split-seq
kits.

The primary differences in arguments to STAR to handle the different
methods revolve around differences in barcode handling. The Parse
Split-seq method has a more complex process to decode the barcode, and
additionally has both a poly-A and random priming barcode for each
cell. The split-seq processing pipeline logically combines the two
actual barcodes into a single merged barcode summing the two reported
values into a single gene score.

The snakemake entry points are
`workflows/process-encode-10x.snakemake` for the 10x chromium pipeline
and `workflows/process-encode-splitseq.snakemake` for the split-seq
pipeline.

They are both driven by similar configuration files and currently
designed to process ENCODE submitted data.

.. code-block:: yaml

    ## Read 1 should contain the real sequence
    read1:
       # list of ENCODE fastq accession IDs
       - ENCFF150FBF
       - ENCFF385IAW
    ## Read 2 should contain the cell barcode and UMI
    read2:
       - ENCFF351VBS
       - ENCFF503CCI
    ## for single nucleus experiments, you will likely want to include
    ## the intronic reads
    ## include_intron: (boolean True or False, defaults to False)
    ## You need to specify the strandedness, for the majority
    ## of protocols including it is likely to be Forward
    stranded: "Forward"
    ## For posting the results back to the ENCODE portal, you need to
    ## specify the ENCODE experiment ID
    experiment_accession: "ENCSR724KET"
    ## store the experiment description for human consumption
    description: "snRNA on human adrenal gland."
    ## Library id is needed to identify replicates.
    #library_accession: "ENCLB002DZK"
    ## Due to the limits of singularity, we need to define a target
    ## directory for the genome indexes to be available at that is
    ## below the analysis directory. The default is "genome"
    ## If should be possible to reuse an extracted version of the
    ## genome indexes by passing
    ## --singularity-args "-B <index dir>:$(pwd)/genome:ro"
    ## to singularity.
    # genome_dir: "genome"
    ## If you provide an accession for the genome index the code will 
    ## attempt to download the genome index from that location
    #genome_accession: "TSTFF984371"
    ## If you only provide a genome_index_url, the code will attempt to use the filename
    ## portion to determine the accession on the portal.
    #genome_index_url: "https://woldlab.caltech.edu/~diane/genome/GRCh38-V29-male-2.7.8a.tar.gz"
    ## If the genome_accession is set the code will pull the assembly and genome_annotation
    ## from the ENCODE portal
    # assembly: "GRCh38"
    # genome_annotation: "V29"
    ## Given the accession it will determine the inclusion_list_url
    ##inclusion_accession: "737K-arc-v1(GEX)"
    ## if there is no inclusion_accession, it will guess the accession from the inclusion_list_url
    #inclusion_list_url: "https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz"
    ## The pipeline now defaults to the container used for my analysis (The version of STAR after 2.7.9a)
    # You'll probably get better performance if you put a copy somewhere local
    #star_container: "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
    ## This is a container with python and scanpy in it. It also defaults to my container.
    # You'll probably get better performance if you put a copy somewhere local
    #scanpy_container: "https://woldlab.caltech.edu/~diane/containers/bullseye-scanpy-1.8.2.sif"
    #
    ## These are resource limit guesses based on my current
    ## experiments and might need to be increased for larger experiments.
    #mem_mb: 65536
    #disk_mb: 51200

This configuration files were generated by a notebook from reading the
encode portal for unprocessed experiments.

https://github.com/ENCODE-AWG/encode-202006-jamboree-detrout-rna-sc-pipeline/blob/master/single-cell-to-process.ipynb
