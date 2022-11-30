Changelog
=========

Release ?????
-------------

Shorten the genome name from the genome, annotation, sex triple to
just a name, that is used as a directory name.

Extend make_trackhub to support multiple genomes in one run.

Add prototype bulk RNA-seq snakemake workflow in
`workflow/sprocess-encode-bulk-rna-star-rsem.snakefile`

Release 1.2.5
-------------

It was decided that splitseq should only submit the unfiltered files
and not STAR's empty drops filtered set. This led to a bunch of
refactoring to remove repeated code and get the rules to properly
fire for attaching the QC matrix to the unfiltered EM count matrix.

The count matrix merging code memory limits were set too low for a few
datasets and needed to be increased.

This is the version used to process single-cell RNA-seq 10x multiome
and splitseq experiments for the ENCODE4 2021 January freeze.

Release 1.2.4
-------------

Use released mex_gene_archive==0.2.1 and encoded-client==0.1.1 packages.

Try to catch cases where a file to be uploaded has already been
uploaded by md5sum. If it has been uploaded provide the uuid and
accession and skip doing a new upload.

Unfortunately, the pipeline isn't always binary reproducable, so
sometimes there will still be duplication.

Use encoded_client to download the genome_index_url in case it hasn't
been released yet.

Release 1.2.3
-------------

Implement the Snakemake split-seq pipeline launched from
workflow/process-encode-splitseq.snakefile. There was some extra
complexity added due to how split-seq needs collapsing two cell
barcodes into a single logical cell. Though the STAR matrices contain
the two barcodes, the real barcodes are merged into their logical
cells when building the mex_gene_archives.

The code to handle the merging is located in
woldrnaseq/splitseq_merger.py

Release 1.2.2
-------------

Additional work on being able to download indexes and inclusion lists
from the ENCODE DCC portal.

Support posting the count violin plot to the DCC.

Changed to PEP517/PEP518 project layout.

Fixed some issues in the unit tests for the condor based pipeline.

Release 1.2.1
-------------

The processing pipeline process-encode-10x.snakefile is handling
processing 10x v2 v3 and GEX+ATAC protocols pretty well, and can post
the results to the DCC, when the lab, award, genome index and
inclusion list accessions are available in the config.yaml file.

Release 1.2.0
-------------

Snakemake is looking to be a much more popular workflow language than
condor_dagman and I started implementing some workflows in snakemake
instead.

The most important of right now is a 10x pipeline using STAR Solo that
will hopefully be useful for processing the ENCODE 4 single cell
RNA-seq experiments.

"workflows/starsolo.snakemake" is a minimal viable product in that
given a configuration file with some ENCODE fastq accession IDs in it
will download them from the ENCODE portal, process them and leave
archives containing the results in our directory.

The next phase of development will to automate generating the
configuration files from ENCODE metadata and preparing to submit the
archived results to the portal


2021 Sept 2
-----------

Renamed project from long-rna-seq-condor to woldlab-rna-seq, because this
was always intended for illumina sized reads, it was only ever long
compared to some especially short ENCODE 3 fragments.

The new url is https://github.com/detrout/wold-rna-seq

And I dropped condor because I've been starting to port to snakemake.


Release 1.1.1
-------------

  * Support generating stranded bigwigs
  * Requires GeorgiScripts 1.1 (available from
    https://github.com/detrout/GeorgiScripts )
  * Generate coverage plots directly from bam files
  * Avoid leaving intermediate bedgraph files in the analysis
    directory.
  * Include the star genome index size in the estimate of STAR's
    memory footprint.
  * Add --merge-lanes to downloader.py to combine fastqs from
    multiple lanes into a single file.

Release 1.1
-----------

  * Requires bokeh ~ 0.12.15
  * Requires GeorgiScripts 1.0 (available from
    https://github.com/detrout/GeorgiScripts )
  * Improved sorting on compute hosts
  * Support more gff/gtf variants.
  * Add gene_coverage_detail plot
  * Convert makersemcsv to use --genome-dir instead of
    needing to provide path to a genome cache file
  * Adjust QC coverage script to to single gene model, with each gene
    normalized to its max expression, and requiring at least 1.0 in
    the STAR normalized bedgraph file
  * Support processing stranded libraries.
    Set 'stranded' attribute for a library in the library metadata file.
    If not present we default to to unstranded.
  * Estimate the disk space needed for a star alignment instead of
    using a hard coded 60G requirement.

Release 1.0
-----------

  * Requires bokeh ~ 0.12.9
  * Generate bam comments for STAR run
  * Supports running report & summary generator in a different
    directory than where the result files are found.
    (See --root argument to madqc and report)
  * Changes to Report

    * Report is heavily modified to now use bokeh plots
    * The genes detected plot was added to the standard report
    * The per library spike in scatter plot was removed
    * The per experiment scatter box plot was removed (at least
      temporarily)

Release 0.9
-----------

This introduces several backwards incompatible changes.

wold-rna-seq-condor now

  * has a setup.py script
  * requires samtools >= 1.3
  * components were moved into a package
  * Support paired end reads
    Libraries have a read_2 column to specify what fastqs are the second read
  * makersemcsv can output either genome or transcriptome summary files
    the filenames have changed slightly to avoid colliding
  * Added support for # as a comment charater in experiment and library files
  * Now require path to UCSC tools (add ucsc_tools_dir to .htsworkflow.ini)
  * long-rna-seq-condors version number is written to the dagman script
  * Uses condor scratch directories for STAR, RSEM, and sorting.
    This should lower NFS load, and leave fewer temp files around that
    need to be cleaned.

Release 0.3
-----------

Adjust madqc.py script to take different arguments.
madqc.py -e now specifies an experiment table like
some of the other scripts. To get the old early
behavior (or if you don't have an experiment table),
you can use -n or --experiment-name to specify
a name for an experiment along with a list of replicates.

Add makersemcsv.py script to read rsem files and
write out a csv file for some column of interest
in various libraries.

Release 0.2
-----------

Date: 2015 Dec 9

This version introduces three new required parameters
so it can be installed on someone elses compute cluster.

The previous version had a number of hard coded
paths in the condor scripts.

So now you'll need to define

  * star_dir
  * rsem_dir
  * georgi_dir

To define the paths where the pipline code expects to find
several pieces of software.

Release 0.1
-----------

Initial release. It works in my hands, and my coworkers who sits
on the other side of the room from me.
