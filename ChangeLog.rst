Changelog
=========

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
