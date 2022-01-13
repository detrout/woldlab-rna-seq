
===================
Long RNA Seq Condor
===================

|ci-test|

Table of Contents
=================

* `Introduction`_
* `Installation`_

  * `Required Dependencies`_

* `Obtaining Indexes`_
* `Configuration`_
* `Processing`_

  * `Running`_

Introduction
============

This started as a port of the ENCODE3 bulk RNA-seq pipeline for processing
Illumina short read sequences from DNA nexus to `HT-Condor`_'s dagman
language. as time progressed I've also added in some of our (Wold Lab)
quality-control steps after the core RNA-Seq STAR/rsem runs.

More recently I added new pipelines to run STAR on 10x style or
split-seq style experiments.

Installation
============

Required Dependencies
---------------------

* HT-Condor
* python3
* numpy
* pandas
* scipy.stats
* jinja2
* bokeh
* GeorgisScripts_

  * gene_coverage_wig_gtf.py
  * SAM_reads_in_genes3_BAM.py
  * SAMstats.py

Obtaining Indexes
-----------------

One convenient thing on the ENCODE DCC portal provides is previously built
indexes and source files for STAR, RSEM, and Tophat.

For example experiments linked to from the ENCODE RNA-Seq evaluation
dataset ENCSR000AJW_ such as ENCSR000AED_ have a visulaization showing
how all the official ENCODE processed files were created.  Between the
two replicates there are cells indicating the human female genome
index files used, which in this case were TopHat ENCFF700TXC_ STAR
ENCFF839KAI_ and RSEM ENCFF826ONU_

(In the following the links are are direct links to download the
files, instead of viewing the DCC metadata about the tar file)

The mouse mm10 / M4 male genomic index files used were:

============  ==========
Download      Index type
============  ==========
ENCFF483PAE_  STAR
ENCFF064YNQ_  RSEM
ENCFF268IGM_  TopHat
============  ==========

Late in ENCODE 3 we swiched to GRCh38 / V24 for the human reference
genome and GENCODE annotation, The index files are.

============  ==========
Download      Index type
============  ==========
ENCFF742NER_  STAR
ENCFF940AZB_  RSEM
ENCFF751OWT_  Tophat
============  ==========

I included the TopHat link although I'm not sure anyone bothered with
TopHat at this point in the project.

Earlier in ENCODE 3 we were using human hg19 reference genome with
GENCODE V19. The male genomic index files used were:

============  ==========
Download      Index type
============  ==========
ENCFF746IZL_  STAR
ENCFF826ONU_  RSEM
ENCFF515BAC_  TopHat
============  ==========

All the tar files containing the various indexes genome version,
annotation version, and chromosomal sex for the supported software are
also available from the following links.

============= ==================================================
Accession     Description
============= ==================================================
ENCSR641UDW_: Tophat indexed reference files with ERCC spike-ins
ENCSR314WMD_: STAR indexed reference files with ERCC spike-ins
ENCSR219BJA_: RSEM indexed reference files with ERCC spike-ins
============= ==================================================

Even though we are currently just running STAR/RSEM the Tophat indexes
include the prebuilt GFF file that contains all of tRNA and spikes
that we need for Wold Lab QC software.

You can also download the components for the GFF files from `GENCODE
reference files used in ENCODE pipelines`_. Though they might still
need to be processed to add the tRNA's and spike-ins commonly used.

To prepare an index directory:

* Pick a set of TopHat, STAR, and RSEM indexes one needs.
* Extract all three ENCFF files into the default out/ directory
* Rename the out/ directory to a "Genome Triple" directory name:
  Genome Build name - Annotation version - Genome sex. For example:
  mm10-M4-male or hg19-V19-male

*NOTE*: I discovered that the human and mouse GFF files are named
differently. The mouse GFF from the ENCFF268IGM_: file was named
gencode.vM4-tRNAs-ERCC.gff while the human GFF file from ENCFF515BAC_
was named annotation.gff. To work around this I needed for the
hg19-V19-male directory I needed to do the following in the
hg19-V19-male directory: ``ln -s annotation.gff
gencode.vV19-tRNAs-ERCC.gff``

Configuration
-------------

Currently this pipeline code needs a bit of help to find the
condor scripts and the genome indexes.

You need to create a configuration file in your home directory
called .htsworkflow.ini with the contents::

  [analysis]
  condor_script_dir=<path to long-rna-seq-condor repository>
  genome_dir=<path to where the root of where the indexes are extracted>
  star_dir=<directory containing STAR executable>
  rsem_dir=<directory containing rsem-calculate-expression>
  georgi_dir=<directory containing "GeorgiScripts" python files>
  ucsc_tools_dir=<directory containing at least bedGraphToBigWig from ucsc tools>

Processing
==========

Running
-------

First, the processing scripts need two files to be created, an
experiments and libraries definition tables.

The experiment file contains:

experiment
  Experiment Name

replicates
  Comma seperated list of library ids.

The more library table contains

library_id
  Library ID (used by experiment file to group related libraries)

genome
  Genome version string. e.g. mm10, hg19

annotation
  Annotation version string. e.g. M4 or V19

sex
  Sex of genome. male, female, uknown

analysis_dir
  Directory relative to the library table where the analysis files
  should be generated. It must exist before the script is run.

reference_prefix
  (Optional) you can use this to override the reference prefix for bedgraph
  bigwig file generation. It defaults to 'chr' but you might want all
  of the references '-' or you may have a genome that uses 'Scaffold'
  for as its referece prefix.

read_1
  Comma seperated list of unix filename globs specifying where to find
  the first read fastq files.

read_2
  (optional) Comma seperated list of unix filename glbos specifying where
  to find the second read (mate pairs) fastq files.

Second, after the definition files are constructed you need to create
the analysis directories. You can do that with this command. You need to
change the `-f 2` to be whatever column you used for analysis_dir. I usually put
analysis_dir as the second column, so I used `-f 2.`::

  tail -n +2 library.tsv  | cut -f 2 | xargs mkdir

Next you can generate the DagMan script to generate the result files with::

  make_dag <list of library.tsv files> > <filename>.dagman
  condor_submit_dag <filename>.dagman

**TODO** Currently the QC summary statistics and report generation are
not integrated into the condor pipeline and need to be run manually

**NOTE** If condor_submit_dag fails it will generate a rescue file
<filename>.dagman.rescue<number>. After investigating the log files
to find the cause of the error you can do::

  condor_submit_dag -autorescue 1 <filename>.dagman

to try to contine.

Fourth, to generate the HDF5 files containing the various pairwise
correlation scores one needs to do::

   madqc -l <library.tsv> -e <output_experiment_name> <list_of_library ids>

**NOTE** Yes. That is currently annoying, the ``make_dag.py`` is
supposed to generate the commands, but it doesn't yet.

Fifth, after all of the experiment correlation scores are generated one can
construct a summary report with::

  qcreport -l <library.tsv> -e <experiment.tsv> > <html filename>

**TODO** Implement a way to specify where the Bokeh JavaScript and CSS is.

Lastly, you probably should delete any bam and bedgraph (.bg) files
you are not planning on using.

Processing Phases
-----------------

Steps for our processing pipeline:

* align-star-se
* sort-samtools
* quant-rsem
* index-samtools
* qc-samstats
* bedgraph-star
* qc-distribution
* qc-coverage
* bedgraph2bigwig

.. references

.. _DNA-Nexus: https://www.dnanexus.com/
.. _HT-Condor: http://research.cs.wisc.edu/htcondor/
.. _ENCODE Long RNA-Seq pipeline: https://github.com/ENCODE-DCC/long-rna-seq-pipeline
.. _ENCSR000AJW: https://www.encodeproject.org/datasets/ENCSR000AJW/
.. _ENCSR000AED: https://www.encodeproject.org/experiments/ENCSR000AED/
.. _ENCSR219BJA: https://www.encodeproject.org/datasets/ENCSR219BJA/
.. _ENCSR641UDW: https://www.encodeproject.org/datasets/ENCSR641UDW/
.. _ENCSR314WMD: https://www.encodeproject.org/datasets/ENCSR314WMD/
.. _ENCSR219BJA: https://www.encodeproject.org/datasets/ENCSR219BJA/
.. _ENCFF268IGM: https://www.encodeproject.org/files/ENCFF268IGM/@@download/ENCFF268IGM.tar.gz
.. _ENCFF483PAE: https://www.encodeproject.org/files/ENCFF483PAE/@@download/ENCFF483PAE.tar.gz
.. _ENCFF064YNQ: https://www.encodeproject.org/files/ENCFF064YNQ/@@download/ENCFF064YNQ.tar.gz
.. _ENCFF700TXC: https://www.encodeproject.org/files/ENCFF700TXC/@@download/ENCFF700TXC.tar.gz
.. _ENCFF839KAI: https://www.encodeproject.org/files/ENCFF839KAI/@@download/ENCFF839KAI.tar.gz
.. _ENCFF515BAC: https://www.encodeproject.org/files/ENCFF515BAC/@@download/ENCFF515BAC.tar.gz
.. _ENCFF746IZL: https://www.encodeproject.org/files/ENCFF746IZL/@@download/ENCFF746IZL.tar.gz
.. _ENCFF826ONU: https://www.encodeproject.org/files/ENCFF826ONU/@@download/ENCFF826ONU.tar.gz
.. _ENCFF742NER: https://www.encodeproject.org/files/ENCFF742NER/@@download/ENCFF742NER.tar.gz
.. _ENCFF940AZB: https://www.encodeproject.org/files/ENCFF940AZB/@@download/ENCFF940AZB.tar.gz
.. _ENCFF751OWT: https://www.encodeproject.org/files/ENCFF751OWT/@@download/ENCFF751OWT.tar.gz
.. _GeorgisScripts: https://github.com/detrout/GeorgiScripts
.. _GENCODE reference files used in ENCODE pipelines: https://www.encodeproject.org/references/ENCSR884DHJ/
.. |ci-test| image:: https://github.com/detrout/woldlab-rna-seq/actions/workflows/ci-test.yml/badge.svg
