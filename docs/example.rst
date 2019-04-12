Example run
===========

Continuing from the :ref:`installation` instructions we will now
attempt an example run.


To test the installation I picked two of our smallest individual
single cell libraries to treat as if they were bulk samples.

.. code-block:: console

    mkdir ~/proj/rna-test
    cd ~/proj/rna-test
    wget https://www.encodeproject.org/files/ENCFF213IBI/@@download/ENCFF213IBI.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF863CYU/@@download/ENCFF863CYU.fastq.gz

    mkdir ENCFF213IBI
    mkdir ENCFF863CYU

Create your :ref:`howto.library.tsv` file. It needs to be a tab
delimited text file.

=========== ============ ====== ========== ====== =====================
library_id  analysis_dir genome annotation sex    read_1
=========== ============ ====== ========== ====== =====================
ENCFF213IBI ENCFF213IBI  mm10   M4         male   ENCFF213IBI.fastq.gz
ENCFF863CYU ENCFF863CYU  mm10   M4         male   ENCFF863CYU.fastq.gz
=========== ============ ====== ========== ====== =====================

.. code-block:: console

    make_dag -l library.tsv > run.dagman
    condor_submit_dag run.dagman


Next create your :ref:`howto.experiment.tsv` file. It also needs to be
a tab delimited text file.

=========== =======================
experiments replicates
=========== =======================
experiment  ENCFF213IBI,ENCFF863CYU
=========== =======================

With that in hand you can move on to generating the quantification
files and reports.

.. code-block:: console

    madqc -l library.tsv -e experiment.tsv
    make_rsem_csv -l library.tsv -e experiment.tsv -q TPM
    qc_report -l library.tsv -e experiment.tsv -o ~/public_html/report-NAME.html
    make_trackhub -l library.tsv -e experiment.tsv \
                  --hub test \
                  --short-name test \
                  --email me@example.edu \
                  --web-root http://example.org/~user/test/ \
                  --bigwig

The web-root variable defines the base URL that can be used to find
the trackhub files test.hub.txt, test.genome.txt, mm10/trackDb.txt,
and the associated bigWig files in the library analysis directories.

The name for the .hub.txt and .genome.txt file comes from the --hub
parameter. The short name is used to label the hub on the left side of
the genome browser. There is a long name parameter which can set a
description that goes above the track on the genome browser.

.. note::

   These are from an archived experiment. but they were replaced
   purely for the convience of the DCC. Instead of having many
   individual fastq libraries they wanted all the libraries to be
   merged into a single fastq file that will need demultiplexing.

   However for the purposes of testing it will be far easier to use
   small input files.
