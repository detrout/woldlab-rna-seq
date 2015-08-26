===================
Long RNA Seq Condor
===================

Introduction
------------

The official `ENCODE Long RNA-Seq pipeline`_ is designed to run on
DNA-Nexus_. This is a conversion to use settings from the official
pipeline for a pipeline suitable for running through HT-Condor_ DagMan.

Installation
------------

One convenient thing on the ENCODE DCC portal provides is previously built
indexes and source files for STAR, RSEM, and Tophat.

For example experiments linked to from the ENCODE RNA-Seq
evaluation dataset ENCSR000AJW_ such as ENCSR000AED_ have a
visulaization showing how all the official ENCODE processed files were
created.  Between the two replicates there are cells indicating the
human female genome index files used, which in this case were TopHat ENCFF700TXC_
STAR ENCFF839KAI_ and RSEM ENCFF826ONU_

The mouse mm10 / M4 male genomic index files used were TopHat
ENCFF268IGM_, STAR ENCFF483PAE_, and RSEM ENCFF064YNQ_.

All the tar files for the supported software are also available from
the following links (with additional metadata describing which
versions of the GENCODE annotations were used). 

ENCSR641UDW_: Tophat indexed reference files with ERCC spike-ins
ENCSR314WMD_: STAR indexed reference files with ERCC spike-ins
ENCSR219BJA_: RSEM indexed reference files with ERCC spike-ins

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





