Introduction
============

This is our RNA-Seq processing pipeline based on the work of the
ENCODE4 RNA-Seq working group. The official pipeline was designed to
run on DNA-Nexus, and is available at `long-rna-seq-pipeline`_. There
is a newer pipeline designed to run with `cromwell`_ on several common
environments available at `rna-seq pipeline`_.

We already had a compute cluster and running `HT-Condor`_ and so I
ported logic and arguments of the official pipeline to use the
HT-Condor's `DAGMan`_ workflow language. We believe that DAGMan offers some
advantages over what is possible with `cromwell`_'s WDL workflow
language.

A later round of experiments led to trying the `snakemake`_ workflow
engine. Those processing pipelines are found in the workflows
directory.

These pipelines are primarily intended to process Illumina style short
reads ranging from 30 to 150 bases in either single or paired end
formats. Our fragment sizes are greater than 200 base pairs and
typically range between 300 to 500 base pairs.

The default pipeline provided by this package is not appropriate for
processing full length RNA reads such as those produced by the Pacifc
Biosystems Sequel or the Oxford Nanopore MinION.

In addition, as time has progressed, we have added a few additional
quality control metrics to this pipeline after the ENCODE 4
development was completed.

Workflow
--------

.. blockdiag::
   :caption: RNA-Seq workflow
   :align: center
   :desctable:

   blockdiag workflow {
     align-STAR -> sort-samtools;
     align-STAR -> index-samtools;
     align-STAR -> bedgraph-star;
     index-samtools -> samstats;
     index-samtools -> distribution;
     sort-samtools -> RSEM;
     bedgraph-star -> coverage;
     bedgraph-star -> bedgraph2bigwig;

     align-STAR [ description = "STAR genomic and transcriptome alignment"];
     sort-samtools [ description = "sort Transcriptome BAM file"];
     index-samtools [ description = "Index Genomic BAM file"];
     bedgraph-star [ description = "Generate bedgraph of reads"];
     samstats [ description = "Calculate simple read statistics"];
     distribution [ description = "Calcualate Exonic/Intronic/Genomic ratios"];
     RSEM [ description = "Generate RNA Genomic and Transcriptome counts"];
     coverage [ description = "Generate 5' to 3' unique transcript coverage"];
     bedgraph2bigwig [ description = "Generate read pileup bigwig"];
   }

.. _long-rna-seq-pipeline: https://github.com/ENCODE-DCC/long-rna-seq-pipeline
.. _HT-Condor: https://research.cs.wisc.edu/htcondor/
.. _rna-seq pipeline: https://github.com/ENCODE-DCC/rna-seq-pipeline   
.. _cromwell: https://software.broadinstitute.org/wdl/
.. _DAGMan: http://research.cs.wisc.edu/htcondor/manual/latest/DAGManApplications.html
.. _snakemake: https://snakemake.readthedocs.io/en/stable/
