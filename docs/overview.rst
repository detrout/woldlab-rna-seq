Introduction
============

This is our RNA-Seq processing pipeline based on the work of the
ENCODE3 RNA-Seq working group. The official pipeline is designed to
run on DNA-Nexus, and is available at `long-rna-seq-pipeline`_.

We already had a compute cluster and running `HT-Condor`_ and so I
ported logic and arguments of the official pipeline to use the
HT-Condor dagman workflow manager.

We have primarily used this software to process Illumina style short
reads ranging from 30 to 150 bases in either single or paired end
formats. Our fragment sizes are greater than 200 base pairs and
typically range between 300 to 500 base pairs.

The name "long" is a historical artifact of ENCODE3, where the pipeline was
long when compared to a different protocol designed to work with sub
200 base pair fragments.

The default pipeline provided by this package is not appropriate for
processing full length RNA reads such as those produced by the Pacifc
Biosystems Sequel or the Oxford Nanopore MinION.

In addition, as time has progressed, this pipeline has a few
additional quality control reports added after ENCODE 3
development completed.

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
   
