Introduction
============

This is our RNA-Seq processing pipeline based on the work of the
ENCODE RNA-Seq working group. The official pipeline is designed to run
on DNA-Nexus, and is available at `long-rna-seq-pipeline`_

We however already had a compute cluster and running `HT-Condor`_ and so I
ported logic and arguments of the official pipeline to use the
HT-Condor dagman workflow manager.

In addition as time has progressed there are a few additional quality
control reports that we added for our own purposes.

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
   
