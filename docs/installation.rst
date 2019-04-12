.. _installation:

Installation
============

I tested these instructions with Debian Stretch 9.8, but it
should work similiarly with Ubuntu 16.x or 18.x.

Assuming you don't have a Linux system available, first install the
base system according to the vendor's instructions.

If not installed, you will also need to install `HT Condor`_. On older
versions the package might be called condor instead. Installing HT
Condor properly can be fairly challenging, but to get a basic test
environment you should be able set up a single host "Personal Condor
Pool" using the instructions in :ref:`Personal Condor`.

Lets install the minimal dependencies for long-rna-seq-condor.

.. code-block:: console

    sudo apt install python3 python3-pip git samtools

First cd to whatever directory you'd like a repository to be checked
out to. For example

.. code-block:: console

     mkdir ~/proj
     cd ~/proj

Next we need to install the project

.. code-block:: console

    python3 -m pip install git+https://github.com/detrout/long-rna-seq-condor.git#egg=woldrnaseq
    git clone https://github.com/detrout/GeorgiScripts


add ~/.local/bin to your PATH. The best solution is to edit your shell
initialzation file `.profile`, but for now you can also do
`export PATH=~/.local/bin:$PATH`

Installing STAR and RSEM
------------------------

If you don't need to exactly match software versions you might be able to do:

.. code-block:: console

   apt install rna-star rsem

Otherwise see the installation instructions for `STAR`_ and `rsem`_ or ask your
system adminstrator where the software might be installed.


Installing UCSC tools
---------------------

Unfortunately UCSC has a license limiting to non-commercial use, so
isn't available to be installed with apt.

We need the bedSort and bedGraphToBigWig utilities, you may, of
course, want to download other utilities as well.

.. code-block:: console

   mkdir ~/proj/ucsc_tools
   cd ~/proj/ucsc_tools
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
   chmod a+x bedSort bedGraphToBigWig

Installing ENCODE DCC prebuilt indexes
--------------------------------------

For ENCODE mm10 M4 we can use these prebuilt indexes provided by
the DCC. The main advantage to this method is you can save several
hours of computer time by not having to build the indexes. The
downside is you are limited to ENCODE's choices and versions. You may
want to use different versions of the supporting software or you might
want to use different annotations than what ENCODE chose to do.

.. code-block:: console

   mkdir ~/proj/genome
   cd ~/proj/genome/
   # The STAR index
   wget https://www.encodeproject.org/files/ENCFF533JRE/@@download/ENCFF533JRE.tar.gz
   tar xavf ENCFF533JRE.tar.gz
   # The RSEM Index
   wget https://www.encodeproject.org/files/ENCFF717NFD/@@download/ENCFF717NFD.tar.gz
   tar xavf ENCFF717NFD.tar.gz
   mv out mm10-M4-male

I found these by looking at a recent RNA-seq experiment for the
species of interst and looked for the accession ids in the circles
that feed into the square program run steps.  (See the bottom of the
following screen shot.)

.. image:: _static/Screenshot_2019-04-08\ ENCSR968QHO\ –\ ENCODE.png

Unfortunately the ENCODE STAR and RSEM indexes lack a "GTF/GFF" file,
So we'll need to build one. For example using the RSEM input file
above we can look at its detail page `ENCFF717NFD`_. From there we can
see the input GTF and fasta files.  `ENCFF001RTP`_ is the accession ID
for the ERCC spike-in set and `ENCFF335FFV`_ is the accession ID for
phiX.

.. image:: _static/Screenshot_2019-04-09\ ENCFF717NFD\ –\ ENCODE.png

.. code-block:: console

    wget https://www.encodeproject.org/files/ENCFF001RTP/@@download/ENCFF001RTP.fasta.gz
    wget https://www.encodeproject.org/files/ENCFF335FFV/@@download/ENCFF335FFV.fasta.gz
    wget https://www.encodeproject.org/files/gencode.vM4.tRNAs/@@download/gencode.vM4.tRNAs.gtf.gz
    wget https://www.encodeproject.org/files/gencode.vM4.annotation/@@download/gencode.vM4.annotation.gtf.gz
    merge_encode_annotations -o mm10-M4-male/gencode.vM4-tRNAs-ERCC.gff \
       gencode.vM4.annotation.gtf.gz \
       gencode.vM4.tRNAs.gtf.gz \
       ENCFF001RTP.fasta.gz \
       ENCFF335FFV.fasta.gz

If you'd like, you might want to delete the downloaded files.

.. code-block:: console

    rm ENCFF533JRE.tar.gz ENCFF717NFD.tar.gz ENCFF001RTP.fasta.gz ENCFF335FFV.fasta.gz \
       gencode.vM4.tRNAs.gtf.gz gencode.vM4.annotation.gtf.gz

See also :ref:`howto.building_indexes`

Configuring Paths
-----------------

edit `~/.htsworkflow.ini` with your favorite editor. If you're logged
into a Linux host and don't already have a favorite editor ``nano`` is a
good starting choice.

We need to add some default paths to find software. Using the paths
defined by the above commands we would create file like the followng. 

.. code-block:: ini

    [analysis]
    genome_dir=~/proj/genome/
    georgi_dir=~/proj/GeorgiScripts/
    ucsc_tools_dir=~/proj/ucsc_tools/
    star_dir=/usr/bin
    rsem_dir=/usr/bin

But all of the paths will need to be adjusted for your environment.

star_dir
  needs to be the directory containing the STAR executable
rsem_dir
  needs to be the directory containing the rsem-calculate-expression>

.. _customizing for your computer:

Customizing for your computer
-----------------------------

The underlying condor submit files have been tuned for our specific
cluster, and may not have appropriate settings for your environment
and workload.

If you installed by checking out from git, you can directly modify the
source and then install. However if you installed via pip you'll need
to find the .condor files to customize them.

The following code block should report the installation directory.

.. code-block:: console

    python3 -c 'import os, woldrnaseq; print(os.path.split(woldrnaseq.__file__)[0])'

The downside is that changes will be replaced on upgrade. Though
hopefully we will come up with a way to allow customizing the
request_cpus, request_memory and request_disk settings.

.. _Personal Condor:
   
Personal Condor
---------------

.. code-block:: console
                
    sudo apt install htcondor

Answer Yes to the question "Manage initial HTCondor configuration
automatically." and answer yes to "Perform a Personal HTCondor
installation."

The submit scripts assume that HT Condor is using dynamic slots, and
so you will also need to edit /etc/condor/condor_config.local (or
another valid condor configuration location) and add the following
lines:

.. code-block:: ini

    SLOT_TYPE_1=auto
    SLOT_TYPE_1_PARTITIONABLE=TRUE
    NUM_SLOTS_TYPE_1=1


.. _HT Condor: https://research.cs.wisc.edu/htcondor/
.. _ENCFF717NFD: https://www.encodeproject.org/files/ENCFF717NFD/
.. _ENCFF001RTP: https://www.encodeproject.org/files/ENCFF001RTP/
.. _ENCFF335FFV: https://www.encodeproject.org/files/ENCFF335FFV/
.. _STAR: https://github.com/alexdobin/STAR
.. _rsem: https://deweylab.github.io/RSEM/
