Installation
============

Assuming you have a system with `HT Condor`_ installed you can install
long rna seq condor with


Setting up system. I tested this with Debian Stretch 9.8, but it
should work similiarly with Ubuntu 16.x or 18.x.

Assuming you don't have a Linux system available first install the
base system according to the vendors instructions.

Next we need to install `HT Condor`_. On older systems it might be
called condor instead. Installing HT Condor properly can be fairly
challenging, but for now follow the defaults for a "Personal Condor
Pool"

.. code-block:: console
                
    sudo apt install htcondor


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

TODO: How do we install STAR, RSEM, UCSC-Tools?

For mm10 M4 we can use these prebuilt indexes provided by the DCC.

.. code-block:: console

   mkdir ~/proj/genome
   cd ~/proj/genome/
   # The STAR index
   wget https://www.encodeproject.org/files/ENCFF533JRE/@@download/ENCFF533JRE.tar.gz
   tar xavf ENCFF533JRE.tar.gz
   # The RSEM Index
   wget https://www.encodeproject.org/files/ENCFF717NFD/@@download/ENCFF717NFD.tar.gz

I found these by looking at a recent experiment of the type I wanted
and looked for the input circles feeding to the program run steps.
(See the bottom of the following screen shot.)

.. image:: _static/Screenshot_2019-04-08\ ENCSR968QHO\ –\ ENCODE.png

To test the installation I picked two of our smallest individual
single cell libraries to treat as if they were bulk samples.

.. code-block:: console

    mkdir ~/proj/rna-test
    cd ~/proj/rna-test
    wget https://www.encodeproject.org/files/ENCFF213IBI/@@download/ENCFF213IBI.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF863CYU/@@download/ENCFF863CYU.fastq.gz

.. note::

   These are from an archived experiment. but they were replaced
   purely for the convience of the DCC. Instead of having many
   individual fastq libraries they wanted all the libraries to be
   merged into a single fastq file that will need demultiplexing.

   However for the purposes of testing it will be far easier to use
   small input files.
   

.. _HT Condor: https://research.cs.wisc.edu/htcondor/
