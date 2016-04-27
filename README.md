# FMRC

Introduction
============

FMRC is a command line tool for correcting errors in DNA short reads from high-throughput sequencing. It uses a Burrows-Wheeler Transform and FM-index to enable a k-mer counting approach for correcting substitution, insertion, and deletion errors. In general, it corrects errors more effectively than other error correction tools, leading to better alignments and de novo assemblies. It takes as input a FASTA or FASTQ file containing reads as well as a multi-string Burrows-Wheeler Transform built from the reads, and outputs a FASTA or FASTQ file containing the corrected reads.

References
==========

> S. Greenstein, J. Holt, and L. McMillan, “Short read error correction using an fm-index,” in Bioinformatics and Biomedicine (BIBM), 2015 IEEE International Conference on. IEEE, 2015, pp. 101–104.

System Requirements
===================

FMRC has been tested using Python 2.7.

Five Python modules are required to run the code.


* cython - Tested with cython 0.20.1

Cython is an optimising static compiler for both the Python programming language and the extended Cython programming
language. Its latest package can be downloaded from
http://cython.org/#download


* numpy - Tested with numpy 1.9.1

NumPy is a scientific computing package for Python. Its latest package can be downloaded from
http://sourceforge.net/projects/numpy/files/


* msbwt - Tested with msbwt 0.2.5

msbwt is used to build the BWT and FM-index that FMRC uses to count _k_-mers. Its latest package can be downloaded from
http://github.com/holtjma/msbwt


* pysam - Tested with pysam 0.8.3

As a wrapper of Samtools, the pysam module facilitates the manipulation of SAM/BAM files in Python. Its latest
package can be downloaded from https://github.com/pysam-developers/pysam


* argparse - Tested with argparse 1.1

The argparse module is used to parse the command line arguments of the module. It has been maintained in Python
Standard Library since Python 2.7.  Its latest package can be downloaded from http://code.google.com/p/argparse/

Installation
============

It is recommended to use easy-install (http://packages.python.org/distribute/easy_install.html) for the
installation.

	easy_install fmrc

Alternatively, users can download the source from https://github.com/sgreenstein/fmrc

and then type:

	easy_install fmrc-master.zip

By default, the package will be installed under the directory of Python dist-packages, and the executable of
fmrc can be found under '/usr/local/bin/'.

If you don't have permission to install it in the system-owned directory, you can install it in locally following
the next steps:

(1) Create a local package directory for python:

	mkdir -p <local_dir>

(2) Add the absolute path of <local_dir> to the environment variable PYTHONPATH:

	export PYTHONPATH=$PYTHONPATH:<local_dir>

(3) Use easy_install to install the package in that directory:

	easy_install -d <local_dir> fmrc-<version>.tar.gz

For example, if you want to install the package under the home directory in
a Linux system, you can type:

	mkdir -p /home/$USER/.local/lib/python/dist-packages/
	export PYTHONPATH=$PYTHONPATH:/home/$USER/.local/lib/python/dist-packages/
	easy_install -d /home/$USER/.local/lib/python/dist-packages/ fmrc-<version>.tar.gz

After installation, fmrc will be located in '/home/$USER/.local/lib/python/dist-packages/'.

Parameter Selection
===================

* _k_-mer Multiplicity Threshold

The number of times a _k_-mer appears in the reads must exceed this threshold in order for the _k_-mer to be trusted.
The optimal setting depends on the coverage. The default threshold of 2 is recommended, but a threshold of 1 may be
tried for low-coverage samples or a higher threshold for higher coverage.

* _k_

The length of the _k_-mers to use. Longer _k_-mers are better able to distinguish between similar sequence, making
it less likely that one portion of the genome will be "corrected" to look like another. If a read has no _k_ consecutive
error-free bases, it cannot be corrected, so shorter _k_-mers are better able to correct reads with many errors. Small
changes to _k_ do not drastically alter the results.

* Read Filtering

This option filters out reads that are suspected to contain errors, none of which can be corrected. This occurs when
there exist no _k_ consecutive error-free bases in the read. This option is recommended.

Detailed Description
===========

FMRC takes two inputs:

(1) FASTQ or FASTA files containing reads

This file contains the reads that you wish to correct. Each read in a FASTA file has the following specification: A single-line description beginning with a ">" symbol, followed by one or more lines containing the sequence of the read. Each read in FASTQ format must have: (1) A with the read identifier preceded by the "@" symbol, (2) the read sequence, (3) a line starting with the "+" that may or may not repeat the read identifier, and (4) the quality values for the read. Further specification of FASTA and FASTQ files can be found at https://en.wikipedia.org/wiki/FASTQ_format and https://en.wikipedia.org/wiki/FASTA_format. FMRC does not take into account pairing for paired-end reads. Multiple FASTA/FASTQ files can be corrected in series.

(2) A Multi-string Burrows-Wheeler Transform (MSBWT) built from all reads from a sample

To make the MSBWT, the msbwt package must be installed using the instructions on http://github.com/holtjma/msbwt. The `cffq` function must be used to build the MSBWT from all the FASTQ/FASTA files for a sample.

After the MSBWT is built using this function, FMRC can be run to output corrected versions of the FASTA/FASTQ files. The directory used as the input argument `outBwtDir` for `msbwt cffq` must be passed as the `bwtDir` argument to FMRC.

The output of FMRC is the corrected reads in a new file that is the same format as the input file. For instance, if the original reads are in FASTQ format, the output will be a new FASTQ file containing the corrected reads.

Example
=======

Say we wish to correct the reads from a sample. The reads are all 100bp and are stored in FASTQ format. The FASTQ file containing the reads is called raw_example.fastq. After installing msbwt and fmrc, we build the MSBWT:

    msbwt cffq -p 4 -u -c bwt uncorrected_example.fastq

This builds a compressed MSBWT in the directory "bwt" using 4 concurrent processes. Now we can run fmrc:

    fmrc --filter -o corrected_example.fastq -p 4 bwt 100 raw_example.fastq

This corrects the reads in raw_example.fastq using 4 concurrent processes. It filters out reads with uncorrectable errors, and outputs the corrected reads to corrected_example.fastq.

Both the input and output FASTQ files used in this example are included in the source, allowing you to verify your installation is working correctly.

