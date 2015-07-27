# FMRC

Introduction
============

FMRC is a tool for correcting errors in short reads from high-throughput sequencing.
For detailed usage, type `fmrc -h`

References
==========

> S. Greenstein, J. Holt, and L. McMillan. "Accurate Error Correction in Short Reads using FMRC," July 2015.

System Requirements
===================

FMRC has been tested using Python 2.7. The package is distributed with both Cython and the
corresponding C files.  We recommend using installing Cython prior to installing this package.

Four Python modules are required to run the code.

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