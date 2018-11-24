pipelign
========

[![GitHub license](https://img.shields.io/github/license/asmmhossain/pipelign.svg)](./LICENSE)
[![Requires.io](https://img.shields.io/requires/github/asmmhossain/pipelign.svg)](https://requires.io/github/asmmhossain/pipelign/requirements/)
[![Travis](https://img.shields.io/travis/asmmhossain/pipelign.svg)](https://travis-ci.org/asmmhossain/pipelign)

Built from [makenew/python-package](https://github.com/makenew/python-package).

Description
-----------

A pipeline for automated multiple sequence alignment, particularly of viral sequences.

Citation
--------

Pipelign: an alignment pipeline for viral sequences. A.S.Md.M. Hossain and S.D.W.Frost, in preparation.

Usage
-----

```sh
usage: pipelign [-h] -i INFILE -o OUTFILE [-t LENTHR] [-a {dna,aa,rna}] [-f]
                [-b] [-z] [-p SIMPER] [-r {J,G}] [-e {P,C}] [-q THREAD]
                [-s MITERATELONG] [-m MITERATEMERGE] -d OUTDIR [-c]
                [-w AMBIGPER] [-n {1,2,3,4,5,6}] [-x]

Pipelign: creates multiple sequence alignment from FASTA formatted sequence file

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --inFile INFILE
                        Input sequence file in FASTA format
  -o OUTFILE, --outFile OUTFILE
                        FASTA formatted output alignment file
  -t LENTHR, --lenThr LENTHR
                        Length threshold for full sequences (default: 0.7)
  -a {dna,aa,rna}, --alphabet {dna,aa,rna}
                        Input sequences can be dna/rna/aa (default: dna)
  -f, --keepOrphans     Add fragments without clusters
  -b, --keepBadSeqs     Add long sequences with too many ambiguous residues
  -z, --mZip            Create zipped intermediate output files
  -p SIMPER, --simPer SIMPER
                        Percent sequence similarity for clustering (default: 0.8)
  -r {J,G}, --run {J,G}
                        Run either (J)oblib/(G)NU parallel version (default: G)
  -e {P,C}, --merge {P,C}
                        Merge using (P)arallel/(C)onsensus strategy  (default: P)
  -q THREAD, --thread THREAD
                        Number of CPU/threads to use (default: 1)
  -s MITERATELONG, --mIterateLong MITERATELONG
                        Number of iterations to refine long alignments (default: 1)
  -m MITERATEMERGE, --mIterateMerge MITERATEMERGE
                        Number of iterations to refine merged alignment (default: 1)
  -d OUTDIR, --outDir OUTDIR
                        Name for output directory to hold intermediate files
  -c, --clearExistingDirectory
                        Remove files from existing output directory
  -w AMBIGPER, --ambigPer AMBIGPER
                        Proportion of ambiguous characters allowed in the long sequences (default: 0.1)
  -n {1,2,3,4,5,6}, --stage {1,2,3,4,5,6}
                        1  Make cluster alignments and HMM of long sequences
                        2  Align long sequences only
                        3  Assign fragments to clusters
                        4  Make cluster alignments with fragments
                        5  Align all sequences
  -x, --excludeClusters
                        Exclude clusters from final alignment
```

In addition, a utility to convert GenBank files into plain FASTA files with the accession as header is included as `gb2fas`.

Dependencies
------------

- MAFFT
- HMMER3
- CD-HIT
- IQTREE

These can be installed e.g. using conda from the bioconda channel:

    $ conda install mafft hmmer cd-hit iqtree -c bioconda

Installation with pip
---------------------

Install it directly using pip with

    $ git clone https://github.com/asmmhossain/pipelign
    $ cd pipelign
    $ pip install .

Installation with setuptools
----------------------------

    $ python3 setup.py install

Development and Testing
-----------------------

### Source Code

The [pipelign source](https://github.com/asmmhossain/pipelign) is hosted
on GitHub. Clone the project with

    $ git clone https://github.com/asmmhossain/pipelign.git

### Requirements

You will need [Python 3](https://www.python.org/) with
[pip](https://pip.pypa.io/).

Install the development dependencies with

    $ pip install -r requirements.devel.txt

### Building a conda package

Installation with conda
-----------------------

First create an environment.

    $ conda create -n pipelign python=3

Activate the environment.

    $ source activate pipelign

Install conda-build:

    $ conda install conda-build

Run conda.

    $ conda-build . -c bioconda

### Tests

Lint code with

    $ python setup.py lint

Run tests with

    $ python setup.py test

Contributing
------------

Please submit and comment on bug reports and feature requests.

To submit a patch:

1.  Fork it (<https://github.com/asmmhossain/pipelign/fork>).
2.  Create your feature branch (`git checkout -b my-new-feature`).
3.  Make changes. Write and run tests.
4.  Commit your changes (`git commit -am 'Add some feature'`).
5.  Push to the branch (`git push origin my-new-feature`).
6.  Create a new Pull Request.

License
-------

This Python package is licensed under the MIT license.

Warranty
--------

This software is provided \"as is\" and without any express or implied
warranties, including, without limitation, the implied warranties of
merchantibility and fitness for a particular purpose.
