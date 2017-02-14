# Pipelign

Semi-automated computational pipeline for aligning virus sequences.

``Pipelign`` can be downloade from github using the following command:

``git clone https://github.com/asmmhossain/pipelign``

Typing ``python3 Pipelign.py -h`` will give a short usage information

Dependencies:

- [Python 3](https://www.python.org/)
- [BioPython](http://biopython.org/wiki/Main_Page): for sequence processing
- [CD-HIT](http://weizhongli-lab.org/cd-hit/): for clustering 'long' sequences
- [IQ-TREE](http://www.cibiv.at/software/iqtree/): for generating a phylogenetic tree of representative 'long' sequences
- [MAFFT](http://mafft.cbrc.jp/alignment/software/): for alignment
- [HMMER](http://hmmer.org/): for matching sequence fragments to clusters
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): for matching sequence fragments to clusters

## Docker

Docker repository for ``Pipelign`` is hosted on [asmmhossain/Pipelign](https://hub.docker.com/r/asmmhossain/pipelign/). Anyone having [Docker](https://www.docker.com/) installed can get ``Pipelign`` using the command:

``docker pull asmmhossain/pipelign``

The following command will show the usage information and command line argument options for running ``Pipelign`` using docker:

``docker run --rm asmmhossain/pipelign Pipelign.py -h``

Then ``Pipelign`` can be run to produce alignments using desired options:

``docker run -i -t --rm -v $(pwd):/data -w /data asmmhossain/pipelign Pipelign.py [options]``

Please make sure the sequence file resides in the current directory from where the docker container is launched. Output file will be stored in the current directory.


