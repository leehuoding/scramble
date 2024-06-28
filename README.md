SCRAMble-MEIs Caller
========

Repository Layout
-----------------
1. `src/cluster_identifier.c` - the application responsible for identifying soft clipped clusters. For how to build see the build section.
2. `src/SCRAMble-MEIs.R` - the analysis of soft-clipped clusters and MEIs events.


Build
-----

Install dependencies (Ubuntu 20.04):

    apt-get update
    apt-get install -y  \
        autoconf \
        autogen \
        build-essential \
        curl \
        libbz2-dev \
        libcurl4-openssl-dev \
        libhts-dev \
        liblzma-dev \
        libncurses5-dev \
        libnss-sss \
        libssl-dev \
        libxml2-dev \
        ncbi-blast+ \
        r-base \
        r-bioc-biostrings \
        r-bioc-rsamtools \
        r-cran-biocmanager \
        r-cran-devtools \
        r-cran-stringr \
        r-cran-optparse \
        zlib1g-dev

Install R packages dependencies:

    Rscript -e "library(devtools); install_github('mhahsler/rBLAST')"

To build the cluster_identifier (estimated install time <5 minutes):

    $ cd src
    $ make

That should be it. It will create an executable named `src/cluster_identifier`.

Running
-------
SCRAMble runs as a two-step process. First `cluster_identifier` is used to generate soft-clipped read cluster consensus
sequences. Second, `SCRAMble-MEIs.R` analyzes the cluster file for likely MEIs. Running SCRAMble on the test bam in the validation directory should take <1 minute for each step.

To run SCRAMble cluster_identifier:

    $ /path/to/scramble/src/cluster_identifier -@ thread \
        /path/to/test.bam > /path/to/test.clusters.txt

    $ /path/to/scramble/src/cluster_identifier -@ thread \
        -R /path/to/genome/human_genome.fasta /path/to/test.cram > /path/to/test.clusters.txt

To run SCRAMble-MEIs(with default settings):

    $ Rscript /path/to/scramble/src/SCRAMble-MEIs.R \
        --mei-refs /path/to/scramble/genome/MEI_consensus_seqs.fa \
        --ref /path/to/genome/human_genome.fasta
        --cluster-file /path/to/test.clusters.txt \
        --out-file /path/to/test.vcf.gz 	\


Output
------
The output of cluster_identifier is a tab delimited text file with clipped cluster consensus sequences.
The columns are as follows:

|      |                                          |
| ---: | ---------------------------------------- |
| 1.   | Coordinate                               |
| 2.   | Side of read where soft-clipped occurred |
| 3.   | Number of reads in cluster               |
| 4.   | Clipped read consensus                   |
| 5.   | Anchored read consensus                  |

Calling `SCRAMble.R` with `--eval-meis` produces a tab delimted file. If a reference `.fa` file is provided, then a VCF is produced as well. The `<out-name>_MEIs.txt` output is a tab delimited text file with MEI calls. If no MEIs are present an output file will still be produced with only the header.
The columns are as follows:

|      |                               |                                                                                              |
| ---: | ----------------------------- | -------------------------------------------------------------------------------------------- |
| 1.   | Insertion                     | Coordinate where MEI insertion occurs (zero-based)                                           |
| 2.   | MEI_Family                    | The consensus sequence to which the clipped sequence aligned best                            |
| 3.   | Insertion_Direction           | Whether MEI is on fwd or rev strand relative to bam reference                                |
| 4.   | Clipped_Reads_In_Cluster      | Number of supporting reads in cluster                                                        |
| 5.   | Alignment_Score               | Pairwise alignment score of clipped read consensus to MEI reference sequence                 |
| 6.   | Alignment_Percent_Length      | Percent of clipped read consensus sequence involved in alignment to MEI reference sequence   |
| 7.   | Alignment_Percent_Identity    | Percent identify of alignment of clipped read consensus sequence with MEI reference sequence |
| 8.   | Clipped_Sequence              | Clipped cluster consensus sequences                                                          |
| 9.   | Clipped_Side                  | Left or right, side of read where soft-clipping ocurred                                      |
| 10.  | Start_In_MEI                  | Left-most position of alignment to MEI reference sequence                                    |
| 11.  | Stop_In_MEI                   | Right-most position of alignment to MEI reference sequence                                   |
| 12.  | polyA_Position                | Position of  polyA clipped read cluster if found                                             |
| 13.  | polyA_Seq                     | Clipped cluster consensus sequences of polyA clipped read cluster if found                   |
| 14.  | polyA_SupportingReads         | Number of supporting reads in polyA clipped read cluster if found                            |
| 15.  | TSD                           | Target site duplication sequence if polyA clipped read cluster found                         |
| 16.  | TSD_length                    | Length of target site duplication if polyA clipped read cluster found                        |

Disclaimers
-----------
In theory, SCRAMble should work well on any MEI reference fasta sequences, however, it has only been tested on the
sequences provided in `/path/to/scramble/genome/MEI_consensus_seqs.fa`.
