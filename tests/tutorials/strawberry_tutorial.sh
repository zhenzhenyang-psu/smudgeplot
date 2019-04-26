#!/bin/bash

mkdir -p strawberry_iinumae && cd strawberry_iinumae
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR013/DRR013884/DRR013884_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR013/DRR013884/DRR013884_2.fastq.gz

mkdir tmp # create a directory for temporary files
ls DRR013884_1.fastq.gz DRR013884_2.fastq.gz > FILES # create a file with both the raw read files
kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmer_counts tmp # run kmc
kmc_dump -ci100 -cx3000 kmer_counts kmer_k21.dump

smudgeplot hetkmers -k 21 -o kmer_pairs < kmer_k21.dump