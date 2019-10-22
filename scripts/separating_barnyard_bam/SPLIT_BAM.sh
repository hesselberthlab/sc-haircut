#! /usr/bin/env bash

#BSUB -J sep_bam
#BSUB -o log/sepbam.%J.%I.out
#BSUB -e log/sepbam.%J.%I.err
#BSUB -n 8
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x
#Split bam files by cell barcodes

DATE=20190819

#Path to repair bam folder
BAM=$HOME/projects/10x_haircut/$DATE/processed_data/bam

#Pathy to mRNA bam
IN=$HOME/projects/10x_haircut/$DATE/JH158_mix_60_mrna/outs/possorted_genome_bam.bam

#Path to list of RNASEH2CKO barcodes as tag strings
CBribo=$HOME/projects/10x_haircut/$DATE/split_bam/barnyard_082019_RNASEH2CKO_cell_tags.tsv

#Path to list of UNGKO barcoes as tag strings
CBuracil=$HOME/projects/10x_haircut/$DATE/split_bam/barnyard_082019_UNGKO_cell_tags.tsv

#Path to header to use when writing BAM files
HEADER=$HOME/projects/10x_haircut/$DATE/processed_data/bam/sep_bam/header.sam

# Make directories
if [ ! -d $BAM/sep_bam ]; then
    mkdir $BAM/sep_bam
fi

# Make and save header
samtools view -H $IN > $HEADER

# Make UNGKO bam
samtools view $IN \
    | grep -f $CBuracil \
    | cat $HEADER - \
    | samtools view -Sb - > $BAM/sep_bam/UNGKO.bam
    
# Make RNASEH2CKO bam
samtools view $IN \
    | grep -f $CBribo \
    | cat $HEADER - \
    | samtools view -Sb - > $BAM/RNASEH2CKO.bam

