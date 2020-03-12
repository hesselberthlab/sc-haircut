#! /usr/bin/env bash

#BSUB -J cov
#BSUB -o log/cov.%J.%I.out
#BSUB -e log/cov.%J.%I.err
##BSUB -R "select[mem>16] span[hosts=1] rusage[mem=16]"

DATE=20181017

##Making directories and such

ANALYSIS=$HOME/projects/10x_haircut/$DATE
index=/beevol/home/ancarr/data-sets/genome/hg38/grcm38_chr.sizes
BAM=$ANALYSIS/geneexp_ampure_37/outs
BED=$ANALYSIS/processed_data/bed_mrna

if [[ ! -d $BED ]]; then
    mkdir -p $BED
fi

# Coverage of BAM output from cellranger
bedtools genomecov -bg -strand - -ibam $BAM/possorted_genome_bam.bam \
    -g $index  \
    > $BED/ampure_37_mrna_neg.bg


