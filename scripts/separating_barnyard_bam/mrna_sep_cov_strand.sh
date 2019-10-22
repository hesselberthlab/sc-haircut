#! /usr/bin/env bash

#BSUB -J cov[1-2]
#BSUB -o log/cov.%J.%I.out
#BSUB -e log/cov.%J.%I.err
##BSUB -R "select[mem>16] span[hosts=1] rusage[mem=16]"

samples=(
RNASEH2CKO
UNGKO)

id=${samples[$(( $LSB_JOBINDEX -1 ))]}

#s=MIX

DATE=20181017

# Path to analysis home
ANALYSIS=$HOME/projects/10x_haircut/$DATE

# Path to genome index file
index=/beevol/home/ancarr/data-sets/genome/hg38/grcm38_chr.sizes

# Path to folder of BAM files
BAM=$ANALYSIS/processed_data/bam

# Path to output bed folder
BED=$ANALYSIS/processed_data/bed/mrna

if [[ ! -d $BED ]]; then
    mkdir -p $BED
fi

# Minus strand coverage
bedtools genomecov -bg -strand - -ibam $bam/$id.sort.bam \
    -g $index.genome  \
    > $bed/$id.mrna.minus.bedgraph


# Plus strand coverage
bedtools genomecov -bg -strand + -ibam $bam/$id.sort.bam \
    -g $index.genome  \
    > $bed/$id.mrna.plus.bedgraph
