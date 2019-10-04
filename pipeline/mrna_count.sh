#! /usr/bin/env bash

#BSUB -J 10xcount[1]
#BSUB -o log/10xcount.%J.%I.out
#BSUB -e log/10xcount.%J.%I.err
#BSUB -n 16
#BSUB -R "select[mem>35] span[hosts=1] rusage[mem=35]"
#BSUB -q normal

set -o nounset -o pipefail -o errexit -x

# path to reference genome
transcriptome="path/to/refdata-cellranger-GRC38-3.0.0/"

# path to fastq files
fastq_path="/path/to/mrna.fastq/"

# List of sample names on fastq files
# For fastq files: pbmc_S1_L001_R*_001.fastq.gz, list the following sample

SAMPLES=(
pbmc
)

sample=${SAMPLES[$(( $LSB_JOBINDEX -1 ))]}


set -x

cellranger count \
    --id=$sample"_mrna" \
    --fastqs=$fastq_path \
    --sample=$sample \
    --localcores=16 \
    --localmem=35 \
    --transcriptome=$transcriptome

