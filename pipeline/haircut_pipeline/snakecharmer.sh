#!/usr/bin/env bash

#BSUB -J haircut 
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q normal
set -o nounset -o pipefail -o errexit -x

args=' 
  -q normal 
  -o {log}.out 
  -e {log}.err 
  -J {params.job_name} 
  -R "{params.memory} span[hosts=1] " 
  -n {threads} ' 
    

#### load necessary programs ####

# If programs are not all in the path then modify code to load 
# the necessary programs

# load modules
. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load samtools/1.5
module load bowtie2/2.3.2
module load R/3.5.0
module load cellranger/3.0.2

# other programs (not in modules)
# umi_tools 
# cutadapt

#### execute snakemake ####

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 72 \
    --resources all_threads=72 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config.yaml 
