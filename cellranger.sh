#!/bin/bash

# Grid Engine options (lines prefixed with #$)
#  job name: -N
#  use the current working directory: -cwd
#  number of cores -pe sharedmem
#  runtime limit: -l h_rt
#  memory limit: -l h_vmem
#$ -N 3DMI
#$ -cwd    
#$ -pe sharedmem 8              
#$ -l h_rt=48:00:00 
#$ -l h_vmem=20G
#$ -e error.txt
#$ -o output.txt

# Initialising the environment modules
. /etc/profile.d/modules.sh
# /exports/applications/apps/user-scripts/

# Loading modules

module load igmm/apps/cellranger/5.0.0
module load igmm/apps/STAR/2.7.1a

# Running cellranger

cellranger count --id=E_MTAB_7376 \
--fastqs=/exports/eddie/scratch/zli23/raw/E_MTAB_7376 \
--sample=MI_day3_TIP \
--transcriptome=/exports/eddie/scratch/zli23/reference/refdata-gex-mm10-2020-A
