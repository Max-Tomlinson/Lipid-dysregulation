#!/bin/bash -l
#SBATCH --output=/scratch/prj/dtr/Groups_WorkSpace/KerrinSmall/Max/EuroBATS/scripts/output/%j_residual.out
#SBATCH --chdir=/scratch/prj/dtr/Groups_WorkSpace/KerrinSmall/Max/EuroBATS
#SBATCH --time=0-48:00
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mail-user=max.tomlinson@kcl.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

source ~/miniconda3/etc/profile.d/conda.sh

conda activate r

list="$(head -n 1 data/pheno/Residual_q5.txt)"

arr=($list)

for i in {1..24}; do Rscript scripts/R/adipose_residual.R genes/EB_F_STAR2016_TMM_5counts25perc_InvNorm_30012019.txt covar/Adipose_TechCovars_19052017.txt pheno/Residual_q5.txt ${arr[i]}; done

