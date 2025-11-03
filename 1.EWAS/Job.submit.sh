#!/bin/bash
#SBATCH --job-name=demo
#SBATCH --time=06:00:00
#SBATCH -p x86_64_cpu
#SBATCH --nodes=1
#SBATCH -o R_%A_%a.out # Standard output
#SBATCH -e R_%A_%a.err # Standard error

source /share/opt/anaconda3/bin/activate r_env # activate the environment

cd /share/users/qicancan/DNAmethylation/demo # change to working directory

Rscript EWAS_demo_script.R \
        --input DNAm.meta.demo.csv \ # input meta data
        --output results.demo.txt \ # name of output file
        --cols CRF,Age,Gender,BMI \ # colunms invovled in the analysis, the first one is the phenotype you are interested
        --verbose
