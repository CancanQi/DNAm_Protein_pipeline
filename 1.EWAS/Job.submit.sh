#!/bin/bash
#SBATCH --job-name=demo
#SBATCH --time=06:00:00
#SBATCH -p x86_64_cpu
#SBATCH --nodes=1
#SBATCH -o R_%A_%a.out # Standard output
#SBATCH -e R_%A_%a.err # Standard error

source /share/opt/anaconda3/bin/activate r_env

cd /share/users/qicancan/DNAmethylation/demo

Rscript EWAS_demo_script.R \
        --input DNAm.meta.demo.csv \
        --output results.demo.txt \
        --cols CRF,Age,Gender,BMI \
        --verbose
