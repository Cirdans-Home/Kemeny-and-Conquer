#!/bin/bash
#SBATCH --job-name=kemeny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=10-0:0:0
#SBATCH --output=test_kemeny.out
#SBATCH --mem=250000
#SBATCH --partition=cl2
#SBATCH --exclusive

module load matlab/R2021a

matlab -batch "batch_toeplitz"


