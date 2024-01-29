#!/bin/bash
#SBATCH --account=def-baillet
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=21
#SBATCH --time=1-00:00

module load StdEnv/2020 StdEnv/2023 r/4.3.1

Rscript model_bootstrapping_cluster.R
