#!/usr/bin/env bash

### TODO:
### modify these options for your system

#SBATCH -p short
#SBATCH --job-name=density
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32gb
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.schneider@colorado.edu
#SBATCH --output=log/density.out
#SBATCH --error=log/density.err

set -e pipefail

python compute_density.py \
    --encodings '/scratch/Users/krsc0813/chr1_22/encodings/' \
    --chrm 8 \
    --out './chrm8_dens/'
