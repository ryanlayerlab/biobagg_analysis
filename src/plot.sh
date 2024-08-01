#!/usr/bin/env bash

### TODO:
### modify these options for your system

#SBATCH -p short
#SBATCH --job-name=plot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.schneider@colorado.edu
#SBATCH --output=log/plot.out
#SBATCH --error=log/plot.err

set -e pipefail
root="/scratch/Users/krsc0813/density/"

python $root"python_scripts/plot_density.py" \
    --ancestry_file $root"1kg_ancestry.tsv" \
    --density_dir $root"data/density_1kg/" \
    --distance_dir $root"data/distance_1kg/" \
    --colors $root"colors.txt" \
    --out ./ \
