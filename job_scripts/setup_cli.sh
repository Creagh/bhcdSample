#!/bin/bash
#SBATCH --account=def-bouchar3
#SBATCH --time=00:10:00
#SBATCH --job-name=installDist
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --nodes 1
#SBATCH --mail-user=creagh@stat.ubc.ca
#SBATCH --mail-type=ALL
module load java/1.8.0_192
newgrp def-bouchar3
source setup-cli.sh