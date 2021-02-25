#!/bin/bash
#SBATCH --account=def-bouchar3
#SBATCH --time=72:00:00
#SBATCH --job-name=bhcd-gw
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=48
#SBATCH --mail-user=creagh@stat.ubc.ca
#SBATCH --mail-type=ALL
module load java/1.8.0_192
module load gcc/7.3.0 r/4.0.0
java -cp build/install/blangBHCD/lib/\* bhcd.BHCD --model.graph.file data/grass_web/grass_web_reduced.csv --experimentConfigs.recordGitInfo false --treatNaNAsNegativeInfinity true --checkIsDAG false --engine PT --engine.nChains 50 --engine.nScans 250000 --engine.initialization COPIES --engine.usePriorSamples true --postProcessor DefaultPostProcessor --model.theta 0.0 --engine.thinning 10 --engine.nPassesPerScan 9