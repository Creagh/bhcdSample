#!/bin/bash
#SBATCH --account=def-bouchar3
#SBATCH --time=00:20:00
#SBATCH --job-name=bhcd
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=2
#SBATCH --mail-user=creagh@stat.ubc.ca
#SBATCH --mail-type=ALL
module load java/1.8.0_192
module load gcc/7.3.0 r/4.0.0
java -cp build/install/blangBHCD/lib/\* bhcd.BHCD --model.graph.file data/example/example.csv --experimentConfigs.recordGitInfo false --checkIsDAG false --engine.usePriorSamples true --engine PT --engine.nChains 2 --engine.nScans 10 --engine.initialization COPIES --treatNaNAsNegativeInfinity true --postProcessor DefaultPostProcessor