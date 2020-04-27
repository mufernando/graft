#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --job-name=testsims
#SBATCH --time=01:00:00
#SBATCH --mem 32G
#SBATCH --output=test.out         ### File in which to store job output
#SBATCH --error=test.err          ### File in which to store job error messages

module use /projects/apps/shared/modulefiles/
module load SLiM/dev

slim -m -t -d N=1000 -d gens=1000 -d "outfile='data/root.trees" recipe.slim
slim -m -t -d N=1000 -d gens=1000 -d "infile='data/root.trees'" -d 'outfile="data/branch1.trees"' recipe.slim
slim -m -t -d N=1000 -d gens=1000 -d "infile='data/root.trees'" -d 'outfile="data/branch2.trees"' recipe.slim

