#!/usr/bin/env bash
#SBATCH --account=bdyrk14
#SBATCH --time=10:00:00
#SBATCH --job-name=continue


# Node resources:

#SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
#SBATCH --gres=gpu:1       # One GPU per node max of 4 (plus 25% of node CPU and RAM per GPU)
#SBATCH --nodes=1

module load hecbiosim
module load gromacs/2022.2


gmx mdrun -cpi md250ns.cpt -deffnm md250ns
