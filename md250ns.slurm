#!/usr/bin/env bash
#SBATCH --account=bdyrk14
#SBATCH --time=1-12:00:00
#SBATCH --job-name=QH-ITE-MD250


# Node resources:

#SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
#SBATCH --gres=gpu:1       # One GPU per node max of 4 (plus 25% of node CPU and RAM per GPU)
#SBATCH --nodes=1

module load hecbiosim
module load gromacs/2022.2


#gmx grompp -f md250ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md250ns.tpr -n index.ndx
gmx mdrun -pme gpu -v -deffnm md250ns
