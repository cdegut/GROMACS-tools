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

simulation_name=md550ns


#gmx grompp -f $simulation_name.mdp -c npt.gro -t npt.cpt -p topol.top -o $simulation_name.tpr -n index.ndx
gmx mdrun -pme gpu -v -deffnm $simulation_name

#do the trajectory alignement

echo "0 1" | gmx trjconv -s $simulation_name.tpr -f $simulation_name.xtc -o $simulation_name_center.xtc -center -pbc mol -ur compact
