#!/usr/bin/env bash
#SBATCH --account=bdyrk14
#SBATCH --time=1:00:00
#SBATCH --job-name=gmxPrep


# Node resources:

#SBATCH --partition=infer    # Choose either "gpu" or "infer" node type
#SBATCH --gres=gpu:1       # One GPU per node max of 4 (plus 25% of node CPU and RAM per GPU)
#SBATCH --nodes=1

module load hecbiosim
module load gromacs/2022.2

##Run NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx
gmx mdrun -v -deffnm nvt

##Run NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index.ndx
gmx mdrun -v -deffnm npt

##Run test 1ns simulation
gmx grompp -f md1ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md1ns.tpr -n index.ndx
##gmx mdrun -update gpu -pme gpu -v -deffnm md1ns update gpu doesn't work with virtual site (4p water)
