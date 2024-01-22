#!/bin/bash
gmx pdb2gmx -ignh -f $1 -o starting_structure.gro
gmx editconf -f starting_structure.gro -o struct_new_box -c -d 1.0 -bt dodecahedron
gmx solvate -cp struct_new_box.gro -cs tip4p -o struct_box_water.gro -p topol.top
gmx grompp -f /local_scratch/clement/gromacs/GROMACS-tools/mdpfiles/ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral
##Run minimization
gmx grompp -f /local_scratch/clement/gromacs/GROMACS-tools/mdpfiles/minim.mdp -c struct_box_water_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
##Run NVT equilibration
gmx grompp -f /local_scratch/clement/gromacs/GROMACS-tools/mdpfiles/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
##Run NPT equilibration
gmx grompp -f /local_scratch/clement/gromacs/GROMACS-tools/mdpfiles/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
##Run test 1ns simulation
gmx grompp -f /local_scratch/clement/gromacs/GROMACS-tools/mdpfiles/md1ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md1ns.tpr
gmx mdrun -pme gpu -v -deffnm md1ns

##cleanup
mkdir run_prep
mv em* run_prep/
mv nvt* run_prep/
mv npt* run_prep/
mv ions.tpr run_prep/
mv struct_* run_prep/
