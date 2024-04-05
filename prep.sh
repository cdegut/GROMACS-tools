#!/bin/bash

mdp_files_folder=/local_scratch/clement/gromacs/GROMACS-tools/mdpfiles


if [ -z "$1" ]; then
    echo "Usage: $0 input.pdb"
    exit 1
fi

gmx pdb2gmx -ignh -f $1 -o starting_structure.gro
gmx editconf -f starting_structure.gro -o struct_new_box -c -d 1.0 -bt dodecahedron
gmx solvate -cp struct_new_box.gro -cs tip4p -o struct_box_water.gro -p topol.top
gmx grompp -f $mdp_files_folder/ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr

echo "Salt Concnetration in M:"
read userInput
gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral -conc $userInput

##Run minimization
gmx grompp -f  $mdp_files_folder/minim.mdp -c struct_box_water_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
##Run NVT equilibration
gmx grompp -f  $mdp_files_folder/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
##Run NPT equilibration
gmx grompp -f  $mdp_files_folder/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
##Run test 1ns simulation
gmx grompp -f  $mdp_files_folder/md1ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md1ns.tpr
gmx mdrun -pme gpu -v -deffnm md1ns

##cleanup
mkdir run_prep
mv em* run_prep/
mv nvt* run_prep/
mv npt* run_prep/
mv ions.tpr run_prep/
mv struct_* run_prep/
mv starting_structure.gro run_prep/.
mv posre.itp run_prep/.
mv *.log run_prep/.
mv *#* run_prep/.
