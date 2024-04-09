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
