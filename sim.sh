#!/bin/bash
mdp_files_folder=/local_scratch/clement/gromacs/GROMACS-tools/mdpfiles
gmx grompp -f $mdp_files_folde/md250ns.mdp -c md1ns.gro -t md1ns.cpt -p topol.top -o md250ns.tpr
gmx mdrun -pme gpu -v -deffnm md250ns
echo "1 1" | gmx trjconv -s md250ns.tpr -f md250ns.xtc -o md250ns_center_po.xtc -center -pbc mol -ur compact
echo "1" | gmx trjconv -s md250ns.tpr -f md250ns_center_po.xtc -o md250ns_po_start.pdb -dump 0
