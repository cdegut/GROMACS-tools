mdp_files_folder=/local_scratch/clement/gromacs/GROMACS-tools/mdpfiles
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
mkdir pre_run
mv em* pre_run/
mv nvt* pre_run/
mv npt* pre_run/
mv ions.tpr pre_run/
mv struct_* pre_run/
mv starting_structure.gro pre_run/.
mv posre.itp pre_run/.
mv *.log pre_run/.
mv *#* pre_run/.
