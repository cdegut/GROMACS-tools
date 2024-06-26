# No ligand
## new version
use the prep.py script.
answer questions
run the .sh file.

## Old (or more granular)
### Prepping the run
Make a folder for the run
Download the latest oplsaam force field in your folder, or symlink to oplsaam.ff
```
ln -s ../../oplsaam.ff/ oplsaam.ff
```
Convert the pdb structure to a .gro file, delete present hydrogens, add new ones.
```
gmx pdb2gmx -ignh -f starting_structure.pdb -o starting_structure.gro
```
Select "1" oplsFF and "1" tip4p water (recomended).\n
Generate a dodecahedron box arround 1A the model
```
gmx editconf -f starting_structure.gro -o struct_new_box -c -d 1.0 -bt dodecahedron
```
Solvate with proper water model
```
gmx solvate -cp struct_new_box.gro -cs tip4p -o struct_box_water.gro -p topol.top
```
Generate ions, this will need to calculate the force field, so need to generate the first tpr file with grompp
```
gmx grompp -f mdpfiles/ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral
```
option -conc ## to set concentration in mol/L

Run slurm script gmxready.slurm
```
sbatch gmxready.slurm
```
or
```
####Run minimization
gmx grompp -f minim.mdp -c struct_box_water_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
##Run NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
##Run NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
##Run test 1ns simulation
gmx grompp -f md1ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md1ns.tpr
gmx mdrun -pme gpu -v -deffnm md1ns
```

# Random things
## Run longer sym
To detach command from terminal:
```
nohup gmx mdrun -pme gpu -v -deffnm md250ns > md250ns.log &; disown
```
or using task-*spooler
```
tsp bash -c "gmx mdrun -pme gpu -v -deffnm md250ns > md250ns.log"
```

## Simple analysis
Center the trajectory
```
gmx trjconv -s md250ns.tpr -f md250ns.xtc -o md250ns_center_po.xtc -center -pbc mol -ur compact
```
Select 1 and 1 to g get only the protein centered on the box, and obtain a much smaller file
Or 1 and 0 to get everything

Extract specific frame:
```
gmx trjconv -s md250ns.tpr -f md250ns_center_po.xtc -o md250ns_start.pdb -dump 0
```
Measure distance to periodic image
```
gmx mindist -s ###.tpr -f ###.xtc -pi
```
## Extend run
```
gmx convert-tpr -s md250ns.tpr -extend 300000 -o md550ns.tpr
gmx mdrun -pme gpu -v -cpi md250ns.cpt -deffnm md550ns -noappend
```
cpt file is the output file with high precision and should be used rather than the gro file with truncated precision

concatenate the rajectories files:
```
gmx trjcat -f *.xtc -o final.xtc
gmx eneconv -f *.edr -o final.edr
```

# Ligand
### Docking
Start by getting ligand, either the structure is simple enough to be directly drawn with  Avogadro.
If structure is complex, may need to use MOPAC to get accurate geometry for the docking
convert to .pdbqt with open babel
Dock to the PDBQT structure with autodock vina (command line)
```
vina –config config.txt
```
Example docking config:
```
receptor = structure.pdbqt
ligand =ITE.pdbqt
center_x = 160
center_y = 165
center_z = 160
size_x = 40
size_y = 40
size_z = 40
```

Extract single molecule from pdbqt (with text editor).
Open with Avogadro2, add hydrogens, export as PDB.
Replace the UNL na min pdb file by a 3 letter name.

Finally itp topology can can be generated with with ligpargen webserver

### Back to runing GROMACS
```
gmx pdb2gmx -ignh -f starting_structure.pdb -o starting_structure.gro
```

Append ITE.gro atoms to starting_structure.gro
example :
```
426PRO     CG 2307  15.227  15.617  14.996
 426PRO    HG1 2308  15.142  15.567  14.978
 426PRO    HG2 2309  15.221  15.660  15.086
 426PRO     CD 2310  15.250  15.720  14.889
 426PRO    HD1 2311  15.207  15.694  14.803
 426PRO    HD2 2312  15.217  15.810  14.918
 426PRO      C 2313  15.562  15.536  14.866
 426PRO     O1 2314  15.532  15.491  14.755
 426PRO     O2 2315  15.626  15.433  14.926
   1LIG    C00    1  15.918  16.400  16.035
   1LIG    O01    2  15.799  16.377  16.027
   1LIG    C02    3  15.999  16.381  16.154
   1LIG    C03    4  15.973  16.279  16.254
   1LIG    C04    5  16.074  16.293  16.355
   1LIG    C05    6  15.875  16.180  16.267
   1LIG    C06    7  16.078  16.210  16.468
   1LIG    N07    8  16.158  16.400  16.316
   1LIG    C08    9  15.980  16.113  16.478
   1LIG    C09   10  15.879  16.098  16.379
   1LIG    C0A   11  16.113  16.452  16.195
```
Increase amount of atom in the first line of the file to corespond to protein + ligand :
```
2345
 284ASN      N    1  16.995  17.030  17.834
 284ASN     H1    2  17.050  17.111  17.855
```
In topol.top file add
```
; Include forcefield parameters
#include "./oplsaam.ff/forcefield.itp"
#include "./LIG.itp"
```

And add at the end of the file add LIG in the \[ molecules \]
```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG                 1
```

dodecahedron box 
```
gmx editconf -f starting_structure.gro -o struct_new_box.gro -c -d 1.0 -bt dodecahedron
```
Solvate
```
gmx solvate -cp  struct_new_box.gro -cs tip4p -o  struct_box_water.gro -p topol.top
```
Ions
```
gmx grompp -f ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral
```
option -conc ## to set concentration in mol/L

If non integer amount of charges : if the difference is small enough and is likely due to rounding error, redistribute the charge on the ligand in the itp file (had to remove 0.003 from ligand)

Run the minimisation separately, it's very short can be run straight on the login server
```
gmx grompp -f mdpfiles/minim.mdp -c solv_ions.gro -p topol.top -o minimisation.tpr
gmx mdrun -v -deffnm minimisation
```

It’s better to put protein and ligand in the same temperature coupling grp, as the more generic grp are protein and non protein:
```
gmx make_ndx -f em.gro -o index.ndx
```
Answer 1 | 19  or whatever protein | ligand number

The grp need to be changed in the mdp files as
```
tc-grps                 = Protein_LIG  Water_and_ions
```
Run slurm script gmxready_ligand.slurm
```
sbatch gmxready_ligand.slurm
```


