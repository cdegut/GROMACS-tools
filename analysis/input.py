import MDAnalysis as mda
from MDAnalysis.auxiliary import EDR
import warnings
import os
import modules.input


### Needed input
#sim_path = "S:\\Work\\gromacs\\SasY\\double"
##sim_path = "S:\\Work\\gromacs\\Cr2RBM-2\\"
sim_path = "D:\\MD\\Shirt\\r3-4"
name= 'md1000ns'
ligand_name = None #set to None if no ligand
coordinate_file_name = f'{name}_po_start.pdb'             ##pdb / gro
trajectory_file_name = f'{name}_center_po.xtc'      ##xtc file
auxiliary_file_name = f'{name}.edr'              ##edr file
reference_structure_path = f'{name}_po_start.pdb'
sim_name = name ##output pdb file will use thise name, if set to "default" will use the auxilliary file name
#############################################################
#############################################################




if sim_path[-1] != "\\":
    sim_path = sim_path + "\\"

if sim_name == 'default':
    sim_name = auxiliary_file_name.split(".")[0]
#MARK: test
def read_files(edr_only = False):
    ## Read the diferents files
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if not edr_only:

            reference_structure = mda.Universe(sim_path + reference_structure_path,)
            if os.path.isfile(sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc"): ## read the aligned file if it exist
                atomistic_system = mda.Universe(sim_path + coordinate_file_name, sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
                print("Aligned file found, reading this instead")
                is_aligned = True
            else:
                atomistic_system = mda.Universe(sim_path + coordinate_file_name, sim_path + trajectory_file_name)
                is_aligned = False
            energy_like_terms = EDR.EDRReader(sim_path + auxiliary_file_name)
            return atomistic_system, reference_structure, energy_like_terms, is_aligned
        
        else:
            return EDR.EDRReader(sim_path + auxiliary_file_name)

def align_traj(is_aligned):
    if not is_aligned:
        atomistic_system = mda.Universe(sim_path + coordinate_file_name, sim_path + trajectory_file_name)
        modules.input.do_trajectory_CAalignement(atomistic_system, sim_path, trajectory_file_name)
        atomistic_system = mda.Universe(sim_path + coordinate_file_name, sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
        is_aligned = True
    else:
        print(f"Trajectory allready aligned, delete {trajectory_file_name.split('.')[0] + '_aligned.xtc'} to rerun alignement")
    return is_aligned

