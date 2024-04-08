import MDAnalysis as mda
from MDAnalysis.auxiliary import EDR
import warnings
import os
import modules.input

### Needed input
#sim_path = "S:\\Work\\gromacs\\Cr_mutants\\APAx2\\"
sim_path = "S:\\Work\\gromacs\\Cr2RBM-2\\"
#sim_path = "D:\\MD\\Cr2RBM-2\\"
ligand_name = None #set to None if no ligand
coordinate_file_name = 'md250ns_po_start.pdb'             ##pdb / gro
trajectory_file_name = 'md250ns_center_po.xtc'      ##xtc file
auxiliary_file_name = 'md250ns.edr'              ##edr file
reference_structure_path = 'md1ns.gro'
sim_name = 'default' ##output pdb file will use thise name, if set to "default" will use the auxilliary file name
#############################################################
#############################################################





if sim_name == 'default':
    sim_name = auxiliary_file_name.split(".")[0]

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

def get_fasta(atomistic_system: mda.Universe):
    fasta = ""
    for i, residue in enumerate(atomistic_system.residues.resnames):
        fasta = fasta + mda.lib.util.convert_aa_code(atomistic_system.residues.resnames[i])
    return fasta

def get_number_line(fasta):
    number_line = ""
    for i in range(0,len(fasta),5):
        i = str(i)
        space = "-" * (5 - len(i))
        number_line = number_line + i + space
    number_line = number_line[1:len(fasta)+1]
    return number_line
