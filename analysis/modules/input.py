import MDAnalysis.analysis.align as align
import MDAnalysis as mda
from .dssp import dssp

def do_trajectory_CAalignement(atomistic_system, sim_path, trajectory_file_name):
    
    average = align.AverageStructure(atomistic_system, atomistic_system, select='protein and name CA', ref_frame=0)
    print("Calculating averaged structure:")
    average.run(verbose=True)
    averaged_ref = average.results.universe
    # Align all structure on the averaged one
    aligner = align.AlignTraj(atomistic_system, 
                                  averaged_ref, 
                                  select='protein and name CA',
                                  in_memory=False,
                                  filename=sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
    
    print("align trajectory on the averaged one, save as " + sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
    aligner.run(verbose=True)

def res_list(atomistic_system):
    residues = []
    for i, residue in enumerate(atomistic_system.residues.resnames):
        residues.append(f"{atomistic_system.residues.resnames[i]}{atomistic_system.residues.resids[i]}")
    print(residues)

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

def show_fasta(atomistic_system):
    atomistic_system.trajectory[0]
    protein = atomistic_system.select_atoms("protein")
    dsspline_start = dssp(protein)
    fasta = get_fasta(atomistic_system)
    number_line = get_number_line(fasta)
    print(f"{number_line}\n{fasta}\n{dsspline_start}")