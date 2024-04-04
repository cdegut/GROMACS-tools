import MDAnalysis.analysis.align as align

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