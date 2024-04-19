# Python script to autogenerate input for GROMACS
## prep.py
Run inside the folder were simulation is to be ran, 
generate all the mdp files for the equilibration and run.
The default mdp option are pulled from dictionaries in 
'''default_parameter.py'''
also generate a .sh file to queue all the computational intensive parts.
## curent_RG.py
get the Rg of the last step of the simulation in folder (can be on runing simulation or not) 
and dump it as a PDB file
