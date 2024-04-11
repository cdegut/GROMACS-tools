import subprocess
import math 
import os

def Rg(filename):
	'''
	Calculates the Radius of Gyration (Rg) of a protein given its .pdb 
	structure file. Returns the Rg integer value in Angstrom.
	'''
	coord = list()
	mass = list()
	Structure = open(filename, 'r')
	for line in Structure:
		try:
			line = line.split()
			x = float(line[6])
			y = float(line[7])
			z = float(line[8])
			coord.append([x, y, z])
			if line[-1] == 'C':
				mass.append(12.0107)
			elif line[-1] == 'O':
				mass.append(15.9994)
			elif line[-1] == 'N':
				mass.append(14.0067)
			elif line[-1] == 'S':
				mass.append(32.065)
		except:
			pass
	xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
	tmass = sum(mass)
	rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
	mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
	rg = math.sqrt(rr / tmass-mm)
	return(round(rg, 3))

def random_walk_Rgs(n_resids, Khun_lenght=False):
    
    typical_amino_acid = 3.6 # Angstrom
    #Khun_lenght = 8.8 # Angstrom

    theta_Rg =  typical_amino_acid * math.sqrt(n_resids/6)
    expanded_Rg = (typical_amino_acid * n_resids**(3/5) )/ math.sqrt(6)
    collapsed_Rg = (typical_amino_acid * n_resids**(1/3) )/ math.sqrt(6)

    if Khun_lenght:
        theta_Rg = math.sqrt((Khun_lenght / typical_amino_acid)) * theta_Rg 
        expanded_Rg = ((Khun_lenght / typical_amino_acid)**(1-3/5)) * expanded_Rg
        collapsed_Rg = ((Khun_lenght / typical_amino_acid)**(1-1/3)) * collapsed_Rg

    return theta_Rg, expanded_Rg, collapsed_Rg


path = "/local_scratch/clement/gromacs/Cr2RBM-heat/"
gmx_check = subprocess.run(f"gmx check -e {path}md50ns.edr", shell=True,  capture_output=True,  text=True)

for line in gmx_check.stderr.split("\n"):
    if line[0:4] == "Last":
        last = line.split()[-1]
last = int(last.split(".")[0])

if os.path.isfile(f"{path}last_state.pdb"):
    subprocess.run(f"rm {path}last_state.pdb", shell=True)
	
extract_pdb = subprocess.run(f"echo \"1 1\" | gmx trjconv -s {path}md50ns.tpr -f {path}md50ns.xtc --center -pbc mol -o {path}last_state.pdb -dump {last}", 
                             capture_output=True,  text=True, shell=True)

current_Rg = Rg(f"{path}last_state.pdb")
print(f"current Rg {current_Rg}")
theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(129, Khun_lenght=False)
print(f"Theta Rg: {theta_Rg}")