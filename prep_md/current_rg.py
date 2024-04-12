import subprocess
import math 
import os

def Rg(filename):
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

def list_txt_files():
    edr_files = [f for f in os.listdir() if os.path.isfile(f) and f.endswith('.edr')]
    return edr_files

def select_file(edr_files):
    print("Please select a file by entering its number:")
    for i, file_name in enumerate(edr_files):
        print(f"{i+1}. {file_name}")
    
    while True:
        try:
            choice = int(input("Enter your choice: "))
            if choice < 1 or choice > len(edr_files):
                print("Invalid choice. Please enter a valid number.")
            else:
                return edr_files[choice - 1]
        except ValueError:
            print("Invalid input. Please enter a number.")
            
def main():      
	edr_files = list_txt_files()
	if not edr_files :
		print("No .pdb files found in the current directory.")
		return
	
	else:
		selected_file = select_file(edr_files)
		print(f"You have selected: {selected_file}")
			
	gmx_check = subprocess.run(f"gmx check -e {selected_file}", shell=True,  capture_output=True,  text=True)

	for line in gmx_check.stderr.split("\n"):
		print(line)
		if line[0:4] == "Last":
			last = line.split()[-1]
	last = int(last.split(".")[0])

	if os.path.isfile(f"last_state.pdb"):
		subprocess.run(f"rm last_state.pdb", shell=True)
	
	name = selected_file.split(".")[0]

	print("Extracting last frame pdb as last_state.pdb")		
	extract_pdb = subprocess.run(f"echo \"1 1\" | gmx trjconv -s  {name}.tpr -f  {name}.xtc --center -pbc mol -o last_state.pdb -dump {last}", 
								capture_output=True,  text=True, shell=True)
	
	current_Rg = Rg(f"last_state.pdb")
	print(f"current Rg {current_Rg}")
	resid_count = subprocess.run("cat pre_run/starting_structure.gro | grep CA | wc -l", shell=True,  capture_output=True,  text=True)
	n_resid = int(resid_count.stdout)
	theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(n_resid, Khun_lenght=False)
	print(f"Theta Rg: {theta_Rg:.3} Collapsed Rg: {collapsed_Rg:.3}")

main()