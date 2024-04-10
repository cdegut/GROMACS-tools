import default_parameter
from pprint import pp
from copy import deepcopy
import subprocess
import os

def read_mdp(mdp_file_name: str) ->dict:
    mdp_dict = {}
    with open(mdp_file_name, 'r') as mdp_file:
        for line in mdp_file:
            if line[0] != ";":
                try:
                    param = line.split()
                    after = " ".join(param[2:])
                    mdp_dict[param[0]] = after.split(";")
                except:
                    pass
    print("{")
    for values in mdp_dict:
        print(f"\"{values}\"\t\t:  {mdp_dict[values]},")
    print("}")

    return mdp_dict
# mdp =  read_mdp('minim.mdp')

def write_mdp(mdp_dict: dict, filename: str):
    with open(filename, 'w') as file:
        for param in mdp_dict:
            file.write(f"{param}\t\t=\t {mdp_dict[param][0]}\t; {mdp_dict[param][1]}\n")


def set_duration(parameter_set: dict, default: int) -> dict:
    dt = parameter_set["dt"][0]
    duration = input(f"Simulation duration in ns (default:{default}ns):\n").strip() or default
    nsteps = (float(duration) / float(dt)) * 1000
    parameter_set["nsteps"][0] = int(nsteps)
    parameter_set["nsteps"][1] = f"{duration}ns => {int(nsteps)}steps * {dt}ps"
    print(f"Set simulation to last {duration}ns => {int(nsteps)} steps * {dt}ps")

def set_temperature(parameter_set: dict, default: int) -> dict:
    tc_grps_default = parameter_set["tc-grps"][0]
    tc_grps = input(f"Coupling groups (default: {tc_grps_default})\n").strip() or tc_grps_default
    ngrps = len(tc_grps.split())
    temperature = input(f"Simulation temperature in K (default:{default}K):\n").strip() or default
    parameter_set["ref_t"][0] = f"{int(temperature)} " * ngrps
    
    tau_t =	parameter_set["tau_t"][0]
    tau_t = tau_t.split(" ")[0]
    parameter_set["tau_t"][0] = f"{tau_t} " * ngrps

    parameter_set["tc-grps"][0] = tc_grps
    print(f"Set {ngrps} coupling groups at {temperature}K")

def copy_temperature(source:dict , destionation: list):
    for parameter_set in destionation:
        parameter_set["ref_t"] = source["ref_t"]
        parameter_set["tau_t"] = source["tau_t"]
        parameter_set["tc-grps"] = source["tc-grps"]

def list_txt_files():
    pdb_files = [f for f in os.listdir() if os.path.isfile(f) and f.endswith('.pdb')]
    return pdb_files

def select_file(pdb_files):
    print("Please select a file by entering its number:")
    for i, file_name in enumerate(pdb_files):
        print(f"{i+1}. {file_name}")
    
    while True:
        try:
            choice = int(input("Enter your choice: "))
            if choice < 1 or choice > len(pdb_files):
                print("Invalid choice. Please enter a valid number.")
            else:
                return pdb_files[choice - 1]
        except ValueError:
            print("Invalid input. Please enter a number.")




def main():
    pdb_files = list_txt_files()
    if not pdb_files :
        print("No .pdb files found in the current directory.")
        return
    else:
        selected_file = select_file(pdb_files)
        print(f"You have selected: {selected_file}")

    nvt_param = deepcopy(default_parameter.nvt)
    npt_param = deepcopy(default_parameter.nvt)
    md_param = deepcopy(default_parameter.md)
    minim_param = deepcopy(default_parameter.minim)


    set_duration(md_param, 250)
    set_temperature(md_param, 300)
    copy_temperature(md_param, [nvt_param, npt_param])


    duration = float(md_param["nsteps"][0]) * float(md_param["dt"][0]) / 1000
    write_mdp(minim_param, "minim.mdp")
    write_mdp(md_param, f"md{int(duration)}ns.mdp")
    write_mdp(nvt_param, f"nvt.mdp")
    write_mdp(npt_param, f"npt.mdp")
    subprocess.run("touch ions.mdp", shell=True)

    ions = input(f"Salt concentration in M\n").strip() or 0

    commands = [f'gmx pdb2gmx -ignh -f {selected_file} -o starting_structure.gro',
    "gmx editconf -f starting_structure.gro -o struct_new_box -c -d 1.0 -bt dodecahedron",
    "gmx solvate -cp struct_new_box.gro -cs tip4p -o struct_box_water.gro -p topol.top",
    "gmx grompp -f ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr",
    f"gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral -conc {ions}" ]

    for command in commands:
        subprocess.run(command, shell=True)

if __name__ == "__main__":
    main()