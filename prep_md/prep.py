import default_parameter
import ff_patch
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
            file.write(f"{param}\t\t=\t{mdp_dict[param][0]}\t;")
            if len(mdp_dict[param]) > 1:
                for i in range(1, len(mdp_dict[param])):
                    file.write(mdp_dict[param][i])
            file.write("\n")


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

def select_ff():
    
    gromacs_path = subprocess.run("echo $GROMACS_DIR",shell=True, stdout = subprocess.PIPE, text=True)
    pathlist = os.listdir(f"{gromacs_path.stdout.rstrip()}/share/gromacs/top/")
    pathlist = [x for x in pathlist if x.endswith('.ff')]
    pathlist.sort()
    ff_dict = {}

    i = 0
    for ff in pathlist:
        with open(f"{gromacs_path.stdout.rstrip()}/share/gromacs/top/{ff}/forcefield.doc") as f:
            ff_dict[i] =  ff , f.readline().rstrip()
            i = i+1

    print("Please select a force field by entering its number:")
    for i in ff_dict:
        print(f"{i}. {ff_dict[i][1]}")
    
    while True:
        try:
            choice = int(input("Enter your choice: "))
            if choice > len(ff_dict)-1:
                print("Invalid choice. Please enter a valid number.")
            else:
                return ff_dict[choice]
        except ValueError:
            print("Invalid input. Please enter a number.")


def apply_ff_patch(ff, parameter_set):
    print(ff)
    if ff in ff_patch.patch_list:
        print(ff_patch.patch_list[ff])
        new_param = parameter_set | ff_patch.patch_list[ff]
        print("patch applied !!")
        return new_param
    else:
        return parameter_set


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

def select_water_model():
    model_list = ["spc216", "tip4p",  "tip5p"]
    water_model = 0
    while int(water_model) not in [1,2,3]:
        water_model = input(f"Choose 1,2 o r3:\nFill box with:\n1: tip3p spc216 (default gromos/amber)\n2: tip4p 216 (default opls)\n3: tip5p216\n").strip() or 2

    return model_list[int(water_model)-1]

def select_box_type():
    box_list = ["triclinic", "cubic", "dodecahedron", "octahedron"]
    sel_box = 0

    while int(sel_box) not in [1,2,3,4]:
        sel_box = input(f"Choose box type: \n1: triclinic \n2: cubic\n3: dodecahedron (default)\n4: octahedron").strip() or 3
    
    box_type = box_list[int(sel_box)-1]

    solvent_distance = input(f"Solvent distance in nm (1.2nm)").strip() or 1.2

    return box_type, solvent_distance


def make_runsh(md_name: str):
    commands_list = [
    "gmx grompp -f  minim.mdp -c struct_box_water_ions.gro -p topol.top -o em.tpr",
    "gmx mdrun -v -deffnm em >> em.out 2>&1",
    "##Run NVT equilibration",
    "gmx grompp -f  nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr",
    "gmx mdrun -v -deffnm nvt >> nvt.out 2>&1",
    "##Run NPT equilibration",
    "gmx grompp -f  npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr",
    "gmx mdrun -v -deffnm npt >> npt.out 2>&1",
    f"gmx grompp -f  {md_name}.mdp -c npt.gro -t npt.cpt -p topol.top -o {md_name}.tpr",
    "##cleanup",
    "mkdir pre_run",
    "mv em* pre_run/",
    "mv nvt* pre_run/",
    "mv npt* pre_run/",
    "mv ions* pre_run/",
    "mv struct_* pre_run/",
    "mv starting_structure.gro pre_run/.",
    "mv posre.itp pre_run/.",
    "mv *.log pre_run/.",
    "mv *#* pre_run/.",
    "##Run proper simulation",
    f"gmx mdrun -pme gpu -v -deffnm {md_name} >> {md_name}.out 2>&1",
    "## generate light protein only files for analysis",
    f"echo \"1 1\" | gmx trjconv -s {md_name}.tpr -f {md_name}.xtc -o {md_name}_center_po.xtc -center -pbc mol -ur compact",
    f"echo \"1\" | gmx trjconv -s {md_name}.tpr -f {md_name}_center_po.xtc -o {md_name}_po_start.pdb -dump 0",
    f"##### Read logfiles using: tail -f {md_name}.out"
    ]

    with open(f"run_{md_name}.sh", 'w') as file:
        for command in commands_list:
            file.write(f"{command}\n")
    subprocess.run(f"chmod +x run_{md_name}.sh", shell=True)    

def run_commands_list(commands):

    for command in commands:
        cmd = subprocess.run(command, shell=True)
        if cmd.returncode != 0:
            status = False
            print(f"Last command exit code {cmd.returncode}")
            return False
        
    return True


def prep_run(minim_param,nvt_param,npt_param,md_param):

    pdb_files = list_txt_files()
    if not pdb_files :
        print("No .pdb files found in the current directory.")
        return
    else:
        selected_file = select_file(pdb_files)
        print(f"You have selected: {selected_file}")

    set_duration(md_param, 250)
    set_temperature(md_param, 300)
    copy_temperature(md_param, [nvt_param, npt_param])


    box_type, solvent_distance = select_box_type()

    force_field = select_ff()

    nvt_param = apply_ff_patch(force_field[0][:-3], nvt_param)
    npt_param = apply_ff_patch(force_field[0][:-3], npt_param)
    md_param = apply_ff_patch(force_field[0][:-3], md_param)
    minim_param = apply_ff_patch(force_field[0][:-3], minim_param)

    box = [f'gmx pdb2gmx -ignh -f {selected_file} -o starting_structure.gro -ff {force_field[0][:-3]}',         
    f"gmx editconf -f starting_structure.gro -o struct_new_box -c -d {solvent_distance} -bt {box_type}"]

    box_run = run_commands_list(box)
    
    if not box_run:
        return
    
    ions = input(f"Salt concentration in M\n").strip() or 0
    water_model = select_water_model()
    solvent = [
    f"gmx solvate -cp struct_new_box.gro -cs {water_model} -o struct_box_water.gro -p topol.top",
    "touch ions.mdp",
    "gmx grompp -f ions.mdp -c struct_box_water.gro -p topol.top -o ions.tpr -maxwarn 2",
    f"gmx genion -s ions.tpr -o struct_box_water_ions.gro -p topol.top -pname NA -nname CL -neutral -conc {ions}" 
    ]

    solvent_run = run_commands_list(solvent)

    if not solvent_run:
        return
    
    duration = float(md_param["nsteps"][0]) * float(md_param["dt"][0]) / 1000
    write_mdp(minim_param, "minim.mdp")
    write_mdp(md_param, f"md{int(duration)}ns.mdp")
    write_mdp(nvt_param, f"nvt.mdp")
    write_mdp(npt_param, f"npt.mdp")
    md_name =  f"md{int(duration)}ns"
    make_runsh(md_name)

    if os.path.isfile(f"run_{md_name}.sh"):
        print(f"run_{md_name}.sh file created, ready to execute")
        
    return {'pdb': selected_file,'box type': box_type, 'solvant distance': solvent_distance, 
            'Force Field': force_field, 'Salt concentration': ions, 'Initial duration':duration}
    

if __name__ == "__main__":

    nvt_param = deepcopy(default_parameter.nvt)
    npt_param = deepcopy(default_parameter.nvt)
    md_param = deepcopy(default_parameter.md)
    minim_param = deepcopy(default_parameter.minim)

    parameters = prep_run(minim_param,nvt_param,npt_param,md_param)
    with open('input_parameters.log', "w") as f:
        for param in parameters:
            f.write(f"{param}   {parameters[param]}\n")