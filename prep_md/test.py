import subprocess
import ff_patch
import default_parameter
import os

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
        print(f"{i}. {ff_dict[i][1]}   ({ff_dict[i][0][:-3]})")
    
    while True:
        try:
            choice = int(input("Enter your choice: "))
            if choice < 1 or choice > len(ff_dict)-1:
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
    
ff = select_ff()
print(apply_ff_patch(ff[0][:-3], default_parameter.nvt))
