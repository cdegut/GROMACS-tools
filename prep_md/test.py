import subprocess
selected_file = "SasY_NMR_double.pdb"
commands = [f'gmx pdb2gmx -ignh -f {selected_file} -o starting_structure.gro'
]

for command in commands:
    cmd = subprocess.run(command, shell=True)
    if cmd.returncode != 0:
        break

