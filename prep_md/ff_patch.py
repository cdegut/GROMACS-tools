CHARMM36 = {
"constraints" : ["h-bonds",],
"cutoff-scheme" :[ "Verlet",],
"vdwtype" : ["cutoff"],
"vdw-modifier" : ["force-switch"],
"rlist" : ["1.2"],
"rvdw" : ["1.2"],
"rvdw-switch" : ["1.0"],
"coulombtype" : ["PME"],
"rcoulomb" : ["1.2"],
"DispCorr" : ["no"]
}

amber99sbildn = {
"constraints"           :  ["h-bonds",],
"vdw-modifier"          :  ["Potential-shift-Verlet"] ,
"coulombtype"           :  ['PME','Particle Mesh Ewald for long-range electrostatics'],
"DispCorr"              :  ["EnerPres"],
"fourierspacing"        :  ['0.125','grid spacing for FFT'] 
}

# from https://gromacs.bioexcel.eu/t/protein-ligand-complex-amber-production-question/1247
amber99sbdisp = {
"constraints"           :  ["h-bonds",],
"vdw-modifier"          :  ["Potential-shift-Verlet"] ,
"coulombtype"           :  ['PME','Particle Mesh Ewald for long-range electrostatics'],
"DispCorr"              :  ["EnerPres"],
"fourierspacing"        :  ['0.125','grid spacing for FFT']
}

#name is folder name without the .ff
patch_list = {
    'CHARMM36': CHARMM36, 
    "amber99sb-ildn" : amber99sbildn, 
    "a99SBdisp" : amber99sbdisp }