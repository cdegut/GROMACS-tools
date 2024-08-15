CHARMM36 = {
"constraints" : "h-bonds",
"cutoff-scheme" : "Verlet",
"vdwtype" : "cutoff",
"vdw-modifier" : "force-switch",
"rlist" : "1.2",
"rvdw" : "1.2",
"rvdw-switch" : "1.0",
"coulombtype" : "PME",
"rcoulomb" : "1.2",
"DispCorr" : "no"
}

amber99sbildn = {
"vdw-modifier":  "Potential-shift" ,
"DispCorr" : "EnerPress" }

patch_list = {'CHARMM36': CHARMM36, "amber99sb-ildn" : amber99sbildn  }

