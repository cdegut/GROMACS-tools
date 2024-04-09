import numpy as np
import pydssp
from unittest.mock import mock_open, patch
import io
import warnings
import MDAnalysis as mda
import matplotlib.pyplot as plt
from typing import Tuple, Union

import os

#MARK: test
##################################
##### DSSP #####
def dssp(protein):
    fake_file = io.StringIO()

    with patch("builtins.open", new_callable=mock_open, create=True) as mock_open_func:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            protein.write(f"tmp.pdb", file_format="pdb")
            # Get the file handle from the mock
            handle = mock_open_func.return_value
            # Extract the written content
            written_content = "".join(call[1][0] for call in handle.write.mock_calls)
            # Write the content to the fake file
            fake_file.write(written_content)

    # Reset the file pointer to the beginning
    fake_file.seek(0)

    # Read data from the fake file
    data = fake_file.read()
    data = data.split("\n")

    backbone = []
    for line in data:
        try:
            if line[13:15] in ["N ","CA", "C ", "O "]:
                backbone.append(line)
        except:
            pass

    add_last = True 
    last = backbone[-1][13:15]
    if last == "O ":
        add_last = False
    elif last == "C ":
        backbone = backbone[:-3]
    elif last == "CA ":
        backbone = backbone[:-2]
    elif last == "N ":
        backbone = backbone[:-1]

    backbone = "\n".join(backbone)
    coord = np.array(pydssp.read_pdbtext(backbone))
    dsspline = ''.join(pydssp.assign(coord))
    
    if add_last:
        dsspline = f"{dsspline}-"
    return dsspline


#MARK: plot_dssp

def calc_dssp_matrix(atomistic_system: mda.Universe, step, start = None, stop = None ):

    protein = atomistic_system.select_atoms("protein")

    if start is None:
        start = 0
    else:
        start = start * 100
    if stop is None:
        stop = len(atomistic_system.trajectory)
    else:
        stop = stop*100


    helix_matrix = np.zeros((int((stop - start)/step)+1, int(len(atomistic_system.residues))), np.float16)
    beta_matrix = np.zeros((int((stop - start)/step)+1, int(len(atomistic_system.residues))), np.float16)

    for  i, ts in enumerate(range(start, stop, step)):
        atomistic_system.trajectory[ts]
        line = dssp(protein)
        helix_value = line.replace("-", "0").replace("H","1").replace("E", "0")
        alpha_as_array = np.array([float(x) for x in helix_value])
        beta_value = line.replace("-", "0").replace("H","0").replace("E", "1")
        beta_as_array = np.array([float(x) for x in beta_value])
        helix_matrix[i] = alpha_as_array
        beta_matrix[i] = beta_as_array
    return helix_matrix, beta_matrix


def plot_dssp_average(atomistic_system: mda.Universe, helix_matrix, beta_matrix) -> Tuple[plt.Figure, plt.Axes]:
    fig= plt.figure(constrained_layout=True)
    fig.set(figwidth=8, figheight=3)
    bar_plot = fig.add_subplot()
    x = np.arange(1,len(atomistic_system.residues)+1, 1)

    bar_plot.bar(x, helix_matrix.mean(0), color="crimson", alpha = 0.7, label = "Alpha conf")
    bar_plot.bar(x, beta_matrix.mean(0), color="yellowgreen", alpha = 0.7, label = "Beta conf")

    bar_plot.set_xlim(1, len(atomistic_system.residues))
    start, end = bar_plot.get_xlim()
    bar_plot.xaxis.set_ticks(np.arange(start-1, end, 5))

    bar_plot.set_xlabel('Residues')
    bar_plot.set_ylabel('Fraction of time in conformation')

    return fig, bar_plot
