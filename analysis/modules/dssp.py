import numpy as np
import pydssp
from unittest.mock import mock_open, patch
import io
import warnings
import os

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