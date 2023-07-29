# Set the CADET binary path_____________________________________________________

import platform
from pathlib import Path
from cadet import Cadet

cadet_bin_path = Path.home() / "codes"/ "cadet" / "cadet"/ "bin"

cadet_path = cadet_bin_path / "cadet-cli"
lwe_path = cadet_bin_path / "createLWE"

if cadet_path.exists() and lwe_path.exists():
    Cadet.cadet_path = cadet_path.as_posix()
elif cadet_path.exists() and not lwe_path.exists():
    print("CADET was found but createLWE.exe was not found. \
           Please make sure that none of the files have been moved.")
else:
    print("CADET could not be found. Please check the bin path")

# Standard imports______________________________________________________________

import os
import subprocess
import time

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from addict import Dict
import json

# Temporary files for simulation objects
import tempfile
tempfile.tempdir = os.path.join(Path.home())

# IPython
from IPython.core.display import display, HTML, clear_output
display(HTML("<style>.container { width:100% !important; }</style>"))

# CADETMatch____________________________________________________________________

from CADETMatch.jupyter import Match
