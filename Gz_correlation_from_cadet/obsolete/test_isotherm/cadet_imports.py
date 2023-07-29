# Set the CADET binary path_____________________________________________________

import platform
from pathlib import Path
from cadet import Cadet

cadet_bin_path = Path.home() / "codes"/ "cadet" / 'cadet' / 'bin'


if platform.system() == 'Windows':
    cadet_path = cadet_bin_path / "cadet-cli.exe"
    lwe_path = cadet_bin_path / "createLWE.exe"
else:
    cadet_path = cadet_bin_path / "cadet-cli"
    lwe_path = cadet_bin_path / "createLWE"

if cadet_path.exists() and lwe_path.exists():
    Cadet.cadet_path = cadet_path.as_posix()
elif cadet_path.exists() and not lwe_path.exists():
    print("CADET was found but createLWE.exe was not found. \
           Please make sure that none of the files have been moved.")
else:
    print("CADET could not be found. Please check the bin path")

# cadet_path = '/home/chase/codes/cadet/cadet/bin/cadet-cli'
# Cadet.cadet_path = cadet_path
# lwe_path = '/home/chase/codes/cadet/cadet/bin/createLWE'

# Standard imports______________________________________________________________

import os
import subprocess
import time

import numpy as np
import scipy
import pandas as pd # I added as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from addict import Dict
import json

# Temporary files for simulation objects
import tempfile
tempfile.tempdir = os.path.join(Path.home())

# IPython

from IPython.display import Image
from ipywidgets import interact, interactive
import ipywidgets as widgets

import copy

# CADETMatch____________________________________________________________________

from CADETMatch.jupyter import Match
