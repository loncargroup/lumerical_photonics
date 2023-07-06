"""This function will be used to run the optical simulation for our side-coupled optomechanics, let's finish this
docstring when we know what we're doing

Michael Haas 062523
"""

import numpy as np
from cavity_shield_functions import make_cavity, simulate_cavity
import phidl.geometry as pg
import pandas as pd
from phidl import quickplot as qp

"""
rows = []
with open('hole_info.csv', 'r') as csv_file:
    reader = csv.reader(csv_file)

    for row in reader:
        rows.append(row)

vals = []
with open('beam_info.csv', 'r') as csv_file:
    reader = csv.reader(csv_file)

    for row in reader:
        vals.append(row)

xpos_list = [float(i) for i in rows[0]]
lat_list = [float(i) for i in rows[1]]
hx_list = [float(i) for i in rows[2]]
hy_list = [float(i) for i in rows[3]]
w0 = float(vals[0][0])
h0 = float(vals[1][0])
hole_val = int(vals[2][0])

if hole_val == 0:
    hole_type = 'ellipse'

rerun_thresh = .94
target_frequency = 406.7e12
source_frequency = 406.7e12
hide = False
beam_length = 30

unit_cells = []

This is what was done for omc_wvgsolver.py, let's do it slightly cleaner using dictionaries. 

Let's actually totally change this, so that we are only making structures in python, never in COMSOL. 
This way, we can ensure that we're accurate in our structure definition

This means that we will still start in COMSOL, but we're only going to be using it to define our structure, then we'll 
do the heavy lifting in python. As such, let's just do all our definitions here.
"""

cp = {  # Hole parameters
    'a': 0.274 * 1e-6,
    'w0': .548 * 1e-6,
    'h0': .137 * 1e-6,
    'cX': 0.683,
    'cY': 1.55,
    'nleft': 6,
    'nright': 6,
    'hole_type': 'ellipse',
    'shift': -1,
    # Taper parameters
    'dA': 0.173,
    'dY': 0,
    'dX': 0.173,
    'cav defect type': 'cubic',
    'ndef': 6,
    # Shield parameters
    'aS': .3 * 1e-6,
    'cSw': .5,
    'source frequency': 406.7e12,
    # Side-couple parameters
    'do_sc': True,
    'gap': .25e-6,
    'min dim': .05 * 1e-6,
    'plot_E': True
}

cavity = make_cavity(cp)
cavity_results = simulate_cavity(cavity, cp)
print(cavity_results)
