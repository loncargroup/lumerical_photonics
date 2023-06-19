## simulate_and_gds.py
'''
Authors: Katie Barajas, Michael Haas
Modified from code by Aaron Day and Chang Jin
'''

#Basic idea:
# 1. Take a cavity parameter array cp, simulate the cavity to get freq, vmode, qe, qi (run_cavity.py)
# 2. Parse into phidl using chang's code (SiC_gds_gen_v5.py, SiC_Overcoupled_gds.ipynb, SiC_cavity_array_code.ipynb)
# 3. Spit out gds of full cavity and a line with parameter values, frequency, qe, qi

import os
import numpy as np
from matplotlib import pyplot as plt
from common_functions import run_cavity


from IPython import get_ipython;
get_ipython().magic('reset -sf')
from phidl import Device
from phidl import quickplot as qp # Rename "quickplot()" to the easier "qp()"
import phidl.geometry as pg
import numpy as np


hide = True
FDTDLoc = 'C:/Program Files/Lumerical/v221/'
FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')

data_dir = f"E:/User Data/michael/lum_results/rib_stuff/"
log_name = data_dir + f"rib_simulate_and_gds.txt"

'''
simulate_and_gds()

Method: 
Take cavity parameter array, cp, simulate the cavity to get freq, vmode, qe, qi by calling run_cavity. 

Return full cavity and a line with parameter values, frequency, qe, qi 
'''
def simulate_and_gds(cp):
    # simulates cavity and returns results in variables freq, vmode, qe, qi
    freq, vmode, qe, qi = run_cavity(cp, data_dir, FDTDLoc)
    output = [freq, vmode, qe, qi]

    #run gds scripting code with phidl
    D = Device() # create blank device and add unit cell to it
    





    return output
