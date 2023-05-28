# sensitivity_study.py
# The goal of this is to figure out given an "optimal design", how much the Q changes as a function of deviation for a
# given parameter. We've written a nice run_cavity function that allows us to sort of 'hide' the lumerical, so now it's
# just a matter of doing a parameter sweep.

# Idea: take some deviation (call it \sigma_p) and vary some parameter over that range
# E.g. if sigma_x = 50 nm, then we'll take n_points over +-50 nm for the hx parameter and run simulations to measure
# freq, vmode, qe, qi

import os
import numpy as np
from matplotlib import pyplot as plt
from common_functions import run_cavity

hide = True
FDTDLoc = 'C:/Program Files/Lumerical/v221/'  # FDTDLoc, I will show you where this is later but this is how you can do it local
FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')

data_dir = f"E:/User Data/michael/lum_results/rib_stuff/"
log_name = data_dir + f"rib_stuff_00.txt"

cp = {'a': 0.2916,
      'cX': 0.67,
      'cY': 1.7,
      'cW': 2,
      'cH': .5,
      'dA': 0.176,
      'dY': 0,
      'dX': 0.176,
      'dA_tap': 0,
      'dY_tap': 0,
      'dX_tap': 0,
      'nleft': 10,
      'nright': 10,
      'ndef': 6,
      'ntaper': 0,
      'hole_type': 'rib',
      'shift': -1,
      'source frequency': 406.7e12}

sigma = .02  # Deviation in um
n_points = 6  # Number of points from hx to hx+-sigma_x
parameter = 'cX'  # Parameter that we're deviating

c_min = (cp['a'] * cp[parameter] - sigma) / cp['a']
c_max = (cp['a'] * cp[parameter] + sigma) / cp['a']
c_range_up = np.linspace(cp[parameter], c_max, n_points)
c_range_down = np.linspace(cp[parameter], c_min, n_points)

freq_lu = []
vmode_lu = []
qe_lu = []
qi_lu = []
print(c_range_up)
for val in c_range_up:
    cp[parameter] = val
    freq, vmode, qe, qi = run_cavity(cp, data_dir, FDTDLoc)
    freq_lu.append(freq)
    vmode_lu.append(vmode)
    qe_lu.append(qe)
    qi_lu.append(qi)
    cp['source frequency'] = freq
    print(freq)
    print(qe)
    print(qi)
    print('')

freq_ld = []
vmode_ld = []
qe_ld = []
qi_ld = []
cp['source frequency'] = 406.7e12
print(c_range_down)
for val in c_range_down:
    cp[parameter] = val
    freq, vmode, qe, qi = run_cavity(cp, data_dir, FDTDLoc)
    freq_ld.append(freq)
    vmode_ld.append(vmode)
    qe_ld.append(qe)
    qi_ld.append(qi)
    cp['source frequency'] = freq
    print(freq)
    print(qe)
    print(qi)
    print('')

sig_lu = np.linspace(0, sigma, n_points)
sig_l = np.linspace(-sigma, 0, n_points)
for x in sig_lu:
    sig_l = np.append(sig_l, x)

for x in qe_lu:
    qe_ld = np.append(qe_ld, x)

for x in qi_lu:
    qi_ld = np.append(qi_ld, x)

for x in vmode_lu:
    vmode_ld = np.append(vmode_ld, x)

plt.plot(sig_l, qi_ld)
plt.show()
