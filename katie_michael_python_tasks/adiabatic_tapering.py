## adiabatic_tapering.py (052223)
# The major goal of this code is to integrate the adiabatic tapering that will fix the scattering issues that we were seeing
# with recent devices.
# We will have 4 major numbers that define our unit cells. n_left, n_def, n_right_, n_taper.

import os

from common_functions import run_cavity

hide = False
FDTDLoc = 'C:/Program Files/Lumerical/v221/'  # FDTDLoc, I will show you where this is later but this is how you can do it local
FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')

data_dir = f"E:/User Data/michael/lum_results/rib_stuff/"
log_name = data_dir + f"rib_stuff_00.txt"

cp = {'a': 0.296,
      'cX': 0.701,
      'cY': 1.7,
      'cW': 2,
      'cH': .5,
      'dA': 0.176,
      'dY': 0,
      'dX': 0.176,
      'dA_tap': 0,
      'dY_tap': 0,
      'dX_tap': 0,
      'nleft': 6,
      'nright': 6,
      'ndef': 6,
      'ntaper': 0,
      'hole_type': 'rib',
      'shift': -1,
      'source frequency': 406.7e12}


freq, vmode, qe, qi = run_cavity(cp, data_dir, FDTDLoc)
print(freq)
print(vmode)
print(qe)
print(qi)
