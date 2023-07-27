# # make_parameter_sweep.py (062023) This is what we're going to use to essentially create fabrication plans Here,
# we're going to define an array of devices that we want to create in a sweep, and then feed the parameters into
# msg_cavity(cp). This way, we can create sets of devices that we have full simulated knowledge over,
# so we can actually trace trends.
import numpy as np
from cavity_functions import make_cavity, simulate_cavity, gds_cavity
from phidl import Device
import phidl.geometry as pg
from phidl import quickplot as qp
import pandas as pd
from dev_layout import dev_map

# There are two mechanisms for creating these sweeps, via brute force or via sensitivity studies. Brute force is that
# we simply define a range that we want to sweep over, for example cX \in (.5,.6,numValsX), a \in (.28,.3,numValsA).

# Sensitivity is way cooler I think, but requires many more things. We essentially start with an optimized cavity
# design, and then calculate its sensitivity to certain parameters (i.e. how much deviation can it have before Q <
# Q_thresh). From this, we can make parameter sweeps based on our known fabrication variances

# For either method, we are going to set up a param_sweep.csv so that we can track the designed parameter values +
# simulated properties for each device that we make. Later on, we can take this .csv and use it as reference for data
# analysis/fabrication offsets (akin to PPLN-type analysis)

# Define the cavity parameter list

cp = {'a': 0.2965 * 1e-6,
      'w0': .6*1e-6,
      'h0': .15*1e-6,
      'cX': 0.701,
      'cY': 1.7,
      'dA': 0.176364,
      'dY': 0,
      'dX': 0.176364,
      'dA_tap': 0,
      'dY_tap': .9,
      'dX_tap': .2,
      'nleft': 10,
      'nright': 3,
      'ndef': 6,
      'ntaper': 3,
      'hole_type': 'rib',
      'dfdA': .0684, # THz / nm
      'shift': -1,
      'target frequency': 406.7e12,
      'cav defect type': 'cubic',
      'min dim': 5e-8,
      'plot_E': False
      }


layout = {
    'layout_index': 0,
    'aper_index': 3,
    'layer': 4,
    'directory': 'G:/Shared drives/SiV-OMC/gds',
    'filename': None
}

do_sim = True

if do_sim:
    layout['filename'] = '071823_ribparametersweep_a-cY_simmed'
    column_data = ["Device x position", "Device y position", "a", "cY", "Frequency", "Scattering Q", "Waveguide Q",
                   "Mode volume"]
else:
    layout['filename'] = '071823_ribparametersweep_a-cY_notSimmed'
    column_data = ["Device x position", "Device y position", "a", "cX"]

dat_loc = layout['directory'] + '/' + layout['filename']
gds_name = dat_loc+'.gds'
df_name = dat_loc+'.csv'
E, pos_list = dev_map(layout)
deviationA = 0.05  # This is a 5% deviation on either side
deviationY = 0.05  # This a
numA = 10
numY = 10
aList = np.linspace(cp['a'] * (1 - deviationA), cp['a'] * (1 + deviationA), numA)
cYList = np.linspace(cp['cY'] * (1 - deviationY), cp['cY'] * (1 + deviationY), numY)
df = pd.DataFrame(data=np.empty((len(pos_list), len(column_data))), columns=column_data)
z = 0
while z < len(pos_list)-1:
    for a in aList:
        cp['source frequency'] = cp['target frequency'] + cp['dfdA'] * (cp['a'] - a)
        for cY in cYList:
            print(z)
            (x, y) = pos_list[z]
            cp['a'] = a
            cp['cY'] = cY
            cavity = make_cavity(cp)
            device = gds_cavity(cavity, cp).move(destination=(x, y)).rotate(45, center=(x, y))
            E.add(device)
            E.write_gds(gds_name)
            if do_sim:
                cav_results = simulate_cavity(cavity, cp)
                freq = cav_results['freq']
                qi = cav_results['qscat']
                qe = cav_results['qwvg']
                vmode = cav_results['vmode']
                df.loc[z, column_data] = [x, y, a*1e6, cY, freq, qi, qe, vmode]
                df.to_csv(df_name)
                cp['source frequency'] = freq
            else:
                df.loc[z, column_data] = [x, y, a*1e6, cY]
                df.to_csv(df_name)

            if z > len(pos_list)-2:
                break
            else:
                z += 1
        if z > len(pos_list)-2:
            break

print(df)
