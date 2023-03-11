# Basic structure
# Parameters: holeType ,a, hx, hy, beamW, beamH
import time
import numpy as np
from matplotlib import pyplot as plt
from wvgsolver import UnitCell
from wvgsolver.engine import LumericalEngine
from holeLibrary import build_hole

FDTDLoc = 'C:/Program Files/Lumerical/v221/'  # FDTDLoc, I will show you where this is later but this is how you can do it local

# The target resonance frequency, in Hz
target_frequency = 406.7e12

# Use level 3 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=4, hide=False, lumerical_path=FDTDLoc)

def band_structure(hole_params):
    hole, cell_size = build_hole(hole_params)
    cell = UnitCell(structures=hole, size=cell_size, engine=engine)
    start = time.time()
    # r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 12), freqs=(0.25e15, 0.7e15, 100000),
    #                    dipole_region=Vec3(0.8, 0, 0), window_pos = 0)
    r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 4), freqs=(0.15e15, 0.75e15, 300000))
    # # # Plot the bandstructure
    r1.show()
    end = time.time()
    print(end - start)


def band_gap(hole_params):
    hole, cell_size = build_hole(hole_params)
    cell = UnitCell(structures=hole, size=cell_size, engine=engine)

    start = time.time()
    r2 = cell.simulate("bandgap", freqs=(0.15e15, 0.75e15, 100000))
    end = time.time()

    diel_freq = r2[0]
    air_freq = r2[1]
    bg = air_freq - diel_freq
    mg = (diel_freq + air_freq) / 2
    bg_mg_rat = bg / mg

    delta_k = .5-hole_params['a']*1e-6*mg/2.998e8

    print('--------')
    print('Elapsed time: %f seconds' % (end - start))
    print('Lower band edge frequency: %f THz' % (diel_freq / 1e12))
    print("Bandgap ratio: %f percent" % (bg_mg_rat * 100))
    print("Midgap frequency: %f THz" % (mg / 1e12))
    print("Delta k: %f " % delta_k)
    print('\n')

    return diel_freq, air_freq, mg, bg_mg_rat, delta_k


aMir = .2942
hx = .625*aMir
hy = 1.76*aMir
w0 = 2.01*aMir
h0 = .503*aMir
a_list = np.linspace(aMir, (1-.173)*aMir, 5)

hole_params = {
        'a': aMir,
        'hx': hx,
        'hy': hy,
        'beam_height': h0,
        'beam_width': w0,
        'ref_index': 2.4028,
        'hole_type': 'tri_fil'
    }

diel_freq_list = a_list
air_freq_list = a_list

band_structure(hole_params)
for i in range(len(a_list)):
    hole_params['a'] = a_list[i]
    hole_params['hx'] = a_list[i]*.625
    diel_freq, air_freq, mg, bg_mg_rat,delta_k = band_gap(hole_params)
    diel_freq_list[i] = diel_freq
    air_freq_list[i] = air_freq


plt.plot(np.linspace(0, .173, 5), diel_freq_list/1e12)
plt.plot(np.linspace(0, .173, 5), air_freq_list/1e12)
plt.show()

