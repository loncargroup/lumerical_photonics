##Make a set of functions
# fitness(hole_params)
# build_hole(hole_params)

import numpy as np
import time
from scipy.optimize import minimize
from holeLibrary import build_hole
from wvgsolver import UnitCell
from wvgsolver.engine import LumericalEngine

FDTDLoc = 'C:/Program Files/Lumerical/v221/'
iter_count = 0
target_frequency = 406.7e12
hide= True
print_text = True
rerun_thresh = .94


def check_hole(hole_params):
    a = hole_params['a'] * 1e-6
    h0 = hole_params['beam_height']*1e-6
    w0 = hole_params['beam_width'] * 1e-6
    hx = hole_params['hx']*1e-6
    hy = hole_params['hy']*1e-6
    min_dim = hole_params['min_dim']
    if hole_params['hole_type'] == 'ellipse':
        if (a-hx) < min_dim or hx < min_dim:
            return True
        if (w0 - hy) < min_dim or hy < min_dim:
            return True
        return False
    elif hole_params['hole_type'] == 'ribAir' or hole_params['hole_type'] == 'ribDiel':
        if (a-hx) < min_dim or hx < min_dim:
            return True
        if (w0 - hy) < min_dim or hy < min_dim:
            print('Fabrication incompatible')
            return True
        return False


def band_gap(hole_params, engine):
    hole, cell_size = build_hole(hole_params)
    cell = UnitCell(structures=hole, size=cell_size, engine=engine)

    start = time.time()
    r2 = cell.simulate("bandgap", freqs=(0.15e15, 0.5e15, 100000))
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


def fitness(hole_parameters):
    global hide
    global iter_count
    global print_text
    global rerun_thresh
    global lat_const

    # aMir = .315
    # hx = .7413 * aMir
    # hy = 1.704 * aMir
    # w0 = 1.865 * aMir
    # h0 = .53 * aMir
    # a_list = np.linspace(aMir, .9 * aMir, 5)

    hole_params = {
        'a': lat_const,
        'hx': hole_parameters[0]*lat_const,
        'hy': hole_parameters[1]*lat_const,
        'beam_width': hole_parameters[2]*lat_const,
        'beam_height': hole_parameters[3] * lat_const,
        'ref_index': 2.4028,
        'hole_type': 'ribAir',
        'min_dim': 8e-8
    }

    fab_flag = check_hole(hole_params)

    if fab_flag:
        witness = 0
        print('Fabrication intolerant')
    else:
        iter_count += 1
        engine = LumericalEngine(mesh_accuracy=4, hide=hide, lumerical_path=FDTDLoc)
        diel_freq, air_freq, mg, bg_mg_rat, delta_k = band_gap(hole_params, engine)

        wavelength_pen = np.exp(-((target_frequency - mg)/target_frequency)**2)
        witness = -1*bg_mg_rat*delta_k
        if witness < -.007:
            lat_const = mg / target_frequency * lat_const

        # if (wavelength_pen < rerun_thresh) and witness < -.1:
        #     #Shift lattice constant to
        #     witness_rerun = fitness(hole_parameters)
        #     print("rerun. Fitness when source is recentered:", witness_rerun)
        #     return witness_rerun

        with open(log_name, "ab") as f:
            f.write(b"\n")
            step_info = np.append(hole_parameters, np.array([lat_const, witness, wavelength_pen, mg, bg_mg_rat, delta_k]))
            np.savetxt(f, step_info.reshape(1, step_info.shape[0]), fmt='%.6f')

    return witness


p0 = np.array([.53, 1.5, 1.8, 1])
# cX, cY, cW, cH
bounds = ((.5, .7), (.7, 2.2), (1, 2.4), (.4, 1.3))
lat_const = .33
log_name = f"C:/Users/LoncarLab/Documents/GitHub/wvgsolver/examples/ribbed_cavity_optimization" \
           f"/optimization_022723_rib_unit/optimal_rib_unit_80nm_limit_00.txt"

with open(log_name, "ab") as f:
    f.write(b"cX    cY  cW cH     fitness    wavelen_pen   mg  bg")

# popt = minimize(fitness, p0, bounds = bounds, method='Nelder-Mead')
# print(popt)
witness = fitness(p0)
