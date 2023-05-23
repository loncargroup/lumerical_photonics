## adiabatic_tapering.py (052223)
# The major goal of this code is to integrate the adiabatic tapering that will fix the scattering issues that we were seeing
# with recent devices.
# We will have 4 major numbers that define our unit cells. n_left, n_def, n_right_, n_taper.

import numpy as np
import matplotlib.pyplot as plt
import os

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine
from holeLibrary import build_hole

iter_count = 0

rerun_thresh = .94
target_frequency = 406.7e12
source_frequency = 406.7e12
hide = False

FDTDLoc = 'C:/Program Files/Lumerical/v221/'   # FDTDLoc, I will show you where this is later but this is how you can do it local
FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')


def cubic_defect(i, ndef, maxdef):
    x = np.abs(i)
    return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)


def quadratic_defect(i, ndef, maxdef):
    x = np.abs(i)
    return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)


def generate_cavity_and_check(pcc_params):
    ndef = pcc_params['ndef']
    nleft = pcc_params['nleft']
    nright = pcc_params['nright']
    ntaper = pcc_params['ntaper']

    a = pcc_params['a']
    w0 = pcc_params['cW'] * a
    hx_scale = pcc_params['cX']
    hy_scale = pcc_params['cY']

    # Make the left mirror region
    left_list = np.ones(nleft)
    lat_list = left_list * a
    hx_list = lat_list * hx_scale
    hy_list = lat_list * hy_scale

    # Make the cavity region
    if pcc_params['cav_defect_type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['cav_defect_type'] == 'quadratic':  ### NOTE THAT QUADRATIC IS NOT RIGHT Now
        defect = quadratic_defect

    defect_index = range(-ndef, ndef)
    lat_def = pcc_params['dA']
    hx_def = pcc_params['dX']
    hy_def = pcc_params['dY']
    lat_list = np.append(lat_list, a * np.array([defect(i, ndef, lat_def) for i in defect_index]))
    hx_list = np.append(hx_list, a * hx_scale * np.array([defect(i, ndef, hx_def) for i in defect_index]))
    hy_list = np.append(hy_list, a * hy_scale * np.array([defect(i, ndef, hy_def) for i in defect_index]))

    #Make the right mirror
    right_list = np.ones(nright)
    lat_list = np.append(lat_list, a * right_list)
    hx_list = np.append(hx_list, a * hx_scale * right_list)
    hy_list = np.append(hy_list, a * hy_scale * right_list)

    #Make the taper region
    if pcc_params['tap_defect_type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['tap_defect_type'] == 'quadratic':  ### NOTE THAT QUADRATIC IS NOT RIGHT Now
        defect = quadratic_defect

    tap_index = range(ntaper, 1, -1)
    lat_def_tap = pcc_params['dA_tap']
    hx_def_tap = pcc_params['dX_tap']
    hy_def_tap = pcc_params['dY_tap']

    lat_list = np.append(lat_list, a * np.array([defect(i, ndef, lat_def_tap) for i in tap_index]))
    hx_list = np.append(hx_list, a * hx_scale * np.array([defect(i, ndef, hx_def_tap) for i in tap_index]))
    hy_list = np.append(hy_list, a * hy_scale * np.array([defect(i, ndef, hy_def_tap) for i in tap_index]))

    return [lat_list, hx_list, hy_list]


def build_cavity(cavity_params):
    global nleft
    global nright
    global ndef
    global ntaper
    global iter_count
    global hide
    iter_count += 1
    # Put this for now, but this might not be the way we make it
    pcc_params = {
        'a': lat_const,
        'cX': cX,
        'cY': cY,
        'cW': cW,
        'cH': cH,
        'dA': lat_def,
        'dY': 0,
        'dX': hx_def,
        'dA_tap': lat_def_tap,
        'dY_tap': hy_def_tap,
        'dX_tap': hx_def_tap,
        'nleft': nleft,
        'nright': nright,
        'ndef': ndef,
        'ntaper': ntaper,
        'ref_index': 2.4028,
        'hole_type': 'rib',
        'effective_index': 1.6,
        'resonance_wavelength': 0.737,
        'cav_defect_type': 'cubic',
        'tap_defect_type': 'cubic',
        'min_dim': 5e-8
    }

    beam_length = 30
    nleft = pcc_params['nleft']
    ndef = pcc_params['ndef']
    nright = pcc_params['nright']
    ntaper = pcc_params['ntaper']
    a_list, hx_list, hy_list = generate_cavity_and_check(pcc_params)
    a = pcc_params['a']
    h0 = pcc_params['cH'] * a
    w0 = pcc_params['cW'] * a

    unit_cells = []

    engine = LumericalEngine(mesh_accuracy=5, hide=hide, lumerical_path=FDTDLoc, working_path="./fsps")

    for (a, hx, hy) in zip(a_list, hx_list, hy_list):
        # Take each unit cell and add it to our model. First update the unit cell parameters according to our cell list

        hole_params = {
            'a': a,
            'hx': hx,
            'hy': hy,
            'beam_width': w0,
            'beam_height': h0,
            'ref_index': 2.4028,
            'hole_type': pcc_params['hole_type']
        }

        hole, cell_size = build_hole(hole_params)
        unit_cells += [UnitCell(structures=hole, size=cell_size, engine=engine)]

    cavity = Cavity1D(
        unit_cells=unit_cells,
        structures=[BoxStructure(Vec3(0), Vec3(beam_length * 1e-6, w0 * 1e-6, h0 * 1e-6),
                                 DielectricMaterial(2.4028, order=2, color="red"))],
        engine=engine,
    )

    cavity_name = np.array2string(cavity_params, prefix='', separator='_')
    cavity_name = cavity_name[1:]
    cavity_name = cavity_name[:-1]
    cavity_name = cavity_name.replace(" ", "")
    file_name = log_name[:-4] + "_" + cavity_name + "_" + str(iter_count)
    file_name = log_name[:-4] + "_rib_00_MidSim.obj"
    cavity.save(file_name)
    # Once we check plot_geom, uncomment this
    # plot_geom(cavity, file_name + "_geom.png", hide)

    return cavity, file_name, engine


def run_cavity(cavity_params):
    global iter_count
    global rerun_thresh
    global target_frequency
    global source_frequency
    global lat_const
    cavity, file_name, engine = build_cavity(cavity_params)
    flag = False # Have this for now, but change when we put fabrication tolerance back in
    if flag:
        print('Fabrication intolerant')
        return 0
    else:
        man_mesh = MeshRegion(BBox(Vec3(0), Vec3(12e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)

        r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=60e-15, analyze_fspan=10e12,
                             analyze_time=600e-15, mesh_regions=[man_mesh], sim_size=Vec3(1.25, 3, 8))
        # r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions=[man_mesh],
        #                     sim_size=Vec3(2, 3, 8))

        # r1 = cavity.simulate("resonance", target_freq=source_frequency)
        print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
            r1["freq"], r1["vmode"],
            1 / (1 / r1["qxmin"] + 1 / r1["qxmax"]),
            1 / (2 / r1["qymax"] + 1 / r1["qzmin"] + 1 / r1["qzmax"])
        ))

        qx = 1 / (1 / r1["qxmin"] + 1 / r1["qxmax"])
        qx1 = r1["qxmin"]
        qx2 = r1["qxmax"]
        qy = 1 / (2 / r1["qymax"])
        qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])

        print("Qx: %f " % qx)
        print("Qy: %f " % qy)
        print("Qz: %f " % qz)

        r1["xyprofile"].show()
        r1["yzprofile"].show()

        # Decide how we want to define qx, it's possible that we only want to define it as 'qxmax'
        # since that's the rightward flux. This would require a new definition of qscat

        cavity = Cavity1D(load_path=file_name, engine=engine)
        r2 = cavity.get_results("resonance")[-1]
        print(r2)
        print(r2['res']["xyprofile"].max_loc())
        print(r2['res']["yzprofile"].max_loc())
        r2["sess_res"].show()

        print(r2)

        qscat = 1 / ((1 / qy) + (1 / qz))
        qtot = 1 / (1 / qx + 1 / qy + 1 / qz)
        vmode = r1["vmode"]
        vmode_copy = vmode
        vmode = 1e6 if vmode < 0.48 else vmode
        qtot_max = 300000
        purcell = qtot / vmode if qtot < qtot_max else qtot_max / vmode
        F = r1["freq"]
        wavelen = (2.99e8 / F) * 1e9

        wavelen_pen = np.exp(-((target_frequency - F) / 4e12) ** 2)
        # qx_pen = (np.exp(-((target_qx - qx) / 120000) ** 2) + np.exp(-((target_qx - qx) / 60000) ** 2)
        #          + np.exp(-((target_qx - qx) / 30000) ** 2) + np.exp(-((target_qx - qx) / 2000) ** 2) / 4)
        qscat_max = 500000
        guidedness = qscat / qx if qscat < qscat_max else qscat_max / qx
        # witness = -1 * purcell * guidedness # * qx_pen
        witness = -1 * qscat / vmode if qscat < qscat_max else -1 * qscat_max / vmode

        # r1["xyprofile"].save(file_name+"_xy.png",title=f"Q = {qtot:.0f} \nQ_scat = {qscat:.04} Qx = {qx:.0f}\nV = {vmode_copy:.3f}")
        # r1["yzprofile"].save(file_name+"_yz.png",title=f"Q = {qtot:.0f} Q_scat = {qscat:.04}\n Qx1 = {qx1:.0f} Qx2 = {qx2:.0f}\nV = {vmode_copy:.3f} "+r"$\lambda$"+f" = {wavelen:.1f}")

        # second condition ensures that we only rerun once
        if ((wavelen_pen < rerun_thresh) and (source_frequency == target_frequency)):
            # shift source frequency to cavity resonance and rerun simulation.
            # (this should help avoid non cavities with artificially low mode volumes)
            source_frequency = F
            witness_rerun = run_cavity(cavity_params)
            print("rerun. Fitness when source is recentered:", witness_rerun)
            source_frequency = target_frequency
            return witness_rerun

        if witness < -.07:
            lat_const = F / target_frequency * lat_const

        with open(log_name, "ab") as f:
            f.write(b"\n")
            step_info = np.append(cavity_params, np.array(
                [witness, wavelen_pen, purcell, r1["qxmin"], r1["qxmax"], qscat, qtot, vmode, vmode_copy, F,
                 lat_const]))
            np.savetxt(f, step_info.reshape(1, step_info.shape[0]), fmt='%.6f')

        return witness


log_name = f"E:/User Data/michael/lum_results/rib_stuff/rib_stuff_00.txt"

lat_const = .287
cX = .625
cY = 1.76
cW = 2.01
cH = .503

lat_def = .1
hx_def = .1
hy_def = 0

lat_def_tap = .2
hx_def_tap = 0
hy_def_tap = .3

nleft = 7
nright = 2
ndef = 6
ntaper = 5

p0 = np.array([.173])
# dA, dX, dY, cX
bounds = ((.05, .25), (-.1, .1), (-.1, .1), (.5, .9))



with open(log_name, "ab") as f:
    f.write(
        b"dA    dX  dY  fitness wavelen_pen purcell qxmin   qxmax   qscat   qtot    vmode   vmode_copy  F lat_const")

# witness = fitness(p0)
run_cavity(p0)
