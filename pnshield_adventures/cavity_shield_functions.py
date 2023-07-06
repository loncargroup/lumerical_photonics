import phidl
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, \
    PolygonStructure, MeshRegion
from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.engine import LumericalEngine
from wvgsolver.utils import BBox
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D
import numpy as np
import phidl.geometry as pg
import phidl.path as pp
from phidl import Path, CrossSection, Device

import os


def build_hole(hole_params):
    a = hole_params['a']
    h0 = hole_params['beam_height']
    w0 = hole_params['beam_width']
    hx = hole_params['hx']
    hy = hole_params['hy']
    hole_type = hole_params['hole_type']
    cell_pos = hole_params['cell_pos']

    if hole_type == 'ellipse':
        # offset the hole to respect the way we define the relevant lattice constant
        cell_hole = CylinderStructure(cell_pos, h0, hx / 2, DielectricMaterial(1, order=1, color="red"),
                                      radius2=hy / 2)
        return cell_hole


def linear_taper(i, ndef, maxdef):
    # We're going to setup point slope for this. (x2,y2) = (ndef,maxdef), (x1,y1) = (0,0)
    return 1 - maxdef * i / ndef


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

    a = pcc_params['a']
    hx_scale = pcc_params['cX']
    hy_scale = pcc_params['cY']

    # Make the left mirror region
    left_list = np.ones(nleft)
    lat_list = left_list * a
    hx_list = lat_list * hx_scale
    hy_list = lat_list * hy_scale

    # Make the cavity region
    if pcc_params['cav defect type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['cav defect type'] == 'quadratic':  ### NOTE THAT QUADRATIC IS NOT RIGHT Now
        defect = quadratic_defect

    defect_index = range(-ndef, ndef + 1)
    lat_def = pcc_params['dA']
    hx_def = pcc_params['dX']
    hy_def = pcc_params['dY']
    lat_list = np.append(lat_list, a * np.array([defect(i, ndef, lat_def) for i in defect_index]))
    hx_list = np.append(hx_list, a * hx_scale * np.array([defect(i, ndef, hx_def) for i in defect_index]))
    hy_list = np.append(hy_list, a * hy_scale * np.array([defect(i, ndef, hy_def) for i in defect_index]))

    # Make the right mirror
    right_list = np.ones(nright)
    lat_list = np.append(lat_list, a * right_list)
    hx_list = np.append(hx_list, a * hx_scale * right_list)
    hy_list = np.append(hy_list, a * hy_scale * right_list)

    w0 = pcc_params['w0']
    aMir = pcc_params['a']
    aDef = aMir * (1 - lat_def)
    hyMir = aMir * hy_scale
    hyDef = hyMir * (1 - hy_def)
    hxMir = aMir * hx_scale
    hxDef = hxMir * (1 - hx_def)
    min_dim = pcc_params['min dim']

    check_vals = [w0 - hyMir, w0 - hyDef, aMir - hxMir, aDef - hxDef]
    min_feature = min(check_vals)
    if min_feature < min_dim:
        return [None, None, None]
    else:
        return [lat_list, hx_list, hy_list]


def make_cavity(cp):
    hide = False
    FDTDLoc = 'C:/Program Files/Lumerical/v221/'  # FDTDLoc, I will show you where this is later but this is how you can do it local
    FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
    FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')
    [lat_list, hx_list, hy_list] = generate_cavity_and_check(cp)

    if lat_list is None:
        return None

    h0 = cp['h0']
    w0 = cp['w0']
    do_sc = cp['do_sc']

    if do_sc:
        sc_pos = Vec3(0, cp['gap'] + w0, 0)
        # [lat_list_sc, hx_list_sc, hy_list_sc,hole_list_sc] = generate_sidecouple_wvg(cp)

    beam_length = 30 * 1e-6
    log_name = "E:/User Data/michael/lum_results"
    engine = LumericalEngine(mesh_accuracy=5, hide=hide, lumerical_path=FDTDLoc, working_path=log_name + "/fsps")
    ref_index = 2.4028
    unit_cells = []
    for i in range(len(lat_list)):
        # Take each unit cell and add it to our model. First update the unit cell parameters according to our cell list
        a = lat_list[i]
        hx = hx_list[i]
        hy = hy_list[i]

        wg_size = Vec3(a, w0, h0)
        cell_box = BoxStructure(Vec3(0), wg_size, DielectricMaterial(ref_index, order=2, color="blue"))

        hole_params = {
            'a': a,
            'hx': hx,
            'hy': hy,
            'beam_width': w0,
            'beam_height': h0,
            'hole_type': cp['hole_type'],
            'shift': cp['shift'],
            'cell_pos': Vec3(0)
        }
        hole = build_hole(hole_params)
        cell_size = Vec3(a, w0, h0)
        if do_sc:
            cell_size = Vec3(a, 2 * w0 + cp['gap'], h0)
            sc_cell_box = BoxStructure(sc_pos, wg_size,
                                       DielectricMaterial(ref_index, order=2, color="blue"))
            if False:
                a_sc = lat_list_sc[i]
                hx_sc = hx_list_sc[i]
                hy_sc = hy_list_sc[i]
                sc_params = {
                    'a': a_sc,
                    'hx': hx_sc,
                    'hy': hy_sc,
                    'beam_width': w0,
                    'beam_height': h0,
                    'hole_type': cp['hole_type'],
                    'shift': cp['shift'],
                    'cell_pos': sc_pos
                }
                sc_hole = build_hole(sc_params)
                unit_cells += [
                    UnitCell(structures=[hole, sc_hole, cell_box, sc_cell_box], size=cell_size, engine=engine)]
            else:
                unit_cells += [
                    UnitCell(structures=[hole, cell_box, sc_cell_box], size=cell_size, engine=engine)]
        else:
            unit_cells += [
                UnitCell(structures=[hole, cell_box], size=cell_size, engine=engine)]

    if do_sc:
        structs = [BoxStructure(Vec3(0), Vec3(beam_length, w0, h0),
                                DielectricMaterial(2.4028, order=2, color="red")),
                   BoxStructure(sc_pos, Vec3(beam_length, w0, h0),
                                DielectricMaterial(2.4028, order=2, color="red"))
                   ]
    else:
        structs = [BoxStructure(Vec3(0), Vec3(beam_length, w0, h0),
                                DielectricMaterial(2.4028, order=2, color="red"))
                   ]

    cavity = Cavity1D(
        unit_cells=unit_cells,
        structures=structs,
        engine=engine,
        center_shift=cp['a']*(1-cp['dA'])
    )
    return cavity


def simulate_cavity(cavity, cp):
    source_frequency = cp['source frequency']
    man_mesh = MeshRegion(BBox(Vec3(0), Vec3(12e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)
    # There are a few different simulation conditions, let's keep this one but see if we need to change later
    r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=60e-15, analyze_fspan=10e12,
                         analyze_time=600e-15, mesh_regions=[man_mesh], TEonly=False, sim_size=Vec3(1.25, 3, 10))
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions=[man_mesh],
    #                     sim_size=Vec3(2, 3, 8))

    # r1 = cavity.simulate("resonance", target_freq=source_frequency)
    if cp['plot_E']:
        r1["xyprofile"].show()
        r1["yzprofile"].show()

    c = 2.99792458e17  # This is the speed of light in nm/s. As such, in order to get wavelength, it's just dividing c/freq
    qleft = r1["qxmin"]
    qright = r1["qxmax"]
    qy = 1 / (2 / r1["qymax"])
    qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
    freq = r1["freq"]
    vmode = r1["vmode"]

    qscat = 1 / (1 / qy + 1 / qz)

    wavelen_pen = np.exp(-((source_frequency - freq) / 4e12) ** 2)
    rerun_thresh = 0.90

    if wavelen_pen < rerun_thresh:
        # shift source frequency to cavity resonance and rerun simulation.
        # (this should help avoid non cavities with artificially low mode volumes)
        print('Source frequency =' + str(source_frequency))
        print('Simmed frequency =' + str(freq))
        print('Wavelen_penalty =' + str(wavelen_pen))
        print('---n')
        cp['source frequency'] = freq
        witness_rerun = simulate_cavity(cavity, cp)
        return witness_rerun

    cavity_results = {
        # Here are the values that we'll optimize over
        'qy': qy,
        'qz': qz,
        'qright': qright,
        'qleft': qleft,
        'vmode': vmode,
        # This is to keep track of wavelengths
        'freq': freq,
        'wl': c / freq,
    }

    print('---')
    return cavity_results
