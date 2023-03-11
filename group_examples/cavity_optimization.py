# Ribbed cavity optimization (Michael Haas 022623)
# Use wvgsolver interface to optimize ribbed cavity designs for AWS project

# Outline:
#  1. Make a set of functions: fitness, build_cavity, check_cavity, plot_geom
#  2. Save results to a file as well as nice plots of the taper functions (modify Erik's code function for this)

#Imports
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import os

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine

from holeLibrary import build_hole

FDTDLoc = 'C:/Program Files/Lumerical/v221/'   # FDTDLoc, I will show you where this is later but this is how you can do it local
FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')

iter_count = 0
nleft = 7
nright = 7
ndef = 6
rerun_thresh = .94
target_frequency = 406.7e12
source_frequency = 406.7e12
hide = False


def plot_geom(cavity, file_name, hide):
    plt.close()
    unit_cells = cavity.get_unit_cells()

    a_s = np.zeros(len(unit_cells))
    hx = np.zeros(len(unit_cells))
    hy = np.zeros(len(unit_cells))
    for i, unit_cell in enumerate(unit_cells):
        size = unit_cell.get_size()
        structures = unit_cell.get_structures()
        hx[i] = structures[1].radius*2
        hy[i] = structures[1].radius2*2
        a_s[i] = size.x

    x_coords = np.copy(a_s)
    x_coords = x_coords[:-1]
    x_coords = np.insert(x_coords,0,0)
    x_coord = 0
    for i, a in enumerate(x_coords):
        x_coord += a
        x_coords[i] = x_coord

    a_coords = np.array([(x+x_coords[i+1])/2 for i, x in enumerate(x_coords[:-1])])

    plt.close()

    plt.vlines(a_coords, 0.18, 1.01, color='gray', alpha=0.3)
    plt.hlines(1, np.min(x_coords), np.max(x_coords), colors='red', linestyle='dashed')
    plt.plot(a_coords, a_s[:-1]/a_s[5], 'ko-',label='a')
    plt.plot(x_coords, hx/a_s[5], 'go-', label='hx')
    plt.plot(x_coords, hy/a_s[5], 'bo-', label='hy')

    plt.ylabel(r"$\mathregular{a}_{\mathregular{nominal}}$",fontsize=14)
    plt.xlabel(r"Unit Cell Position ($\mu$m)", fontsize=12)

    plt.legend(loc=8)

    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.tight_layout()
    plt.savefig(file_name)
    if not hide:
        plt.show()
    plt.close()


def cubic_defect(i, nleft, nright, ndef, maxdef):
    if nleft < i < nleft + 2 * ndef - 1:
        x = np.abs(i - nleft - ndef)
        return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)
    else:
        return 1


def quadratic_defect(i, nleft, nright, ndef, maxdef):
    if nleft < i < nleft + 2 * ndef - 1:
        x = np.abs(i - nleft - ndef)
        return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)
    else:
        return 1


def generate_cavity_and_check(pcc_params):
    ndef = pcc_params['ndef']
    nleft = pcc_params['nleft']
    nright = pcc_params['nright']
    ntot = nleft + 2 * ndef + nright
    hole_index = np.arange(ntot)

    if pcc_params['defect_type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['defect_type'] == 'quadratic':
        defect = quadratic_defect

    a = pcc_params['a']
    lat_def = pcc_params['dA']
    lat_def_list = np.array([defect(i, nleft, nright, ndef, lat_def) for i in hole_index])
    lat_list = a * lat_def_list

    hx_scale = pcc_params['cX']
    hx_def = pcc_params['dX']
    hx_def_list = np.array([defect(i, nleft, nright, ndef, hx_def) for i in hole_index])
    hx_list = hx_scale * a * hx_def_list

    hy_scale = pcc_params['cY']
    hy_def = pcc_params['dY']
    hy_def_list = np.array([defect(i, nleft, nright, ndef, hy_def) for i in hole_index])
    hy_list = hy_scale * a * hy_def_list

    # Check fabrication feasibility
    # Check mirror cell
    w0 = pcc_params['cW']*a
    a_mirr = a
    hx_mirr = hx_list[0]
    hy_mirr = hy_list[0]

    a_def = lat_list[nleft + ndef]
    hx_def = hx_list[nleft + ndef]
    hy_def = hy_list[nleft + ndef]
    min_dim = pcc_params['min_dim']
    flag = False
    if a_mirr-hx_mirr < min_dim or hx_mirr < min_dim or (w0 - hy_mirr) < min_dim or hy_mirr < min_dim \
            or (a_def-hx_def) < min_dim or hx_def < min_dim or (w0 - hy_def) < min_dim or hy_def < min_dim:
        # If any of the cavity hole dimensions drop below a minimum fabrication tolerance, throw it out
        flag = True

    return [lat_list, hx_list, hy_list, flag]


def build_cavity(cavity_params):
    global nleft
    global nright
    global ndef
    global iter_count
    global hide
    iter_count += 1
    #Put this for now, but this might not be the way we make it
    lat_def, hx_scale, hy_scale = cavity_params
    hx_def = lat_def
    pcc_params = {
        'a': lat_const,
        'cX': hx_scale,
        'cY': hy_scale,
        'cW': 2,
        'cH': .503,
        'dA': lat_def,
        'dY': 0,
        'dX': hx_def,
        'nleft': nleft,
        'nright': nright,
        'ndef': ndef,
        'ref_index': 2.4028,
        'hole_type': 'ribAir',
        'effective_index': 1.6,
        'resonance_wavelength': 0.737,
        'defect_type': 'cubic',
        'min_dim': 5e-8
    }

    beam_length = 30
    nleft = pcc_params['nleft']
    ndef = pcc_params['ndef']
    nright = pcc_params['nright']
    a_list, hx_list, hy_list, flag = generate_cavity_and_check(pcc_params)
    a = pcc_params['a']
    h0 = pcc_params['cH'] * a
    w0 = pcc_params['cW'] * a

    if flag:
        return None, None, flag
    unit_cells = []

    engine = LumericalEngine(mesh_accuracy=5, hide=hide, lumerical_path=FDTDLoc, working_path="./fsps")

    for i in range(len(a_list)):
        # Take each unit cell and add it to our model. First update the unit cell parameters according to our cell list

        hole_params = {
            'a': a_list[i],
            'hx': hx_list[i],
            'hy': hy_list[i],
            'beam_width': w0,
            'beam_height': h0,
            'ref_index': 2.4028,
            'hole_type': pcc_params['hole_type']
        }

        hole, cell_size = build_hole(hole_params)
        unit_cells += [UnitCell(structures=hole, size=cell_size, engine=engine)]

    cavity = Cavity1D(
        unit_cells=unit_cells,
        structures=[BoxStructure(Vec3(0), Vec3(beam_length*1e-6, w0*1e-6, h0*1e-6),
                                 DielectricMaterial(2.4028, order=2, color="red"))],
        engine=engine,
    )

    cavity_name = np.array2string(cavity_params, prefix='', separator='_')
    cavity_name = cavity_name[1:]
    cavity_name = cavity_name[:-1]
    cavity_name = cavity_name.replace(" ", "")
    file_name = log_name[:-4] + "_" + cavity_name + "_" + str(iter_count)
    cavity.save(file_name + ".obj")
    #Once we check plot_geom, uncomment this
    #plot_geom(cavity, file_name + "_geom.png", hide)

    return cavity, file_name, flag


def fitness(cavity_params):
    global iter_count
    global rerun_thresh
    global target_frequency
    global source_frequency
    global lat_const
    cavity, file_name, flag = build_cavity(cavity_params)

    if flag:
        print('Fabrication intolerant')
        return 0
    else:
        man_mesh = MeshRegion(BBox(Vec3(0), Vec3(12e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)

        r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=60e-15,
                              analyze_time=600e-15, mesh_regions=[man_mesh], sim_size=Vec3(1.25, 3, 8))
        #r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions=[man_mesh], sim_size=Vec3(2, 3, 5.5))

        #r1 = cavity.simulate("resonance", target_freq=source_frequency)
        # r2 = cavity.get_results("resonance")[-1]
        # print(r2)
        # print(r2['res']["xyprofile"].max_loc())
        # print(r2['res']["yzprofile"].max_loc())
        # r2["xyprofile"].show()

        # Decide how we want to define qx, it's possible that we only want to define it as 'qxmax'
        # since that's the rightward flux. This would require a new definition of qscat

        qx = 1 / (1 / r1["qxmin"] + 1 / r1["qxmax"])
        qx1 = r1["qxmin"]
        qx2 = r1["qxmax"]
        qy = 1 / (2 / r1["qymax"])
        qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
        qscat = 1 / ((1 / qy) + (1 / qz))
        qtot = 1 / (1 / qx + 1 / qy + 1 / qz)
        vmode = r1["vmode"]
        vmode_copy = vmode
        vmode = 1e6 if vmode < 0.48 else vmode
        qtot_max = 300000
        purcell = qtot / vmode if qtot < qtot_max else qtot_max / vmode
        F = r1["freq"]
        wavelen = (2.99e8 / F) * 1e9

        print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
            r1["freq"], r1["vmode"],
            1 / (1 / r1["qxmin"] + 1 / r1["qxmax"]),
            1 / (2 / r1["qymax"] + 1 / r1["qzmin"] + 1 / r1["qzmax"])
        ))
        print(purcell)

        wavelen_pen = np.exp(-((target_frequency - F) / 4e12) ** 2)
        #qx_pen = (np.exp(-((target_qx - qx) / 120000) ** 2) + np.exp(-((target_qx - qx) / 60000) ** 2)
        #          + np.exp(-((target_qx - qx) / 30000) ** 2) + np.exp(-((target_qx - qx) / 2000) ** 2) / 4)
        qscat_max = 500000
        guidedness = qscat / qx if qscat < qscat_max else qscat_max / qx
        #witness = -1 * purcell * guidedness # * qx_pen
        witness = -1 * abs(qscat) / vmode if qscat < qscat_max else -1 * abs(qscat_max) / vmode
        if witness < -.07:
            lat_const = F / target_frequency * lat_const

        with open(log_name, "ab") as f:
            f.write(b"\n")
            step_info = np.append(cavity_params,np.array([witness, wavelen_pen, purcell, r1["qxmin"], r1["qxmax"], qscat, qtot, vmode, vmode_copy, F, lat_const]))
            np.savetxt(f, step_info.reshape(1, step_info.shape[0]), fmt='%.6f')
        # r1["xyprofile"].save(file_name+"_xy.png",title=f"Q = {qtot:.0f} \nQ_scat = {qscat:.04} Qx = {qx:.0f}\nV = {vmode_copy:.3f}")
        # r1["yzprofile"].save(file_name+"_yz.png",title=f"Q = {qtot:.0f} Q_scat = {qscat:.04}\n Qx1 = {qx1:.0f} Qx2 = {qx2:.0f}\nV = {vmode_copy:.3f} "+r"$\lambda$"+f" = {wavelen:.1f}")

        # second condition ensures that we only rerun once
        if (wavelen_pen < rerun_thresh) and (source_frequency == target_frequency):
            # shift source frequency to cavity resonance and rerun simulation.
            # (this should help avoid non cavities with artificially low mode volumes)
            source_frequency = F
            witness_rerun = fitness(cavity_params)
            print("rerun. Fitness when source is recentered:",witness_rerun)
            source_frequency = target_frequency
            return witness_rerun

        return witness

log_name = f"C:/Users/LoncarLab/Documents/GitHub/wvgsolver/examples/ribbed_cavity_optimization" \
           f"/optimization_022723_rib_cavity/optimal_rib_cav_04.txt"

lat_const = .287
p0 = np.array([.173, .625, 1.76])
# dA, cX, cY
bounds = ((.05, .25),  (.5, .9), (1, 2))

with open(log_name, "ab") as f:
    f.write(b"dA    cX  cY  fitness wavelen_pen purcell qxmin   qxmax   qscat   qtot    vmode   vmode_copy  F lat_const")

# witness = fitness(p0)
popt = minimize(fitness, p0, bounds=bounds, method='Nelder-Mead')
print(popt)
