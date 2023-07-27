from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, \
    PolygonStructure
from wvgsolver import Cavity1D, UnitCell, Vec3
import numpy as np
from wvgsolver.engine import LumericalEngine
import os


def build_hole(hole_params):
    a = hole_params['a']
    h0 = hole_params['beam_height']
    w0 = hole_params['beam_width']
    hx = hole_params['hx']
    hy = hole_params['hy']
    hole_type = hole_params['hole_type']

    if hole_type == 'ellipse':
        # offset the hole to respect the way we define the relevant lattice constant
        cell_hole = CylinderStructure(Vec3(0), h0, hx / 2, DielectricMaterial(1, order=1, color="red"),
                                      radius2=hy / 2)
        return cell_hole
    elif hole_type == 'rib':
        shift = hole_params['shift']
        coef = hy / (a - hx) ** 2
        npoints = 20
        bot = (w0 - hy / 2) / 2 + shift * hy / 4
        # creates the unit cell
        rib_up_verts = []
        for s in np.linspace(-a / 2, a / 2, num=npoints):
            x = s
            y = w0 / 2
            rib_up_verts.append((x, y))

        for s in np.linspace(a / 2, hx / 2, num=npoints):
            x = s
            y = -coef * (s - a / 2) ** 2 + bot + hy / 2
            rib_up_verts.append((x, y))

        for s in np.linspace(hx / 2, hx - a / 2, num=npoints):
            x = s
            y = bot + coef * (s - hx + a / 2) ** 2
            rib_up_verts.append((x, y))

        for s in np.linspace(a / 2 - hx, -hx / 2, num=npoints):
            x = s
            y = bot + coef * (s + hx - a / 2) ** 2
            rib_up_verts.append((x, y))

        for s in np.linspace(-hx / 2, -a / 2, num=npoints):
            x = s
            y = -coef * (s + a / 2) ** 2 + bot + hy / 2
            rib_up_verts.append((x, y))
        # rib_up_verts.append((-a / 2, w0 / 2))
        # rib_up_verts.append((-a / 2, w0))
        # rib_up_verts.append((a / 2, w0))
        # rib_up_verts.append((a / 2, w0 / 2))
        rib_up = PolygonStructure(pos=Vec3(0), verts=rib_up_verts, height=h0,
                                  material=DielectricMaterial(1, order=1))

        rib_down_verts = []
        for (x, y) in rib_up_verts:
            rib_down_verts.append((x, -y))
        rib_down = PolygonStructure(pos=Vec3(0), verts=rib_down_verts, height=h0,
                                    material=DielectricMaterial(1, order=1))

        return [rib_up, rib_down]


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
    ntaper = pcc_params['ntaper']

    a = pcc_params['a']
    hx = pcc_params['hx']
    hy = pcc_params['hy']

    # Make the left mirror region
    left_list = np.ones(nleft)
    lat_list = left_list * a
    hx_list = left_list * hx
    hy_list = left_list * hy

    # Make the cavity region
    if pcc_params['cav defect type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['cav defect type'] == 'quadratic':  ### NOTE THAT QUADRATIC IS NOT RIGHT Now
        defect = quadratic_defect

    defect_index = range(-ndef, ndef)
    lat_def = pcc_params['dA']
    hx_def = pcc_params['dX']
    hy_def = pcc_params['dY']
    lat_list = np.append(lat_list, a * np.array([defect(i, ndef, lat_def) for i in defect_index]))
    hx_list = np.append(hx_list, hx * np.array([defect(i, ndef, hx_def) for i in defect_index]))
    hy_list = np.append(hy_list, hy * np.array([defect(i, ndef, hy_def) for i in defect_index]))

    # Make the right mirror
    right_list = np.ones(nright)
    lat_list = np.append(lat_list, a * right_list)
    hx_list = np.append(hx_list, hx * right_list)
    hy_list = np.append(hy_list, hy * right_list)

    if ntaper > 0:
        # Make the taper region
        defect = linear_taper

        tap_index = range(0, ntaper, 1)
        lat_def_tap = pcc_params['dA_tap']
        hx_def_tap = pcc_params['dX_tap']
        hy_def_tap = pcc_params['dY_tap']

        lat_list = np.append(lat_list, a * np.array([defect(i, ntaper, lat_def_tap) for i in tap_index]))
        hx_list = np.append(hx_list, hx * np.array([defect(i, ntaper, hx_def_tap) for i in tap_index]))
        hy_list = np.append(hy_list, hy * np.array([defect(i, ntaper, hy_def_tap) for i in tap_index]))

    w0 = pcc_params['w0']
    aMir = pcc_params['a']
    aDef = aMir * (1 - lat_def)
    hyMir = hy
    hyDef = hyMir * (1 - hy_def)
    hxMir = hx
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
        return None ,False

    h0 = cp['h0']
    w0 = cp['w0']
    beam_length = 30 * 1e-6
    log_name = "E:/User Data/michael/lum_results"
    engine = LumericalEngine(mesh_accuracy=5, hide=hide, lumerical_path=FDTDLoc, working_path=log_name + "/fsps")
    ref_index = 2.4028
    unit_cells = []
    for (a, hx, hy) in zip(lat_list, hx_list, hy_list):
        # Take each unit cell and add it to our model. First update the unit cell parameters according to our cell list
        cell_size = Vec3(a, w0, h0)
        cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(ref_index, order=2, color="blue"))
        hole_params = {
            'a': a,
            'hx': hx,
            'hy': hy,
            'beam_width': w0,
            'beam_height': h0,
            'hole_type': cp['hole_type'],
            'shift': cp['shift']
        }

        hole = build_hole(hole_params)
        if cp['hole_type'] == 'rib':
            hole.append(cell_box)
            unit_cells += [UnitCell(structures=hole, size=cell_size, engine=engine)]
        else:
            unit_cells += [UnitCell(structures=[hole, cell_box], size=cell_size, engine=engine)]

    cav_shift = 0
    cavity = Cavity1D(
        unit_cells=unit_cells,
        structures=[BoxStructure(Vec3(0), Vec3(beam_length, w0, h0),
                                 DielectricMaterial(2.4028, order=2, color="red"))],
        engine=engine,
        center_shift=cav_shift
    )
    return cavity, True
