from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, \
    PolygonStructure
import gdspy
import numpy as np
import matplotlib.pyplot as plt
import os

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine


def build_hole(hole_params):
    a = hole_params['a'] * 1e-6
    h0 = hole_params['beam_height'] * 1e-6
    w0 = hole_params['beam_width'] * 1e-6
    hx = hole_params['hx'] * 1e-6
    hy = hole_params['hy'] * 1e-6
    ref_index = hole_params['ref_index']
    hole_type = hole_params['hole_type']
    cell_size = Vec3(a, w0, h0)
    if hole_type == 'ellipse':
        cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(ref_index, order=2, color="blue"))

        # offset the hole to respect the way we define the relevant lattice constant
        cell_hole = CylinderStructure(Vec3(0), h0, hx / 2, DielectricMaterial(1, order=1, color="red"),
                                      radius2=hy / 2)
        return [cell_box, cell_hole], cell_size
    elif hole_type == 'ribDiel':
        rib_verts = []
        npoints = 40
        coef = hy / (a - hx) ** 2
        bot = (w0 - hy) / 2
        # Describes the upper part of the rib, goes from top left edge of rib.

        # creates the unit cell
        for s in np.linspace(-a / 2, -hx / 2, num=npoints):
            x = s
            y = -coef * (s + a / 2) ** 2 + w0 / 2
            rib_verts.append((x, y))

        for s in np.linspace(-hx / 2, a / 2 - hx, num=npoints):
            x = s
            y = bot + coef * (s + hx - a / 2) ** 2
            rib_verts.append((x, y))

        for s in np.linspace(hx - a / 2, hx / 2, num=npoints):
            x = s
            y = bot + coef * (s - hx + a / 2) ** 2
            rib_verts.append((x, y))

        for s in np.linspace(hx / 2, a / 2, num=npoints):
            x = s
            y = -coef * (s - a / 2) ** 2 + w0 / 2
            rib_verts.append((x, y))

        # This continues the polygon in a counterclockwise fashion, starting at bottom right edge of rib and moves left.

        # creates the inverted unit
        for s in np.linspace(a / 2, hx / 2, num=npoints):
            x = s
            y = (-coef * (s - a / 2) ** 2 + w0 / 2) * -1
            rib_verts.append((x, y))

        for s in np.linspace(hx / 2, hx - a / 2, num=npoints):
            x = s
            y = (bot + coef * (s - hx + a / 2) ** 2) * -1
            rib_verts.append((x, y))

        for s in np.linspace(a / 2 - hx, -hx / 2, num=npoints):
            x = s
            y = (bot + coef * (s + hx - a / 2) ** 2) * -1
            rib_verts.append((x, y))

        for s in np.linspace(-hx / 2, -a / 2, num=npoints):
            x = s
            y = (-coef * (s + a / 2) ** 2 + w0 / 2) * -1
            rib_verts.append((x, y))

        rib = PolygonStructure(pos=Vec3(0), verts=rib_verts, height=h0, material=DielectricMaterial(ref_index, order=1))
        return [rib], cell_size
    elif hole_type == 'rib':
        npoints = 40
        coef = hy / (a - hx) ** 2
        bot = (w0 - hy) / 2
        # creates the unit cell
        rib_up_verts = []
        for s in np.linspace(-a / 2, a / 2, num=npoints):
            x = s
            y = w0 / 2
            rib_up_verts.append((x, y))

        for s in np.linspace(a / 2, hx / 2, num=npoints):
            x = s
            y = -coef * (s - a / 2) ** 2 + w0 / 2
            rib_up_verts.append((x, y))

        for s in np.linspace(hx / 2, hx - a / 2, num=npoints):
            x = s
            y = bot + coef * (s - hx + a / 2) ** 2
            rib_up_verts.append((x, y))

        for s in np.linspace(a / 2 - hx, -hx / 2, num=npoints):
            x = s
            y = bot + coef * (s + hx - a / 2) ** 2
            rib_up_verts.append((x, y))

        for s in np.linspace(-a / 2, -hx / 2, num=npoints):
            x = s
            y = -coef * (s + a / 2) ** 2 + w0 / 2
            rib_up_verts.append((x, y))
        rib_up = PolygonStructure(pos=Vec3(0), verts=rib_up_verts, height=h0,
                                  material=DielectricMaterial(1, order=1))

        rib_down_verts = []
        for s in np.linspace(-a / 2, a / 2, num=npoints):
            x = s
            y = w0 / 2
            rib_down_verts.append((x, -y))

        for s in np.linspace(a / 2, hx / 2, num=npoints):
            x = s
            y = -coef * (s - a / 2) ** 2 + w0 / 2
            rib_down_verts.append((x, -y))

        for s in np.linspace(hx / 2, hx - a / 2, num=npoints):
            x = s
            y = bot + coef * (s - hx + a / 2) ** 2
            rib_down_verts.append((x, -y))

        for s in np.linspace(a / 2 - hx, -hx / 2, num=npoints):
            x = s
            y = bot + coef * (s + hx - a / 2) ** 2
            rib_down_verts.append((x, -y))

        for s in np.linspace(-a / 2, -hx / 2, num=npoints):
            x = s
            y = -coef * (s + a / 2) ** 2 + w0 / 2
            rib_down_verts.append((x, -y))
        rib_down = PolygonStructure(pos=Vec3(0), verts=rib_down_verts, height=h0,
                                    material=DielectricMaterial(1, order=1))

        cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(ref_index, order=2, color="blue"))
        return [rib_up, rib_down, cell_box], cell_size
    elif hole_type == 'tri_fil':
        r1 = 0

        c = gdspy.Curve(-hx / 2, 0).L(0, hy / 2, hx / 2, 0, 0, -hy / 2, hx / 2, 0)
        print(c.get_points())
        poly = gdspy.Polygon(c.get_points())
        if r1 > 0:
            poly.fillet(r1)

        verts = []
        for i in poly.polygons[0]:
            verts.append((i[0], i[1]))
        print(verts)
        tri = PolygonStructure(pos=Vec3(0), verts=verts, height=h0,
                               material=DielectricMaterial(1, order=1))
        cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(ref_index, order=2, color="blue"))
        return [tri, cell_box], cell_size


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

    # Make the right mirror
    right_list = np.ones(nright)
    lat_list = np.append(lat_list, a * right_list)
    hx_list = np.append(hx_list, a * hx_scale * right_list)
    hy_list = np.append(hy_list, a * hy_scale * right_list)

    if ntaper > 0:
        # Make the taper region
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


def build_cavity(cp, log_name, FDTDLoc):
    hide = False
    # Put this for now, but this might not be the way we make it
    pcc_params = {
        'ref_index': 2.4028,
        'effective_index': 1.6,
        'resonance_wavelength': 0.737,
        'cav_defect_type': 'cubic',
        'tap_defect_type': 'cubic',
        'min_dim': 5e-8
    }

    pcc_params.update(cp)

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

    engine = LumericalEngine(mesh_accuracy=5, hide=hide, lumerical_path=FDTDLoc, working_path=log_name + "/fsps")

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

    cavity_name = 'test'
    cavity_name = cavity_name[1:]
    cavity_name = cavity_name[:-1]
    cavity_name = cavity_name.replace(" ", "")
    file_name = log_name[:-4] + "_" + cavity_name
    file_name = log_name[:-4] + "_rib_00_MidSim.obj"
    cavity.save(file_name)
    # Once we check plot_geom, uncomment this
    # plot_geom(cavity, file_name + "_geom.png", hide)

    return cavity, file_name, engine
