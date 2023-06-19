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
from phidl import quickplot as qp
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
        rib_up_verts.append((-a / 2, w0 / 2))
        rib_up_verts.append((-a / 2, w0))
        rib_up_verts.append((a / 2, w0))
        rib_up_verts.append((a / 2, w0 / 2))
        rib_up = PolygonStructure(pos=Vec3(0), verts=rib_up_verts, height=h0,
                                  material=DielectricMaterial(1, order=1))

        rib_down_verts = []
        for (x, y) in rib_up_verts:
            rib_down_verts.append((x, -y))
        rib_down = PolygonStructure(pos=Vec3(0), verts=rib_down_verts, height=h0,
                                    material=DielectricMaterial(1, order=1))

        return [rib_up, rib_down]
    # elif hole_type == 'tri_fil':
    #     r1 = 0
    #
    #     c = gdspy.Curve(-hx / 2, 0).L(0, hy / 2, hx / 2, 0, 0, -hy / 2, hx / 2, 0)
    #     print(c.get_points())
    #     poly = gdspy.Polygon(c.get_points())
    #     if r1 > 0:
    #         poly.fillet(r1)
    #
    #     verts = []
    #     for i in poly.polygons[0]:
    #         verts.append((i[0], i[1]))
    #     print(verts)
    #     tri = PolygonStructure(pos=Vec3(0), verts=verts, height=h0,
    #                            material=DielectricMaterial(1, order=1))
    #     return tri


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

        lat_list = np.append(lat_list, a * np.array([defect(i, ntaper, lat_def_tap) for i in tap_index]))
        hx_list = np.append(hx_list, a * hx_scale * np.array([defect(i, ntaper, hx_def_tap) for i in tap_index]))
        hy_list = np.append(hy_list, a * hy_scale * np.array([defect(i, ntaper, hy_def_tap) for i in tap_index]))

    return [lat_list, hx_list, hy_list]


def make_cavity(cp):
    hide = False
    FDTDLoc = 'C:/Program Files/Lumerical/v221/'  # FDTDLoc, I will show you where this is later but this is how you can do it local
    FDTDexeLoc = os.path.join(FDTDLoc, 'bin/fdtd-solutions')
    FDTDmpiLoc = os.path.join(FDTDLoc, 'bin/fdtd-engine-ompi-lcl')
    [lat_list, hx_list, hy_list] = generate_cavity_and_check(cp)
    a = cp['a']
    h0 = cp['cH'] * a
    w0 = cp['cW'] * a
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
    return cavity


def simulate_cavity(cavity, cp):
    source_frequency = cp['source frequency']
    man_mesh = MeshRegion(BBox(Vec3(0), Vec3(12e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)

    # There are a few different simulation conditions, let's keep this one but see if we need to change later
    r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=60e-15, analyze_fspan=10e12,
                         analyze_time=600e-15, mesh_regions=[man_mesh], sim_size=Vec3(1.25, 3, 8))
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions=[man_mesh],
    #                     sim_size=Vec3(2, 3, 8))

    # r1 = cavity.simulate("resonance", target_freq=source_frequency)

    qleft = r1["qxmin"]
    qright = r1["qxmax"]
    qy = 1 / (2 / r1["qymax"])
    qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
    return qz


def gds_cavity(cavity, cp):
    filename = 'SiGcenterCavity_design6.gds'
    parser = DielectricExtrusionFaceGDSParser(cavity, invert=True)
    parser.save(filename)
    c = pg.import_gds(filename, cellname=None, flatten=False)
    qp(c)

    # Let's build up the waveguide that we're going to etch from
    back_length = 3
    front_length = 3
    cav_length = c.xmax - c.xmin
    beam_width = cp['a'] * cp['cW'] * 1e6

    dev_length = cav_length+back_length+front_length
    wvg = phidl.geometry.straight((dev_length, beam_width)).move(
        destination=(c.xmin - back_length, -beam_width / 2))
    add_tether(beam_width)
    cav = pg.boolean(wvg, c, 'not')
    cav.write_gds(filename)

def add_tether(beam_width):
    trench_width = 4
    gauss_length = 6
    gauss_width = 1.5*beam_width
    sigma = .1

    def gaussian(t):
        x = t-1/2
        return gauss_width*np.exp(-x**2/sigma)+beam_width*(1-np.exp(-x**2/sigma))

    P = pp.straight(gauss_length)
    X = CrossSection()
    X.add(width=gaussian, offset=0, layer=0)
    D = P.extrude(X)
    qp(D)
