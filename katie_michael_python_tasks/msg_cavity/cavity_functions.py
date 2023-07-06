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
        defect = linear_taper

        tap_index = range(0, ntaper, 1)
        lat_def_tap = pcc_params['dA_tap']
        hx_def_tap = pcc_params['dX_tap']
        hy_def_tap = pcc_params['dY_tap']

        lat_list = np.append(lat_list, a * np.array([defect(i, ntaper, lat_def_tap) for i in tap_index]))
        hx_list = np.append(hx_list, a * hx_scale * np.array([defect(i, ntaper, hx_def_tap) for i in tap_index]))
        hy_list = np.append(hy_list, a * hy_scale * np.array([defect(i, ntaper, hy_def_tap) for i in tap_index]))

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
                         analyze_time=600e-15, mesh_regions=[man_mesh], sim_size=Vec3(1.25, 3, 10))
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

    qperp = 1 / (1 / qy + 1 / qz)
    qscat = 1 / (1 / qperp + 1 / qleft)
    qwvg = qright

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
        'qscat': qscat,
        'qwvg': qwvg,
        'qleft': qleft,
        'vmode': vmode,
        # This is to keep track of wavelengths
        'freq': freq,
        'wl': c / freq,
    }

    print('---')
    return cavity_results



def gds_cavity(cavity, cp):
    filename = 'SiGcenterCavity_design6.gds'
    parser = DielectricExtrusionFaceGDSParser(cavity, invert=True)
    parser.save(filename)
    c = pg.import_gds(filename, cellname=None, flatten=False)

    # Let's build up the waveguide that we're going to etch from
    back_length = 3
    front_length = 3
    taper_length = 30
    cav_length = c.xmax - c.xmin
    beam_width = cp['w0'] * 1e6
    dev_length = cav_length + back_length + front_length

    wvg = phidl.geometry.straight((dev_length, beam_width)).move(
        destination=(c.xmin - back_length, -beam_width / 2))
    cav = pg.boolean(wvg, c, 'not')
    cav.add_port('cav_right_port', midpoint=(c.xmax + front_length, 0), orientation=0)

    wvg_right = phidl.geometry.straight((front_length, beam_width)).movey(destination=-beam_width / 2)
    wvg_right.add_port('wvg_right_port', midpoint=(front_length, 0), orientation=0)
    wvg_right.add_port('wvg_left_port', midpoint=(0, 0), orientation=180)

    taper = pg.taper(taper_length, width1=beam_width, width2=.06)
    taper.add_port('taper_left_port', midpoint=(0, 0), orientation=180)

    D = Device()
    tether_dev = add_tether(beam_width)
    tether_ref = D.add_ref(tether_dev)
    wvg_ref = D.add_ref(wvg_right)
    taper_ref = D.add_ref(taper)

    cav_ref = D.add_ref(cav)
    tether_ref.connect(port=tether_ref.ports['left_port'], destination=cav_ref.ports['cav_right_port'])
    wvg_ref.connect(port=wvg_ref.ports['wvg_left_port'], destination=tether_ref.ports['right_port'])
    taper_ref.connect(port=taper_ref.ports['taper_left_port'], destination=wvg_ref.ports['wvg_right_port'])
    trench = pg.rectangle((D.xmax - D.xmin + 10, D.ymax - D.ymin)).move(destination=(D.xmin, D.ymin))
    cav_invert = pg.boolean(trench, D, 'not')
    return cav_invert


def add_tether(beam_width):
    trench_width = 4
    gauss_length = 6
    gauss_width = 1.5 * beam_width
    vteth_width = .1
    sigma = .05

    def gaussian(t):
        x = t - 1 / 2
        return gauss_width * np.exp(-x ** 2 / sigma) + beam_width * (1 - np.exp(-x ** 2 / sigma))

    P = pp.straight(gauss_length)
    X = CrossSection()
    X.add(width=gaussian, offset=0, layer=0)
    D = P.extrude(X)

    vert_teth = pp.straight(trench_width).rotate(90).movex(gauss_length / 2).movey(-trench_width / 2)
    X = CrossSection().add(width=vteth_width)
    D = pg.boolean(D, vert_teth.extrude(X), 'or')
    D = pg.union(D)

    D.add_port('left_port', midpoint=(0, 0), orientation=180)
    D.add_port('right_port', midpoint=(gauss_length, 0), orientation=0)
    return D
