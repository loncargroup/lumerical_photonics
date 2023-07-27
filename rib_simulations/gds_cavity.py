from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D
import phidl.geometry as pg
import phidl.path as pp
from phidl import Path, CrossSection, Device
import numpy as np


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

    wvg = pg.straight((dev_length, beam_width)).move(
        destination=(c.xmin - back_length, -beam_width / 2))
    cav = pg.boolean(wvg, c, 'not')
    cav.add_port('cav_right_port', midpoint=(c.xmax + front_length, 0), orientation=0)

    wvg_right = pg.straight((front_length, beam_width)).movey(destination=-beam_width / 2)
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