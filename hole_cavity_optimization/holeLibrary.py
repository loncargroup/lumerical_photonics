# holeLibrary: this is a set of code that takes in hole parameters and spits out a unit cell with those properties
# Inputs: a, hx, hy, w0, h0 in units of um
# Outputs: wvgsolver geometric structure

from wvgsolver import Vec3
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, \
    PolygonStructure
import numpy as np
import gdspy


def build_hole(hole_params):
    a = hole_params['a']*1e-6
    h0 = hole_params['beam_height']*1e-6
    w0 = hole_params['beam_width']*1e-6
    hx = hole_params['hx']*1e-6
    hy = hole_params['hy']*1e-6
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
    elif hole_type == 'ribAir':
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
        return [rib_up,rib_down,cell_box], cell_size
    elif hole_type == 'tri_fil':
        r1 = 0

        c = gdspy.Curve(-hx/2, 0).L(0, hy/2, hx/2, 0, 0, -hy/2, hx/2, 0)
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

