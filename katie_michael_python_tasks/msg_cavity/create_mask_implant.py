# Goal of this code is to create standardized mask implantation that we can use for
# future omcs and rib designs

import phidl.geometry as pg
from phidl import Device, CrossSection
import phidl.path as pp
import phidl
import numpy as np
import csv
import sys
import os

path = os.path.abspath("/Users/michaelhaas/iCloud Drive (Archive)/Documents/GitHub/phidl-work")
sys.path.append(path)
import cavity_functions as cf


def make_aperture(x, y, width, height, mask_angle):
    return pg.rectangle((width, height)).move(destination=(x - width / 2, y - height / 2)).rotate(mask_angle, (x, y))


def make_marker(x0, y0):
    smallSq = 1.75
    largeSq = 1.975

    mark = Device()
    big_rect: Device = pg.rectangle((5, 5)).move(destination=(x0 - 2.5, y0 - 2.5))
    mark.add(pg.rectangle((largeSq, largeSq)).move(destination=(x0 - 2, y0 - 2)))
    mark.add(pg.rectangle((smallSq, smallSq)).move(destination=(x0 - 2, y0 + 2 - smallSq)))
    mark.add(pg.rectangle((largeSq, largeSq)).move(destination=(x0 + 2 - largeSq, y0 + 2 - largeSq)))
    mark.add(pg.rectangle((smallSq, smallSq)).move(destination=(x0 + 2 - smallSq, y0 - 2)))
    return mark, big_rect


def make_tri(x0, y0):
    tri = pg.taper(20, width1=20, width2=0).rotate(90).move(destination=(x0, y0))
    bigRect = pg.rectangle((25, 25)).move(destination=(x0 - 12.5, y0 - 2.5))
    return tri, bigRect


def cavity(lat_const, lat_def, hx):
    global z
    z += 1
    hx_def = lat_def
    pcc_params = {
        'a': lat_const,
        'cX': hx,
        'cY': 1.4,
        'cW': 2,
        'cH': .5,
        'dA': 0,
        'dY': 0,
        'dX': 0,

        'nleft': 20,
        'nright': 7,
        'ndef': 6,
        'hole_type': 'ellipse',
        'defect_type': 'cubic',
        'min_dim': 1e-7
    }
    cav = cf.make_cavity(pcc_params)
    # text = pg.text(str(z),size= 4).move(destination=(0,5))
    # return pg.boolean(cav,text,'or')
    return cav


# Let's make an omc
# Cavity Constants
W = .54
R = .175

Wm = .24  # Will vary this later, keep for testing
Nright = 20  # Will vary this later, keep for testing

theta = .2 * np.pi
Ndef = 6
Nleft = 20

# Tether Constants
W_g = 2
L_L = 1
L_R = 4
W_T = 1.3 * W
L_O = 5 * W
L_T = .6 * W
W_Tri = 1
L_Tri = 2


## Length Functions
def bD(x, d0, dMax):
    dS = (dMax - d0) / Ndef ** 2
    return (dS * x ** 2 + d0)


def lenLeft(d0, dMax):
    sum = L_an(d0) / 2
    for i in range(1, Ndef):
        sum += L_an(bD(i, d0, dMax))
    sum += (Nleft) * L_an(bD(Ndef, d0, dMax))
    return sum


def lenRight(d0, dMax, Nright):
    sum = L_an(d0) / 2
    for i in range(1, Ndef):
        sum += L_an(bD(i, d0, dMax))
    sum += (Nright) * L_an(bD(Ndef, d0, dMax))
    return sum


def L_an(d):
    return (d + 2 * ((Wm / 2 - R * np.sin(theta)) * np.tan(theta) + R * (1 - np.cos(theta))))


## Construction Functions


def bridge(x0, y0, d):
    x1 = d / 2 + R * (1 - np.cos(theta))
    y1 = R * np.sin(theta)
    x2 = L_an(d) / 2
    y2 = Wm / 2
    y3 = W / 2

    lat_const = L_an(d)
    top_block = pg.rectangle((lat_const, y3 - y2)).move(destination=(-lat_const / 2, y2))
    top_taper = pg.taper(y2 - y1, width1=lat_const, width2=2 * x1)
    center_block = pg.rectangle()

    c3 = gdspy.Curve(x1 + x0, y1 + y0).L(x2 + x0, y2 + y0, x2 + x0, y3 + y0, -x2 + x0, y3 + y0, -x2 + x0, y2 + y0,
                                         -x1 + x0, y1 + y0)  # Upper Polygon
    c4 = gdspy.Curve(x1 + x0, -y1 + y0).L(x2 + x0, -y2 + y0, x2 + x0, -y3 + y0, -x2 + x0, -y3 + y0, -x2 + x0, -y2 + y0,
                                          -x1 + x0, -y1 + y0)  # Lower Polygon
    p4 = gdspy.Polygon(c4.get_points())
    p3 = gdspy.Polygon(c3.get_points())
    x = pg.boolean(gdspy.Rectangle((-x1 + x0, -y1 + y0), (x1 + x0, y1 + y0)),
                   pg.boolean(
                       gdspy.Round((d / 2 + R + x0, y0), R, initial_angle=np.pi - theta, final_angle=np.pi + theta,
                                   tolerance=.0001),
                       gdspy.Round((-d / 2 - R + x0, y0), R, initial_angle=-theta, final_angle=theta, tolerance=.0001),
                       "or"), "not")  # Bridge Construction
    return pg.boolean(pg.boolean(p4, p3, "or"), x, "or")


def cavity(d0, dMax, Nright):
    a = bridge(0, 0, d0)
    x0 = 0
    for i in range(1, Ndef + 1):
        x0 += (L_an(bD(i - 1, d0, dMax)) + L_an(bD(i, d0, dMax))) / 2
        a = pg.boolean(a, bridge(x0, 0, bD(i, d0, dMax)), "or")
        a = pg.boolean(a, bridge(-x0, 0, bD(i, d0, dMax)), "or")
    xR = x0
    xL = x0
    for i in range(Nright):
        xR += L_an(bD(Ndef, d0, dMax))
        a = pg.boolean(a, bridge(xR, 0, bD(Ndef, d0, dMax)), "or")
    for i in range(Nleft):
        xL += L_an(bD(Ndef, d0, dMax))
        a = pg.boolean(a, bridge(-xL, 0, bD(Ndef, d0, dMax)), "or")
    return a


def taper2(d0, dMax, Nright):
    xc = lenRight(d0, dMax, Nright)
    yc = W / 2
    x1 = xc + L_R
    y1 = yc
    x2 = x1 + (L_O - L_T) / 2
    y2 = y1 + (W_T - W) / 2
    x3 = x2
    y3 = y2 + W_g - (W_T + W_Tri) / 2
    x4 = x3 - (L_Tri - L_T) / 2
    y4 = y3 + W_Tri
    x5 = x4 + L_Tri
    y5 = y4
    x6 = x5 - (L_Tri - L_T) / 2
    y6 = y5 - W_Tri
    x7 = x6
    y7 = y6 - W_g + (W_T + W_Tri) / 2
    x8 = x7 + (L_O - L_T) / 2
    y8 = y7 - (W_T - W) / 2
    x9 = x8 + 30
    y9 = .03
    x10 = x9
    y10 = -y9
    x11 = x8
    y11 = -y8
    x12 = x7
    y12 = -y7
    x13 = x6
    y13 = -y6
    x14 = x5
    y14 = -y5
    x15 = x4
    y15 = -y4
    x16 = x3
    y16 = -y3
    x17 = x2
    y17 = -y2
    x18 = x1
    y18 = -y1
    x19 = xc
    y19 = -yc
    c = gdspy.Curve(xc, yc).L(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9
                              , x10, y10, x11, y11, x12, y12, x13, y13, x14, y14, x15, y15, x16, y16, x17, y17, x18,
                              y18, x19, y19)

    c2 = gdspy.Curve(xc + 2 * L_O, yc).L(x1 + 2 * L_O, y1, x2 + 2 * L_O, y2, x3 + 2 * L_O, y3, x4 + 2 * L_O, y4,
                                         x5 + 2 * L_O, y5, x6 + 2 * L_O, y6, x7 + 2 * L_O, y7, x8 + 2 * L_O, y8,
                                         x9 + 2 * L_O, y9
                                         , x10 + 2 * L_O, y10, x11 + 2 * L_O, y11, x12 + 2 * L_O, y12, x13 + 2 * L_O,
                                         y13, x14 + 2 * L_O, y14, x15 + 2 * L_O, y15, x16 + 2 * L_O, y16, x17 + 2 * L_O,
                                         y17, x18 + 2 * L_O, y18, x19 + 2 * L_O, y19)
    p1 = gdspy.Polygon(c.get_points())
    p2 = gdspy.Polygon(c2.get_points())

    return pg.boolean(p1, p2, 'or')


## Assembler Functions
def device(xc, yc, d0, dMax, Nright):
    a = pg.boolean(cavity(d0, dMax, Nright),
                   gdspy.Rectangle((-lenLeft(d0, dMax) - L_L, -W / 2), (-lenLeft(d0, dMax), W / 2)), "or")
    b = taper2(d0, dMax, Nright)
    c = pg.boolean(a, b, "or")
    c.translate(xc, yc)
    return c


def text2(xc, yc, z, rot=np.pi / 4):
    w = 1
    h = 5
    if z > 9:
        xc = xc - 5
    text = pg.text(str(z), 6).move(destination=(xc - 3, yc - 3))
    text.rotate(rot, center=(xc - 3, yc - 3))
    return text


def deviceInverse2(xc, yc, d0, dMax, Nright, z, rot=np.pi / 4):
    tLeft = lenLeft(d0, dMax) + L_L
    tRight = lenRight(d0, dMax, Nright) + L_R + 2 * L_O + 30 + 30

    a = device(0, 0, d0, dMax, Nright)
    trench_length = tRight - tLeft
    trench_height = 2 * W_g
    rect = pg.rectangle((trench_length, trench_height)).move(destination=(tLeft, -W_g))
    c = pg.boolean(rect, a, "not")
    c.rotate(rot)
    c.translate(xc, yc)
    if int(z) < 10:
        c = pg.boolean(c, text2(xc - lenLeft(d0, dMax) - 10, yc - 4, z, rot), "or", layer=1)
    else:
        c = pg.boolean(c, text2(xc - lenLeft(d0, dMax) - 15, yc - 4, z, rot), "or", layer=1)
    return c


lat_const_sweep = [.275]
lat_def_sweep = np.linspace(.16, .19, 3)
hx_sweep = np.linspace(.65, .7, 3)
E = Device()

dose = 1.5e12
mask_angle = 45
mask_angle_rad = mask_angle * np.pi / 180
depth = 100
filename = '052823_maskDev_' + str(int(dose / 1e11)) + 'e11_' + str(mask_angle) + 'deg_' + str(depth)

array_spacer_y = 15
array_spacer_x = 100
subarray_spacer_y = 3
marker_spacer = 420

subarray_row_num = 5
mask_size = [.05, .065, .08, .095, .11]

if mask_angle == 45:
    array_col_num = 5
    array_row_num = 12
    max_row_vals = 36

    marker_x = -130
    marker_y = 10
else:
    array_col_num = 4
    array_row_num = 35
    marker_x = -20
    marker_y = -20

parameterList = [['Chip Name = ' + filename, 'Dose' + str(dose)],
                 ['Column', 'Row', 'Subarray row', 'X Position', 'Y Postion', 'Width', 'Height']]
z = 0
parameterList = np.array(parameterList)
for i in range(array_col_num):
    x_pos = i * array_spacer_x
    if mask_angle == 45:  # Set this differently so that we can maximise device packing
        mid_val = (array_col_num - 1) / 2
        rank = abs(i - mid_val)

        array_row_vals = int(-(max_row_vals - array_row_num) / mid_val * rank + max_row_vals)
        y_start = (max_row_vals - array_row_num) / mid_val * (rank - mid_val) * array_spacer_y / 2
    else:
        array_row_vals = array_row_num
        y_start = 0

    for j in range(array_row_vals):  # Make array rows
        y_pos = j * array_spacer_y + y_start
        z += 1
        for k in range(subarray_row_num):  # Make subarray rows
            # Make the appropriate rotations
            x = x_pos * np.cos(mask_angle_rad) - (y_pos + k * subarray_spacer_y) * np.sin(mask_angle_rad)
            y = x_pos * np.sin(mask_angle_rad) + (y_pos + k * subarray_spacer_y) * np.cos(mask_angle_rad)
            E.add(make_aperture(x, y, mask_size[k], mask_size[k], mask_angle))  # Actually make the aperture
            if k == 1:
                E.add(cavity(.25, .16, .6).move(destination=(x, y)).rotate(mask_angle, center=(x, y)))
            parameterList = np.append(parameterList, [i, j, k, int(x * 1e3), int(y * 1e3), int(mask_size[k] * 1e3),
                                                      int(mask_size[k] * 1e3)], axis=0)

print(z)
for i in [0, 1]:
    for j in [0, 1]:
        mark, bigRect = make_marker(marker_x + i * marker_spacer, marker_y + j * marker_spacer)
        inverse = pg.boolean(bigRect, mark, 'not', layer=1)
        E.add(inverse)

for i in [.35, .65]:
    for j in [0, 1]:
        mark, bigRect = make_marker(marker_x + i * marker_spacer, marker_y + j * (marker_spacer + 55) - 35)
        inverse = pg.boolean(bigRect, mark, 'not', layer=1)
        E.add(inverse)
triangle, bigRect = make_tri(marker_x, marker_y - 45)
inverse = pg.boolean(bigRect, triangle, 'not', layer=1)
E.add(inverse)
triangle, bigRect = make_tri(marker_x + marker_spacer, marker_y - 45)
inverse = pg.boolean(bigRect, triangle, 'not', layer=1)
E.add(inverse)

bigAper_x = marker_x + 60
big_Aper_y = marker_y - 35

# label = pg.text(str(1),size = 30).move(destination=(marker_x+marker_spacer-30,marker_y-25))
# E.add(label)

E.add(make_aperture(bigAper_x, big_Aper_y, 10, 10, 0))

E.add(make_aperture(bigAper_x - 12, big_Aper_y, .06, .06, 0))
E.add(make_aperture(bigAper_x + 12, big_Aper_y, .06, .06, 0))
E.add(make_aperture(bigAper_x, big_Aper_y + 12, .06, .06, 0))

E.add(make_aperture(bigAper_x - 10, big_Aper_y, .1, .1, 0))
E.add(make_aperture(bigAper_x + 10, big_Aper_y, .1, .1, 0))
E.add(make_aperture(bigAper_x, big_Aper_y + 10, .1, .1, 0))

E.add(make_aperture(bigAper_x - 8, big_Aper_y, .2, .2, 0))
E.add(make_aperture(bigAper_x + 8, big_Aper_y, .2, .2, 0))
E.add(make_aperture(bigAper_x, big_Aper_y + 8, .2, .2, 0))

E.move(origin=(-marker_x - 210, -marker_y - 200))
parameterList = np.asarray(parameterList)
# with open(filename+'_params.csv', 'w',newline='') as file:
#     mywriter = csv.writer(file, delimiter=',')
#     mywriter.writerows(parameterList)  

E.write_gds(filename + '.gds')
