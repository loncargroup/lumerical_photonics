import numpy as np
import gdspy
import math
from scipy import special
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import time

W = .54
theta = .2 * np.pi
Ndef = 6
Nright = 20
Nleft = 20

# Tether Parameters
W_g = 6 * W
L_L = 1
L_R = 4
W_T = 1.3 * W
L_O = 5 * W
L_T = .2
W_Tri = 1
L_Tri = 2


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
    global R
    return (d + 2 * ((Wm / 2 - R * np.sin(theta)) * np.tan(theta) + R * (1 - np.cos(theta))))


## Construction Functions

def bridge(x0, y0, d):
    x1 = d / 2 + R * (1 - np.cos(theta))
    y1 = R * np.sin(theta)
    x2 = L_an(d) / 2
    y2 = Wm / 2
    y3 = W / 2
    # Geometry must be placed in cells.
    c3 = gdspy.Curve(x1 + x0, y1 + y0).L(x2 + x0, y2 + y0, x2 + x0, y3 + y0, -x2 + x0, y3 + y0, -x2 + x0, y2 + y0,
                                         -x1 + x0, y1 + y0)  # Upper Polygon
    c4 = gdspy.Curve(x1 + x0, -y1 + y0).L(x2 + x0, -y2 + y0, x2 + x0, -y3 + y0, -x2 + x0, -y3 + y0, -x2 + x0, -y2 + y0,
                                          -x1 + x0, -y1 + y0)  # Lower Polygon
    p4 = gdspy.Polygon(c4.get_points())
    p3 = gdspy.Polygon(c3.get_points())
    x = gdspy.boolean(gdspy.Rectangle((-x1 + x0, -y1 + y0), (x1 + x0, y1 + y0)),
                      gdspy.boolean(
                          gdspy.Round((d / 2 + R + x0, y0), R, initial_angle=np.pi - theta, final_angle=np.pi + theta,
                                      tolerance=.0001),
                          gdspy.Round((-d / 2 - R + x0, y0), R, initial_angle=-theta, final_angle=theta,
                                      tolerance=.0001), "or"), "not")  # Bridge Construction
    return gdspy.boolean(gdspy.boolean(p4, p3, "or"), x, "or")


def cavity(d0, dMax, Nright):
    a = bridge(0, 0, d0)
    x0 = 0
    for i in range(1, Ndef + 1):
        x0 += (L_an(bD(i - 1, d0, dMax)) + L_an(bD(i, d0, dMax))) / 2
        a = gdspy.boolean(a, bridge(x0, 0, bD(i, d0, dMax)), "or")
        a = gdspy.boolean(a, bridge(-x0, 0, bD(i, d0, dMax)), "or")
    xR = x0
    xL = x0
    for i in range(Nright):
        xR += L_an(bD(Ndef, d0, dMax))
        a = gdspy.boolean(a, bridge(xR, 0, bD(Ndef, d0, dMax)), "or")
    for i in range(Nleft):
        xL += L_an(bD(Ndef, d0, dMax))
        a = gdspy.boolean(a, bridge(-xL, 0, bD(Ndef, d0, dMax)), "or")
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

    return gdspy.boolean(p1, p2, 'or')


## Assembler Functions
def device(xc, yc, d0, dMax, Nright):
    a = gdspy.boolean(cavity(d0, dMax, Nright),
                      gdspy.Rectangle((-lenLeft(d0, dMax) - L_L, -W / 2), (-lenLeft(d0, dMax), W / 2)), "or")
    b = taper2(d0, dMax, Nright)
    c = gdspy.boolean(a, b, "or")
    c.translate(xc, yc)
    return c


def text2(xc, yc, z, rot=np.pi / 4):
    text = gdspy.Text(str(z), 4, (0, 0))
    text.rotate(rot)
    text.translate(xc, yc)
    return text


def deviceInverse2(xc, yc, d0, dMax, Nright, Wm0,  z, rot=np.pi / 4, R0 = .1):
    global R
    global Wm
    R = R0
    Wm = Wm0
    tLeft = lenLeft(d0, dMax) + L_L
    tRight = lenRight(d0, dMax, Nright) + L_R + 2 * L_O + 30 + 30

    a = device(0, 0, d0, dMax, Nright)
    rect = gdspy.Rectangle((-tLeft, -W_g), (tRight, W_g))
    c = gdspy.boolean(rect, a, "not")
    c.rotate(rot)
    c.translate(xc, yc)
    if z < 100:
        return gdspy.boolean(c, text2(xc - lenLeft(d0, dMax) - 6, yc - 16, z, rot), "or", layer=1)
    else:
        return gdspy.boolean(c, text2(xc - lenLeft(d0, dMax) - 8, yc - 19, z, rot), "or", layer=1)