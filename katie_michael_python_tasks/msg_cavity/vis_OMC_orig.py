import os
import numpy as np
import gdspy
import math
from scipy import special

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import time

start = time.time()
print("hello")

# The GDSII file is called a library, which contains multiple cells.
lib = gdspy.GdsLibrary()
cell = lib.new_cell('FIRST')
#Device Parameters
W = .54
Wm = .24
theta = .2*np.pi
Ndef = 6
Nright = 20
Nleft = 20

#Tether Parameters
W_g = 6*W 
L_L = 1
L_R = 4 
W_T = 1.3*W
L_O = 5*W
L_T = .6*W
W_Tri = 1
L_Tri = 2

#Length Functions

def bD(x, d0, dMax):
    dS = (dMax-d0)/Ndef**2
    return (dS*x**2 + d0)
def lenLeft(d0,dMax,R,Wm,Nleft=Nleft):
    sum = L_an(d0,R,Wm)/2
    for i in range(1,Ndef):
        sum+= L_an(bD(i,d0,dMax),R,Wm)
    sum+= (Nleft)*L_an(bD(Ndef, d0,dMax),R,Wm)
    print(sum)
    return sum
def lenRight(d0,dMax,R, Wm,Nright):
    sum = L_an(d0,R,Wm)/2
    for i in range(1,Ndef):
        sum += L_an(bD(i,d0,dMax),R,Wm)
    sum += (Nright)* L_an(bD(Ndef, d0,dMax),R,Wm)
    return sum
    
def L_an(d , R,Wm):
    return (d + 2*( (Wm/2 - R*np.sin(theta))*np.tan(theta) + R*(1-np.cos(theta))) )

#Construction Functions
def bridge(x0, y0, d, R,Wm):
    x1 = d/2 +R*(1-np.cos(theta))
    y1 = R*np.sin(theta)
    x2 = L_an(d,R,Wm)/2
    y2 = Wm/2
    y3 = W/2
    # Geometry must be placed in cells.
    c3 = gdspy.Curve(x1+x0,y1 + y0).L(x2+x0,y2+y0, x2+x0,y3+y0, -x2+x0,y3+y0,-x2+x0,y2+y0,-x1+x0,y1+y0) #Upper Polygon
    c4 = gdspy.Curve(x1+x0,-y1+y0).L(x2+x0,-y2+y0,x2+x0,-y3+y0, -x2+x0,-y3+y0,-x2+x0,-y2+y0,-x1+x0,-y1+y0) #Lower Polygon
    p4=gdspy.Polygon(c4.get_points())
    p3 = gdspy.Polygon(c3.get_points())
    x = gdspy.boolean(gdspy.Rectangle((-x1+x0, -y1+y0), (x1+x0, y1+y0)), 
    gdspy.boolean( gdspy.Round(( d/2 + R+x0, y0), R, initial_angle=np.pi-theta,final_angle=np.pi+theta,tolerance=.0001) ,
    gdspy.Round(( -d/2 - R+x0, y0), R, initial_angle=-theta,final_angle=theta,tolerance = .0001),"or" ),"not") #Bridge Construction
    return gdspy.boolean(gdspy.boolean(p4,p3,"or"),x,"or")
    return p4
def cavity(d0,dMax,R,Wm, Nright):
    a = bridge(0, 0, d0, R,Wm)
    x0=0
    for i in range(1,Ndef+1):
        x0 += (L_an(bD(i-1,d0,dMax),R,Wm) + L_an(bD(i,d0,dMax),R,Wm))/2
        a= gdspy.boolean(a, bridge(x0,0, bD(i,d0,dMax) ,R,Wm),"or" )
        a= gdspy.boolean(a, bridge(-x0,0, bD(i,d0,dMax) ,R,Wm),"or" )
    xR = x0
    xL = x0
    for i in range(Nright):
        xR += L_an(bD(Ndef,d0,dMax),R,Wm)
        a= gdspy.boolean(a, bridge(xR,0, bD(Ndef,d0,dMax) ,R,Wm),"or" )
    for i in range(Nleft):
        xL += L_an(bD(Ndef,d0,dMax),R,Wm)
        a= gdspy.boolean(a, bridge(-xL,0, bD(Ndef,d0,dMax) ,R,Wm),"or" )
    return a
def taper(d0,dMax,R,Nright):
    xc = lenRight(d0,dMax,R,Nright)
    yc = W/2
    x1 = xc + L_R
    y1 = yc
    x2 = x1+(L_O- L_T)/2
    y2 = y1+ (W_T-W)/2
    x3 = x2+0
    y3 = W_g

    x4 = x3+ L_T
    y4 = y3
    x5 = x4
    y5 = y2
    x6 = x5+ (L_O- L_T)/2
    y6 = y5-(W_T-W)/2
    x7 = x6 + 30
    y7 = y6-W/2 +.03
    x8=x7
    y8 = y7-.06
    x9 = x8-30
    y9 = y8-W/2 +.03
    x10 = x9-(L_O- L_T)/2
    y10 = y9-(W_T-W)/2
    x11 = x10
    y11= -W_g
    x12 = x11 - L_T
    y12 = y11
    x13 = x12
    y13 = y10
    x14=x1
    y14 = y9
    x15 = xc
    y15 =yc-W


    c = gdspy.Curve(xc,yc).L(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,y15)
    return gdspy.Polygon(c.get_points())


def taper2(xc):
    yc = W/2
    x1 = xc+ L_R
    y1 = yc
    x2 = x1+ (L_O-L_T)/2
    y2 = y1+ (W_T-W)/2
    x3 = x2
    y3 = y2+ W_g-(W_T+W_Tri)/2
    x4 = x3-(L_Tri - L_T)/2
    y4 = y3+W_Tri
    x5 = x4+L_Tri
    y5 = y4
    x6 = x5- (L_Tri - L_T)/2
    y6 = y5-W_Tri
    x7 = x6
    y7 =y6-W_g+(W_T+W_Tri)/2
    x8 = x7+ (L_O-L_T)/2
    y8 = y7 - (W_T-W)/2
    x9 = x8+30
    y9 = .03
    x10 = x9
    y10 = -y9
    x11 = x8
    y11 = -y8
    x12 = x7
    y12 = -y7
    x13 = x6
    y13=-y6
    x14 = x5
    y14 = -y5
    x15 = x4
    y15 = -y4
    x16 = x3
    y16 = -y3
    x17 = x2
    y17=-y2
    x18 = x1
    y18=-y1
    x19 = xc
    y19 = -yc
    c = gdspy.Curve(xc,yc).L(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9
    ,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,y15,x16,y16,x17,y17,x18,y18,x19,y19)
    return gdspy.Polygon(c.get_points())



def device(xc,yc,d0,dMax,R,Wm,Nright):
    tL = lenLeft(d0,dMax,R,Wm) + L_L
    a = gdspy.boolean(cavity(d0,dMax,R,Wm,Nright),   gdspy.Rectangle(( -tL ,-W/2  ),(-lenLeft(d0,dMax,R,Wm),W/2   ) ),"or")
    xc = lenRight(d0,dMax,R,Wm,Nright)
    b = taper2(xc)
    c= gdspy.boolean(a,b,"or")
    c.translate(xc,yc)
    return c


def deviceInverse(xc,yc,d0,dMax,R,Nright,rot = np.pi/4):
    tLeft = lenLeft(d0,dMax,R)+L_L
    tRight = lenRight(d0,dMax,R,Nright)+L_R+L_O+30+10
    
    a= device(0,0,d0,dMax,R,Nright)
    rect = gdspy.Rectangle((-tLeft,-W_g), (tRight,W_g))
    c= gdspy.boolean(rect,a, "not")
    c.rotate(rot)
    c.translate(xc,yc)
    c = gdspy.boolean(c,text(xc,yc,d0,R,Nright),"or")
    return c
def text(xc,yc, d0,R,Nright,rot = np.pi/4):
    text = gdspy.Text(str(int(d0*1000))+ "_" + str(int(R*1000)) + "_" + str(Nright) , 6, (0, 0))
    text.rotate(rot)
    text.translate(xc,yc+13)
    return text

def text2(xc,yc, z,rot = np.pi/4):
    text = gdspy.Text(str(z) , 6, (0, 0))
    text.rotate(rot)
    text.translate(xc,yc+13)
    return text
def deviceInverse2(xc,yc,d0,dMax,R,Wm, Nright,z,rot = np.pi/4):
    tLeft = lenLeft(d0,dMax,R,Wm)
    tRight = lenRight(d0,dMax,R,Wm,Nright)+L_R+L_O+30+10
    print(tLeft)
    a= device(0,0,d0,dMax,R,Wm,Nright)
    rect = gdspy.Rectangle((-tLeft,-W_g), (tRight,W_g))
    c= gdspy.boolean(rect,a, "not")
    c.rotate(rot)
    c.translate(xc,yc)
    c = gdspy.boolean(c,text2(xc,yc,z),"or")
    return c
# z=0
# a= deviceInverse2(0,0,.1,.16,.1,20,0)
# cell.add(a)

# a= bridge(.05,.1,.1,.1)

# cell.add(a)

dArr = np.linspace(.05,.15,5)
dM = [.0707, .095,.131, .16,.191]
Wm = .24
rArr = np.linspace(.1,.2,5)

xStep = 40
yStep = 45
xOff = 20
coaX = 240
coaY= 240
z=0
#For rotated smample
# for x in [0,1]:
#     for y in [0,1]:
#        for i in range(len(dArr)):
#             for j in range(len(dArr)):
#                 if (j%2 == 0):
#                     a = deviceInverse2(xStep*i+coaX*x+xOff,yStep*j+coaY*y, dArr[i],dM[i],rArr[j],Wm,5*(x+1),z)
#                     cell.add(a)
#                     print(x,y,i,j)
#                     z+=1
#                 else:
#                     a = deviceInverse2(xStep*i+coaX*x,yStep*j+coaY*y, dArr[i],dM[i],rArr[j],Wm,5*(x+1),z)
#                     cell.add(a)
#                     print(x,y,i,j)
#                     z+=1
a= device(0,0,dArr[0],dM[0],rArr[0],Wm,10)            
# a = bridge(0,0, dArr[0],rArr[0])
cell.add(a)
end = time.time()
print(end - start)

#Save the library in a file called 'first.gds'.
# lib.write_gds('chip2_layout.gds')
# Optionally, save an image of the cell as SVG.
#cell.write_svg('vis_OMC.svg')

# Display all cells using the internal viewer.
gdspy.LayoutViewer()