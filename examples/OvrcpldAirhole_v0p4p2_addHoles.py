"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser

import numpy as np
import os


hole_radii = np.loadtxt(os.path.join(os.path.curdir, 'examples/v0p4p2/holeStruct.txt'), dtype=float, usecols=(0,1), unpack=False)
hole_radii /= 2
lattice_constants = np.loadtxt(os.path.join(os.path.curdir, 'examples/v0p4p2/periodStruct.txt'), dtype=float, usecols=(0), unpack=False)
# Unit cells are triangular prisms
beam_width = 0.480e-6
apex_half_angle = 50*np.pi/180
lattice_constant = 0.270e-6
beam_height = beam_width / 2 * np.tan(np.pi - apex_half_angle)
# Radius of the air holes in the cells
# hole_radius = 0.5*0.120e-6
# The length of the cavity beam
beam_length = 10e-6
# The target resonance frequency, in Hz
target_frequency = 407e12

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=4, hide=False)

unit_cells = []
for ((radius_x,radius_y), lattice_constant) in zip(hole_radii,lattice_constants):
  cell_box = TriStructure(Vec3(0), Vec3(beam_width, apex_half_angle, lattice_constant), 
                         DielectricMaterial(2.4028, order=2), rot_angles=(np.pi/2, np.pi/2, 0))
  
  cell_hole = CylinderStructure(Vec3(0), beam_height, radius_x, DielectricMaterial(1, order=1), radius2=radius_y)
  unit_cells += [UnitCell(structures=[ cell_box, cell_hole ], size=Vec3(lattice_constant), engine=engine)]

cavity = Cavity1D(
  unit_cells=unit_cells,
  structures=[TriStructure(Vec3(0), Vec3(beam_width, apex_half_angle, beam_length), 
                        DielectricMaterial(2.4028, order=2), rot_angles=(np.pi/2, np.pi/2, 0))],
  engine=engine
)

parsed = DielectricExtrusionFaceGDSParser(cavity, invert=False)
parsed.show()

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("v0p4p2p2.obj")

r1 = cavity.simulate("resonance", target_freq=target_frequency)

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
  r1["freq"], r1["vmode"],
  1/(1/r1["qxmin"] + 1/r1["qxmax"]),
  1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

# r2 = cavity.simulate("quasipotential", target_freq=target_frequency)

# # Plot the quasipotential
# r2.show()

# r3 = cavity.simulate("guidedness", target_freq=target_frequency)

# print("Guidedness: %f" % r3)
