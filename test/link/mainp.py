#!/usr/bin/env python

import numpy as np
from plumed import *
from tools import *
import sys

# Read XYZ
filename= "methane.xyz"
atom_type, pos = read_xyz(filename)
num_atoms=pos.shape[0]
# Define a few things
step=40
timestep=0.001
masses=np.ones(num_atoms)
charges=np.zeros(num_atoms)
box=np.diag(10*np.ones(3))
forces=np.zeros([num_atoms,3])
virial=np.zeros([3,3])
temp=3.

# Intialize class instance
plumed = Plumed()

# Pre-init settings
plumed.cmd("setMDEngine","none")
plumed.cmd("setPlumedDat","")
#plumed.cmd("setLogFile","plumed.log")
plumed.cmd("setNatoms",num_atoms)
plumed.cmd("setTimestep",timestep)
plumed.cmd("setKbT",temp)

# Init
plumed.cmd("init")

# Post-init settings
plumed.cmd("setStep",step)
# This is new. I had to change PlumedMain.cpp to do this.
plumed.cmd("action","DISTANCE LABEL=l6 ATOMS=3,4")
plumed.cmd("action","PRINT ARG=l6 FILE=COLVAR STRIDE=1")
plumed.cmd("setMasses",masses)
plumed.cmd("setCharges",charges)
plumed.cmd("setBox",box)
plumed.cmd("setPositions",pos)
plumed.cmd("setForces",forces)
plumed.cmd("setVirial",virial)

# Calculate
plumed.cmd("calc")

# TO DO: Write wrapper for actions. For instance "DISTANCE LABEL=l6 ATOMS=3,4" = Distance(label="l6",atoms)
# TO DO: Wrappers to get values back from actions. Very important!
# TO DO: Perhaps have two options. plumed.setMasses(masses) = plumed.cmd("setMasses",masses) . The first would be more pythonic...
# IDEAS AND TO DO: Retrieve results of calculation directly in plumed. Check how forces, etc. are passed back to the md code (I think that they are pointers). Remember that the bias is passed back to gromacs for replica exchange so we can also passed it to python (perhaps do the same for CVs?). Pass things back according to a label. Check how to delete actions.
