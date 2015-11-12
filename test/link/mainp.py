#!/usr/bin/env python

import numpy as np
from plumed import *
import sys

n=4
step=40
timestep=0.001
masses=np.ones(n)
charges=np.zeros(n)
box=np.diag(10*np.ones(3))
pos=np.zeros([n,3])
pos[:,0]=np.arange(n)
forces=np.zeros([n,3])
virial=np.zeros([3,3])
temp=3.


plumed = Plumed()

plumed.cmd("setMDEngine","none")
plumed.cmd("setPlumedDat","plumed.dat")
#plumed.cmd("setLogFile","plumed.log")
plumed.cmd("setNatoms",n)
plumed.cmd("setTimestep",timestep)
plumed.cmd("setKbT",temp)

plumed.cmd("init")

plumed.cmd("setStep",step)
# Write wrapper for actions. Wrappers to get values back from actions.
plumed.cmd("action","DISTANCE LABEL=l6 ATOMS=3,4")
plumed.cmd("action","PRINT ARG=l6 FILE=COLVAR2 STRIDE=1")

plumed.cmd("setMasses",masses)

# To do. Perhaps have the two options.
#plumed.setMasses(masses)

plumed.cmd("setCharges",charges)
plumed.cmd("setBox",box)
plumed.cmd("setPositions",pos)
plumed.cmd("setForces",forces)
plumed.cmd("setVirial",virial)

plumed.cmd("calc")

# Retrieve results of calculation directly in plumed. Check how forces, etc. are passed back to the md code. Remember that the bias is passed back to gromacs for replica exchange so we can also passed it to python. Pass things back according to a label. 
# Check how to delete actions.
