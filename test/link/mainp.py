#!/usr/bin/env python

import numpy as np
from plumed import *
import sys

n=4
step=40
masses=np.ones(n)
charges=np.zeros(n)
box=np.diag(10*np.ones(3))
pos=np.zeros([n,3])
pos[:,0]=np.arange(n)
forces=np.zeros([n,3])
virial=np.zeros([3,3])


plumed = Plumed()

plumed.cmd("setMDEngine","none")
plumed.cmd("setPlumedDat","plumed.dat")
#plumed.cmd("setLogFile","plumed.log")
plumed.cmd("setNatoms",n)
plumed.cmd("init")

plumed.cmd("action","DISTANCE LABEL=l6 ATOMS=3,4")
plumed.cmd("action","PRINT ARG=l6 FILE=COLVAR2 STRIDE=1")
plumed.cmd("setStep",step)
plumed.cmd("setMasses",masses)
plumed.cmd("setCharges",charges)
plumed.cmd("setBox",box)
plumed.cmd("setPositions",pos)
plumed.cmd("setForces",forces)
plumed.cmd("setVirial",virial)

plumed.cmd("calc")
