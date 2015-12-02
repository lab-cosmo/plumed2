import ctypes as ct
import numpy
from sys import platform as _platform

# Find operating system
if _platform == "linux" or _platform == "linux2":
  _libplumed = ct.CDLL('libplumed.so')
elif _platform == "darwin":
  # OS X
  _libplumed = ct.CDLL('libplumed.dylib')

class Plumed:
  def __init__(self):
    self.p = _libplumed.plumed_create()
  def __del__(self):
    _libplumed.plumed_finalize(self.p)
  def cmd(self, key, val=None):
    # Handle different data types
    if (type(val) is int):
      value=ct.byref(ct.c_int(val))
    elif (type(val) is float):
      value=ct.byref(ct.c_double(val))
    # In case someone crazy defines a long
    #elif (type(val) is long):
      #value=ct.byref(ct.c_double(float(val)))
    elif (type(val) is numpy.ndarray):
      # Make linear array
      val=val.reshape(-1)
      value=val.ctypes.data
    elif (type(val) is str):
      value=bytes(val.encode(encoding='UTF-8',errors='strict'))
    # To do: Add the possibility to input tuples or lists.
    else:
      value=val
    key=bytes(key.encode(encoding='UTF-8',errors='strict'))
    _libplumed.plumed_cmd(self.p, key, value)
