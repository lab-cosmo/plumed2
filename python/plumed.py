import ctypes as ct

import numpy

_libplumed = ct.CDLL('libplumed.so')

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
        elif (type(val) is long):
            value=ct.byref(ct.c_double(float(val)))
        elif (type(val) is numpy.ndarray):
            # What if there is integer array?
            # Make linear array
            val=val.reshape(-1)
            value=val.ctypes.data
        else:
            value=val
        # Add the possibility to input tuples or lists.
       _libplumed.plumed_cmd(self.p, key, value)
