"""
To call fortran subroutine

"""
import ctypes
import numpy as np
import time
from cProfile import Profile
import sys

def lebedev(num,x,y,z,w):
    path = sys.path[0]+"/lebedev.so"
    f = ctypes.CDLL(path)
    num_c = ctypes.byref(ctypes.c_int32(num)) 
    f.ld_by_order_.argtypes = [     ctypes.POINTER(ctypes.c_int32),   
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    np.ctypeslib.ndpointer(dtype=np.float64)
                                    ]
    f.ld_by_order_.restype = ctypes.c_void_p
    f.ld_by_order_(num_c,x,y,z,w)
    return 0 

def main():
    num = 2030
    x = np.zeros(num)
    y = np.zeros(num)
    z = np.zeros(num)
    w = np.zeros(num)
    call_fortran(num,x,y,z,w)

if __name__=="__main__":
  main()
