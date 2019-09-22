import ctypes
import numpy as np

def call_fortran(A):
  B=np.array(A,dtype=np.int32)
  BT=B.flatten("F")
  f=np.ctypeslib.load_library("libfort2.so",".")
  f.rensashirabe_.argtypes=[np.ctypeslib.ndpointer(dtype=np.int32)]
  f.rensashirabe_.restype=ctypes.c_int32
  a=f.rensashirabe_(BT)
  return a





