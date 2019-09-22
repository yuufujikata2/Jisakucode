"""
To call fortran subroutine

"""
import ctypes
import numpy as np
import time
from cProfile import Profile
import sys

class Fort():
  def __init__(self):
    self.path=sys.path[0]+"/ctrlsucell.so"
    self.f=ctypes.CDLL(self.path)
    #self.f=ctypes.CDLL("/home/harada/Jisakucode/ctrlclusterf.so")
    self.f.cluster_.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                             ctypes.POINTER(ctypes.c_int32),
                             np.ctypeslib.ndpointer(dtype=np.float64),
                             ctypes.POINTER(ctypes.c_double),
                             np.ctypeslib.ndpointer(dtype=np.int32),
                             ctypes.POINTER(ctypes.c_int32),
                             ctypes.POINTER(ctypes.c_int32),
                             ctypes.POINTER(ctypes.c_int32)
                             ]
    self.f.cluster_.restype=ctypes.c_void_p

  def call_fortran(self,Plat,nbas,cart,alat,alliz,sa,sb,sc):
    Pl2=np.array(Plat)
    Pl3=Pl2.flatten("F")
    nb2=ctypes.byref(ctypes.c_int32(nbas))
    ca2=np.array(cart)
    ca3=ca2.flatten("F")
    al2=ctypes.byref(ctypes.c_double(alat))
    iz2=alliz.flatten("F")
    sa2=ctypes.byref(ctypes.c_int32(sa))
    sb2=ctypes.byref(ctypes.c_int32(sb))
    sc2=ctypes.byref(ctypes.c_int32(sc))
    self.f.cluster_(Pl3,nb2,ca3,al2,iz2,sa2,sb2,sc2)
    return 0 

def main():
  clradius=6.
  Plat=np.array([[1., 0., -1.],[0.,2., 0.],[0.,0.,1.]])
  origincart=np.array([0. , 0. ,0.])
  fort=Fort()
  fort.call_fortran(clradius,Plat,origincart)

if __name__=="__main__":
  main()
