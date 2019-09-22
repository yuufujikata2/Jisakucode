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
    self.path=sys.path[0]+"/ctrlxpdcluster.so"
    self.f=ctypes.CDLL(self.path)
    #self.f=ctypes.CDLL("/home/harada/Jisakucode/ctrlclusterf.so")
    self.f.cluster_.argtypes=[ctypes.POINTER(ctypes.c_float),
                             np.ctypeslib.ndpointer(dtype=np.float32),
                             ctypes.POINTER(ctypes.c_int32),
                             np.ctypeslib.ndpointer(dtype=np.float32),
                             ctypes.POINTER(ctypes.c_float),
                             np.ctypeslib.ndpointer(dtype=np.int32),
                             np.ctypeslib.ndpointer(dtype=np.float32)
                             ]
    self.f.cluster_.restype=ctypes.c_void_p

  def call_fortran(self,clradius,Plat,nbas,cart,alat,alliz,miller):
    Pl2=np.array(Plat ,dtype = np.float32)
    Pl3=Pl2.flatten("F")
    cl2=ctypes.byref(ctypes.c_float(clradius))
    nb2=ctypes.byref(ctypes.c_int32(nbas))
    ca2=np.array(cart,dtype = np.float32)
    ca3=ca2.flatten("F")
    al2=ctypes.byref(ctypes.c_float(alat))
    iz2=alliz.flatten("F")
    miller2=np.array(miller,dtype = np.float32)
    miller3=miller2.flatten("F")
    self.f.cluster_(cl2,Pl3,nb2,ca3,al2,iz2,miller3)
    return 0 

def main():
  clradius=6.
  Plat=np.array([[1., 0., -1.],[0.,2., 0.],[0.,0.,1.]])
  origincart=np.array([0. , 0. ,0.])
  fort=Fort()
  fort.call_fortran(clradius,Plat,origincart)

if __name__=="__main__":
  main()
