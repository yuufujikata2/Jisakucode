import ctypes
import numpy as np
import sys

def inv_remat(n,mat):
    path = sys.path[0]+"/Inv/inv_remat.so"
    f = ctypes.CDLL(path)
    f.inv_remat.argtypes = [ ctypes.c_int32,
                          np.ctypeslib.ndpointer(dtype=np.float64)
                      ]

    f.inv_remat.restype = ctypes.c_int32
    f.inv_remat(ctypes.c_int32(n),mat)
    return  

def main():
    a = 0.
    b = 0.
    e = 0.
    eb1 = 0.
    eb2 = 0.
    g = np.zeros((2))
    gfac = np.zeros((2))
    l = 0
    nod = 0
    nr = 0
    nre = 0
    nrep = 0
    rofi = np.zeros((2))
    slo = 0.
    v = np.zeros((2))
    val = 0.
    iz = 0
    rseq(a, b, e, eb1, eb2, g, gfac, l, nod, nr, nre, rofi, slo, v, val, iz)


if __name__ == "__main__" :
    main()
