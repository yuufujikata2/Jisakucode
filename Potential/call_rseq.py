import ctypes
import numpy as np
import sys

def rseq(a, b, e, eb1, eb2, g, gfac, l, nod, nr, nre, rofi, slo, v, val, iz):
    path = sys.path[0]+"/rseq.so"
    f = ctypes.CDLL(path)
    f.rseq.argtypes = [ ctypes.c_double,
                        ctypes.c_double,
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.c_double,
                        ctypes.c_double,
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.c_int32,
                        ctypes.c_int32,
                        ctypes.c_int32,
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_double),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.c_int32
                      ]

    f.rseq.restype = ctypes.c_void_p
    e = ctypes.c_double(e)
    nre = ctypes.c_int32(nre)
    val = ctypes.c_double(val)
    slo = ctypes.c_double(slo)
    f.rseq(a, b, ctypes.byref(e), eb1, eb2, g, gfac, l, nod, nr, ctypes.byref(nre), rofi, ctypes.byref(slo), v, ctypes.byref(val), iz)
    e = e.value
    nre = nre.value
    val = val.value
    slo = slo.value
    return a, b, e, eb1, eb2, g, gfac, l, nod, nr, nre, rofi, slo, v, val, iz

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
