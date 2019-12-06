import ctypes
import numpy as np
import sys

def besneu(modified,x,lvalsh):
    path = sys.path[0]+"/Besneu/besneu.so"
    f = ctypes.CDLL(path)
    f.besneu.argtypes = [ ctypes.c_int32,
                        ctypes.c_double,
                        ctypes.c_int32,
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                      ]

    f.besneu.restype = ctypes.c_void_p
    modified = ctypes.c_int32(modified)
    x = ctypes.c_double(x)
    lvalsh = ctypes.c_int32(lvalsh)
    jl = ctypes.c_double(0)
    nl = ctypes.c_double(0)
    jlp = ctypes.c_double(0)
    nlp = ctypes.c_double(0)
    f.besneu(modified,x,lvalsh,ctypes.byref(jl),ctypes.byref(nl),ctypes.byref(jlp),ctypes.byref(nlp))
    jl = jl.value
    nl = nl.value
    jlp = jlp.value
    nlp = nlp.value
    return  jl,nl,jlp,nlp

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
