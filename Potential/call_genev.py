import ctypes
import numpy as np
import sys


def solve_genev(n,lambda_mat,qmetric_mat,alphalong,betalong,revec):
    path = sys.path[0] + "/LAPACK/lapack.so"
    f = ctypes.CDLL(path)
    JOBVL = ctypes.c_wchar_p("N")
    JOBVR = ctypes.c_wchar_p("V")
    LDVL = 1
    INFO = ctypes.c_int32(0)
    ALPHAI = np.zeros(n)
    VL = np.zeros(LDVL * n)
    LWORKFAC = 40
    LWORK = n * LWORKFAC
    WORK = np.zeros(LWORK)
    LWORK = ctypes.c_int32(n * LWORKFAC)
    n = ctypes.c_int32(n)
    LDVL = ctypes.c_int32(LDVL)

    f.dggev_.argtypes = [ ctypes.POINTER(ctypes.c_char),
                        ctypes.POINTER(ctypes.c_char),
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_int32),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        ctypes.POINTER(ctypes.c_int32),
                        ctypes.POINTER(ctypes.c_int32)
                      ]

    f.dggev_.restype = ctypes.c_void_p
    f.dggev_(JOBVL, JOBVR, ctypes.byref(n), lambda_mat, ctypes.byref(n), qmetric_mat, ctypes.byref(n), alphalong, ALPHAI, betalong, VL, ctypes.byref(LDVL), revec, ctypes.byref(n), WORK, ctypes.byref(LWORK), ctypes.byref(INFO) ) 
    info = info.value
    return  info

if __name__ == "__main__" :
    main()
