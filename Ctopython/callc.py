import ctypes
import numpy as np
import sys
def call_c(xp,yp):
    path = sys.path[0]+"/test.so"
    c = ctypes.CDLL(path)
    c.add.argtypes = [ctypes.POINTER(ctypes.c_int32),
                       ctypes.POINTER(ctypes.c_int32)]
    c.add.restype = ctypes.c_int32

    result = c.add(xp,yp)

    return result

def main ():
    x = ctypes.c_int32(3)
    y = ctypes.c_int32(4)
    xp = ctypes.pointer(x)
    xpp = ctypes.pointer(xp)
    yp = ctypes.pointer(y)
    result  = call_c(xp,yp)

    print ("x = ",x.value,"y = ", y.value ,"result = " ,result, "xpp = ",xpp[0][0],"y2 = ",yp[0])



if __name__ == "__main__":
    main()
