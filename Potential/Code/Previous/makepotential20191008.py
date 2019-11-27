import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def makepotential(x,y,z,pottype="cubic",potbottom=-1.,potshow_f=False):
    if pottype =="cubic" :
        V = cubic(x,y,z,potbottom)
    elif pottype == "cylinder" :
        V = cylinder(x,y,z,potbottom)
    elif pottype == "flat" :
        V = flat(x,y,z,potbottom)
    else :
        print("error: There is no pottype entered")
        sys.exit()

    if potshow_f :
        potentialshow(x,y,z,V,potbottom)

    return V

def flat(x,y,z,a):
    V = np.zeros((len(x),len(y),len(z)))
    V += a
    return V

def cubic(x,y,z,a):
    xx, yy, zz = np.meshgrid(x,y,z)
    V = np.where( (abs(xx) <= x[-1] / 2) & (abs(yy) <= y[-1] / 2) & (abs(zz) <= z[-1] / 2), a,0.)
    return V

def cylinder(x,y,z,a):
    xx, yy, zz = np.meshgrid(x,y,z)
    rr = np.sqrt(xx**2 + yy**2)
    V = np.where((rr <= 0.3 * x[-1]) & (abs(zz) <= 0.7 * z[-1]), a, 0.)
    return V

def potentialshow(x,y,z,V,a):
    xx, yy, zz = np.meshgrid(x,y,z)
    xxindex = np.where(V == a )
    xx2 = xx[xxindex]
    yyindex = np.where(V == a )
    yy2 = yy[yyindex]
    zzindex = np.where(V == a )
    zz2 = zz[zzindex]

    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    ax.scatter(xx2,yy2,zz2)
    plt.show()



