import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab

def makepotential(xx,yy,zz,pot_region,pot_type="cubic",pot_bottom=-1.,pot_show_f=False):
    if pot_type =="cubic" :
        V = cubic(xx,yy,zz,pot_region,pot_bottom)
    elif pot_type == "cylinder" :
        V = cylinder(xx,yy,zz,pot_region,pot_bottom)
    elif pot_type == "flat" :
        V = flat(xx,yy,zz,pot_bottom)
    elif pot_type == "sinproduct" :
        V = sinproduct(xx,yy,zz,pot_bottom)
    elif pot_type == "cosproduct" :
        V = cosproduct(xx,yy,zz,pot_bottom)
    else :
        print("error: There is no pottype entered")
        sys.exit()

    if pot_show_f :
        if len(xx) * len(yy) * len(zz) > 500000 :
            print("error: memory over")
            sys.exit()
        potentialshow(xx,yy,zz,V,pot_bottom)

    return V

def flat(xx,yy,zz,a):
    V = np.zeros((len(xx),len(yy),len(zz)))
    V += a
    return V

def cubic(xx,yy,zz,pot_region,a):
    V = np.where( (abs(xx) <= pot_region[0]) & (abs(yy) <= pot_region[1]) & (abs(zz) <= pot_region[2]), a,0.)
    return V

def cylinder(xx,yy,zz,pot_region,a):
    rr = np.sqrt(xx**2 + yy**2)
    V = np.where((rr <= 0.3 * pot_region[0]) & (abs(zz) <= 0.7 * pot_region[2]), a, 0.)
    return V

def sinproduct(xx,yy,zz,a):
    V = np.sin(xx) * np.sin(yy) * np.sin(zz) * a
    return V

def cosproduct(xx,yy,zz,a):
    V = np.cos(xx-0.5) * np.cos(yy) * np.cos(zz) * a
    return V

def potentialshow(xx,yy,zz,V,a):
    xxindex = np.where(V == a )
    xx2 = xx[xxindex]
    yyindex = np.where(V == a )
    yy2 = yy[yyindex]
    zzindex = np.where(V == a )
    zz2 = zz[zzindex]
    xxindex0 = np.where(V == 0. ) 
    xx0 = xx[xxindex0]
    yyindex0 = np.where(V == 0. )
    yy0 = yy[yyindex0]
    zzindex0 = np.where(V == 0. )
    zz0 = zz[zzindex0]

    """
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    ax.plot(xx2,yy2,zz2,"o",color = "blue")
    ax.plot(xx0[::2][::2][::2],yy0[::2][::2][::2],zz0[::2][::2][::2],"o",color = "green",ms = 2, mew = 0.1)
    plt.show()
    """

    mlab.points3d(V,scale_factor= 0.4)
    mlab.show()




