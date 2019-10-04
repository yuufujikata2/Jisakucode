import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from surfaceintegral import surfaceintegral
from call_rseq import rseq
from makepotential import makepotential
from basis import Basis
from integrate import integrate
from scipy import interpolate
from mayavi import mlab

EPSVAL = 1.e-20

def main():

    # make environment
    region = 2
    nr = 201
    gridpx = 100
    gridpy = 100
    gridpz = 100
    x,y,z = grid(gridpx,gridpy,gridpz,region)
    xx, yy, zz = np.meshgrid(x,y,z)

    #make mesh

    #linear mesh
    #r =np.linspace(0.0001,region,nr)

    #log mesh
    a = np.log(2) / (nr - 1) 
    b = region / (np.e**(a * ( nr - 1)) - 1)
    rofi = np.array([b * (np.e**(a * i) - 1) for i in range(nr)])

  
    # make potential
    V = makepotential(x,y,z,pottype="cubic",potbottom=-1,potshow_f=False)

    # surface integral
    V_radial = surfaceintegral(x,y,z,rofi,V,method="lebedev_py",potshow_f=False)
    vofi = np.array (V_radial)  # select method of surface integral

    # make basis
    
    node_open = 2
    node_close = 2
    LMAX = 3

    all_basis = []

    for lvalsh in range (LMAX):
        l_basis = []
        # for open channel
        val = 1.
        slo = 0.
        for node in range(node_open):
            basis = Basis(nr)
            emin = -10.
            emax = 100.
            basis.make_basis(a,b,emin,emax,lvalsh,node,nr,rofi,slo,vofi,val)
            l_basis.append(basis)

        # for close channel
        val = 0.
        slo = -1.
        for node in range(node_close):
            basis = Basis(nr)
            emin = -10.
            emax = 100.
            basis.make_basis(a,b,emin,emax,lvalsh,node,nr,rofi,slo,vofi,val)
            l_basis.append(basis)

        all_basis.append(l_basis)

    with open ("wavefunc.dat", mode = "w") as fw_w :
        for i in range(nr):
            fw_w.write("{:>13.8f}".format(rofi[i]))
            for l_basis in all_basis:
                for nyu_basis in l_basis:
                    fw_w.write("{:>13.8f}".format(nyu_basis.g[i]))
            fw_w.write("\n")

    hsmat = np.zeros((LMAX,LMAX,node_open + node_close,node_open + node_close),dtype = np.float64)
    lmat = np.zeros((LMAX,LMAX,node_open + node_close,node_open + node_close), dtype = np.float64)
    qmat = np.zeros((LMAX,LMAX,node_open + node_close,node_open + node_close), dtype = np.float64)

    for l1 in range (LMAX):
        for l2 in range (LMAX):
            if l1 != l2 :
                continue
            for n1 in range (node_open + node_close):
                for n2 in range (node_open + node_close):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                    hsmat[l1][l2][n1][n2] = integrate(all_basis[l1][n1].g[:nr] * all_basis[l2][n2].g[:nr],rofi,nr) * all_basis[l1][n1].e
                    lmat[l1][l2][n1][n2] = all_basis[l1][n1].val * all_basis[l2][n2].slo
                    qmat[l1][l2][n1][n2] = all_basis[l1][n1].val * all_basis[l2][n2].val
    print ("\nhsmat")
    print (hsmat)
    print ("\nlmat")
    print (lmat)
    print ("\nqmat")
    print (qmat)


    #make not spherical potential
    my_radial_interfanc = interpolate.interp1d(rofi, V_radial)
    V_ang = np.where(np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] , V - my_radial_interfanc(np.sqrt(xx * xx + yy * yy + zz * zz)),0. )
    #WARING!!!!!!!!!!!!!!!!!!!!!!
    """
    Fujikata rewrote ~/.local/lib/python3.6/site-packages/scipy/interpolate/interpolate.py line 690~702
    To avoid exit with error "A value in x_new is below the interpolation range."
    """
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!

    with open ("V_ang.dat",mode = "w") as fw_a:
        for i in range(len(xx)):
            for j in range(len(yy)):
                for k in range(len(zz)):
                    fw_a.write("{:>13.8f}".format(xx[j][i][k]))
                    fw_a.write("{:>13.8f}".format(yy[j][i][k]))
                    fw_a.write("{:>13.8f}".format(zz[j][i][k]))
                    fw_a.write("{:>13.8f}\n".format(V_ang[j][i][k]))

    
    obj = mlab.volume_slice(V_ang)
    #mayavi.mlab.plot3d(xx,yy,V_ang)
    #mlab.contour3d(V_ang)
    mlab.show()
    




def grid (nx,ny,nz,region):
    x = np.linspace(-region,region,nx)
    y = np.linspace(-region,region,ny)
    z = np.linspace(-region,region,nz)
    return x,y,z


def fitting(r,V,d):
    w = np.polyfit(r,V,d)
    print(w)
    xs = np.linspace(r[0],r[-1],300)
    ys = np.polyval(w,xs)
    return xs,ys
    
if __name__ =="__main__":
    main()
