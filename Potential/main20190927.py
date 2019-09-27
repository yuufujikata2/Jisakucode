import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from surfaceintegral import surfaceintegral
from call_rseq import rseq
from makepotential import makepotential

EPSVAL = 1.e-20

def main():

    # make environment
    region = 1
    nr = 201
    gridpx = 100
    gridpy = 100
    gridpz = 100
    x,y,z = grid(gridpx,gridpy,gridpz,region)


    #make mesh

    #linear mesh
    #r =np.linspace(0.0001,region,nr)

    #log mesh
    a = np.log(2) / (nr - 1) 
    b = region / (np.e**(a * ( nr - 1)) - 1)
    r = np.array([b * (np.e**(a * i) - 1) for i in range(nr)])

  
    # make potential
    V = makepotential(x,y,z,pottype="cubic",potbottom=-1,potshow_f=False)

    # surface integral
    V_radial = surfaceintegral(x,y,z,r,V,method="rasen",potshow_f=False)

    # solve radial shreduinger equation
    energy = 0.
    emin = -10.
    emax = 100.
    gtry = np.zeros((2 * nr),dtype = np.float64)
    gfac = np.zeros((nr),dtype = np.float64)
    lvalsh = 0 # for s orbital
    node = 1 # for minimum energy orbital
    nre = 0
    slo = -1.  # for close
    val = EPSVAL
    vofi = np.array (V_radial)  # select method of surface integral 
    iz = 0

    energy, gtry, gfac, nre, slo, val = rseq(a, b, energy, emin, emax, gtry, gfac, lvalsh, node, nr, nre, r, slo, vofi, val, iz)

    print ("E = ", energy)
    print ("slope = ", slo)
    print ("value = ", val)
    with open ("wavefunc.dat", mode = "w") as fw_w :
        fw_w.write("# rfoi, wavefunction\n")
        for i in range (nr):
            fw_w.write("{:>13.8f}".format(r[i]))
            fw_w.write("{:>13.8f}".format(gtry[i]))
            fw_w.write("\n")

    energy = 0.
    emin = -10.
    emax = 100.
    gtry2 = np.zeros((2 * nr),dtype = np.float64)
    gfac = np.zeros((nr),dtype = np.float64)
    lvalsh = 0 # for s orbital
    node = 0 # for minimum energy orbital
    nre = 0
    slo = -1.  # for close
    val = EPSVAL
    vofi = np.array (V_radial)  # select method of surface integral 
    iz = 0

    energy, gtry2, gfac, nre, slo, val = rseq(a, b, energy, emin, emax, gtry2, gfac, lvalsh, node, nr, nre, r, slo, vofi, val, iz)

    print ("E = ", energy)
    print ("slope = ", slo)
    print ("value = ", val)
    with open ("wavefunc2.dat", mode = "w") as fw_w :
        fw_w.write("# rfoi, wavefunction\n")
        for i in range (nr):
            fw_w.write("{:>13.8f}".format(r[i]))
            fw_w.write("{:>13.8f}".format(gtry[i] * gtry2[i]))
            fw_w.write("\n")


    gall = 0.

    for i in range (nr):
        gall += gtry[i] * gtry2[i]
    print ("gall = ",gall)




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
