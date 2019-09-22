import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RegularGridInterpolator
import time
import rasen
import sys
from surface_int_rasen import SI_rasen
from surface_int_lebedev_f import SI_lebedev_f
from surface_int_lebedev_py import SI_lebedev_py
from call_rseq import rseq

EPSVAL = 1.e-20

def main():
    # file open
    fw_r = open("v_r.dat",mode="w")
    fw_lf = open("v_lf.dat",mode="w")
    fw_lp = open("v_lp.dat",mode="w")


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
    r = np.array([b * (np.e**(a * ( i )) - 1) for i in range(nr )])

  
    # make potential
    V = cubic(x,y,z)
    #V = cylinder(x,y,z)
    #V = flat(x,y,z,-1)
    #potentialshow(x,y,z,V)
    #sys.exit()


    # make interpolate function
    my_inter_func = RegularGridInterpolator((x, y, z), V)


    # make surface integral object
    surface_integral_r = SI_rasen(5000)
    surface_integral_lf = SI_lebedev_f(-1)
    surface_integral_lp = SI_lebedev_py()


    # array to contain result of surface integral
    V_all_r = []
    V_all_lf = []
    V_all_lp = []


    # do surface integral 
    for dr in r:

        """
        newV = np.zeros(npoints)
        for i in range(npoints):
            if abs(coord[i][0]) <= 0.5 and abs(coord[i][1]) <= 0.5 and abs(coord[i][2]) <= 0.5:
                newV[i] = 1.
            else:
                newV[i] = 0.
        """

        integral_result_r = surface_integral_r.integral(my_inter_func, dr)
        integral_result_lf = surface_integral_lf.integral(my_inter_func, dr)
        integral_result_lp = surface_integral_lp.integral(my_inter_func, dr)
        V_all_r.append(integral_result_r)
        V_all_lf.append(integral_result_lf)
        V_all_lp.append(integral_result_lp)
        fw_r.write("{:>11.6f}".format(dr))
        fw_r.write("{:>11.6f}".format(integral_result_r))
        fw_r.write("\n")
        fw_lf.write("{:>11.6f}".format(dr))
        fw_lf.write("{:>11.6f}".format(integral_result_lf))
        fw_lf.write("\n")
        fw_lp.write("{:>11.6f}".format(dr))
        fw_lp.write("{:>11.6f}".format(integral_result_lp))
        fw_lp.write("\n")
        
                
    fw_r.close()
    fw_lf.close()
    fw_lp.close()
    #plt.vlines(np.sqrt(2) / 2, 0,1)
    #plt.vlines(np.sqrt(3) / 2, 0,1)
    plt.plot(r,V_all_r,marker=".")
    plt.plot(r,V_all_lf,marker=".")
    plt.plot(r,V_all_lp,marker=".")
    #xs,ys = fitting(r,V_all,2)
    #plt.plot(xs,ys)
    plt.show()


    # solve radial shreduinger equation
    energy = 0.
    emin = -10.
    emax = 100.
    gtry = np.zeros((2 * nr),dtype = np.float64)
    gfac = np.zeros((nr),dtype = np.float64)
    lvalsh = 1 # for s orbital
    node = 1 # for minimum energy orbital
    nre = 0
    slo = -1.  # for close
    val = EPSVAL
    vofi = np.array (V_all_r)  # select method of surface integral 
    iz = 0

    a, b, energy, emin, emax, gtry, gfac, lvalsh, node, nr, nre, r, slo, vofi, val, iz = rseq(a, b, energy, emin, emax, gtry, gfac, lvalsh, node, nr, nre, r, slo, vofi, val, iz)

    print ("E = ", energy)
    print ("slop = ", slo)
    print ("val = ", val)
    with open ("wavefunc.dat", mode = "w") as fw_w :
        fw_w.write("# rfoi, wavefunction\n")
        for i in range (nr):
            fw_w.write("{:>13.8f}".format(r[i]))
            fw_w.write("{:>13.8f}".format(gtry[i]))
            fw_w.write("\n")




def grid (nx,ny,nz,region):
    x = np.linspace(-region,region,nx)
    y = np.linspace(-region,region,ny)
    z = np.linspace(-region,region,nz)
    return x,y,z

def flat(x,y,z,a):
    V = np.zeros((len(x),len(y),len(z)))
    V += a
    return V

def cubic(x,y,z):
    xx, yy, zz = np.meshgrid(x,y,z)
    V = np.where( (abs(xx) <= 0.5) & (abs(yy) <= 0.5) & (abs(zz) <= 0.5), -1.,0.)
    return V

def cylinder(x,y,z):
    xx, yy, zz = np.meshgrid(x,y,z)
    rr = np.sqrt(xx**2 + yy**2)
    V = np.where((rr <= 0.3) & (abs(zz) <= 0.7), -1., 0.)
    return V

def fitting(r,V,d):
    w = np.polyfit(r,V,d)
    print(w)
    xs = np.linspace(r[0],r[-1],300)
    ys = np.polyval(w,xs)
    return xs,ys

def potentialshow(x,y,z,V):
    xx, yy, zz = np.meshgrid(x,y,z)
    xxindex = np.where(V == 1. )
    xx2 = xx[xxindex]
    yyindex = np.where(V == 1. )
    yy2 = yy[yyindex]
    zzindex = np.where(V == 1. )
    zz2 = zz[zzindex]

    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    ax.scatter(xx2,yy2,zz2)
    plt.show()
    
def polargrid(dr):
    ntheta = 100
    nphi = 200
    coordx = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    coordy = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    coordz = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    number = 0
    for phi in range (nphi):
        if phi == 0 :
            for theta in range (ntheta + 1):
                coordx[number] = dr * np.sin(np.pi / ntheta * theta) * np.cos(2 * np.pi / nphi *  phi )
                coordy[number] = dr * np.sin(np.pi / ntheta * theta) * np.sin(2 * np.pi / nphi *  phi )
                coordz[number] = dr * np.cos(np.pi / ntheta * theta) 
                number += 1

        else:
            for theta in range (1,ntheta):
                coordx[number] = dr * np.sin(np.pi / ntheta * theta) * np.cos(2 * np.pi / nphi *  phi )
                coordy[number] = dr * np.sin(np.pi / ntheta * theta) * np.sin(2 * np.pi / nphi *  phi )
                coordz[number] = dr * np.cos(np.pi / ntheta * theta) 
                number += 1
    coord = np.array([coordx,coordy,coordz])
    coord = coord.T
    return coord,number-1
   


if __name__ =="__main__":
    main()
