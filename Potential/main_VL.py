import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from surfaceintegral import surfaceintegral
from call_rseq import rseq
from call_lebedev import lebedev
from makepotential import makepotential
from basis import Basis
from integrate import integrate
from scipy import interpolate
from scipy.integrate import tplquad
from scipy.integrate import simps,cumtrapz
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator
from scipy.special import sph_harm
from sympy.physics.quantum.cg import Wigner3j

EPSVAL = 1.e-20
lebedev_num_list = (6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810)

def main():
    # make environment
    pot_region = (2/np.sqrt(3),2/np.sqrt(3),2/np.sqrt(3))
    radius = np.sqrt(pot_region[0] **2 + pot_region[1] **2 + pot_region[2] **2 )
    region = (radius,radius,radius)
    nr = 201
    gridpx = 200
    gridpy = 200 
    gridpz = 200 
    x,y,z = grid(gridpx,gridpy,gridpz,region)
    xx, yy, zz = np.meshgrid(x,y,z)

    #make mesh

    #linear mesh
    #r =np.linspace(0.0001,region,nr)

    #log mesh
    a = np.log(2) / (nr - 1) 
    b = radius / (np.e**(a * ( nr - 1)) - 1)
    rofi = np.array([b * (np.e**(a * i) - 1) for i in range(nr)])

    # make potential
    V = makepotential(xx,yy,zz,pot_region,pottype="cubic",potbottom=-1,potshow_f=False)

    # surface integral
    V_radial = surfaceintegral(x,y,z,rofi,V,method="lebedev_py",potshow_f=False)
    vofi = np.array (V_radial)  # select method of surface integral

    # make basis
    
    node_open = 1
    node_close = 2
    LMAX = 4

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
        fw_w.write("#r,   l, node, open or close = ") 
        for l_basis in all_basis:
            for nyu_basis in l_basis:
                fw_w.write(str(nyu_basis.l))
                fw_w.write(str(nyu_basis.node))
                if nyu_basis.open:
                    fw_w.write("open")
                else:
                    fw_w.write("close")
                fw_w.write("    ")
        fw_w.write("\n")

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
                        sys.exit()
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
    my_radial_interfunc = interpolate.interp1d(rofi, V_radial)

    #V_ang = np.where(np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] , V - my_radial_interfunc(np.sqrt(xx * xx + yy * yy + zz * zz)),0. )
    #my_V_ang_inter_func = RegularGridInterpolator((x, y, z), V_ang)
    my_V_inter_func = RegularGridInterpolator((x, y, z), V)

    #WARING!!!!!!!!!!!!!!!!!!!!!!
    """
    Fujikata rewrote ~/.local/lib/python3.6/site-packages/scipy/interpolate/interpolate.py line 690~702
    To avoid exit with error "A value in x_new is below the interpolation range."
    """
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!


    #mayavi.mlab.plot3d(xx,yy,V_ang)
#    mlab.contour3d(V_ang,color = (1,1,1),opacity = 0.1)
#    obj = mlab.volume_slice(V_ang)
#    mlab.show()


    #for V_L
    fw_umat_vl = open("umat_vl.dat",mode="w")
    for LMAX_k in range(3,10):
        t1 = time.time()
        umat = np.zeros((node_open + node_close,node_open + node_close,LMAX,LMAX,2 * LMAX + 1,2 * LMAX + 1), dtype = np.complex64)
        #LMAX_k = 7
        #fw_umat_vl.write(str(lebnum))
        fw_umat_vl.write(str(LMAX_k))
        igridnr = 200
        leb_r = np.linspace(0,radius,igridnr)
        lebedev_num = lebedev_num_list[-8]
    
        V_L = np.zeros((LMAX_k, 2 * LMAX_k + 1, igridnr), dtype = np.complex64)
        leb_x =np.zeros(lebedev_num)
        leb_y =np.zeros(lebedev_num)
        leb_z =np.zeros(lebedev_num)
        leb_w =np.zeros(lebedev_num)
    
        lebedev(lebedev_num,leb_x,leb_y,leb_z,leb_w)
    
        theta = np.arccos(leb_z)
        phi = np.where( leb_x **2 + leb_y **2 != 0. , np.arccos(leb_x / np.sqrt(leb_x **2 + leb_y **2)), 0.)
    
        for l1 in range(LMAX_k):
            for m1 in range(-l1,l1+1):
                for i in range(igridnr):
                    V_L[l1][m1][i] = np.sum(my_V_inter_func((leb_r[i] * leb_x, leb_r[i] * leb_y, leb_r[i] * leb_z)) * leb_w * sph_harm(m1,l1,theta,phi).conjugate())
                #print("l = ",l1,"m = ", m1)
                #plt.plot(leb_r,V_L[l1][m1].real,marker=".")
                #plt.show()
    
        g_ln = np.zeros((node_open + node_close,LMAX,igridnr),dtype = np.float64)
        for n1 in range (node_open + node_close):
            for l1 in range (LMAX):
                my_radial_g_inter_func = interpolate.interp1d(rofi,all_basis[l1][n1].g[:nr])
                g_ln[n1][l1] = my_radial_g_inter_func(leb_r)
        for n1 in range (node_open + node_close):
            for n2 in range (node_open + node_close):
                for l1 in range (LMAX):
                    for l2 in range (LMAX):
                        for m1 in range(-l1,l1+1):
                            for m2 in range(-l2,l2+1):
                                for k in range(1,LMAX_k):
                                    for q in range(-k,k+1):
                                        umat[n1][n2][l1][l2][m1][m2] += simps(g_ln[n1][l1] * V_L[k][q] * g_ln[n2][l2], leb_r) * (-1) **(-m1) * np.sqrt((2 * l1 + 1) * (2 * l2 +1)) *Wigner3j(l1,0,k,0,l2,0).doit() * Wigner3j(l1,-m1,k,q,l2,m2).doit() * np.sqrt((2 * k + 1) / (4 * np.pi))
                                #print(umat[n1][n2][l1][l2][m1][m2])
                                fw_umat_vl.write("{:>15.8f}".format(umat[n1][n2][l1][l2][m1][m2].real))
    
        t2 = time.time()
        #print("number of lebedev grid = ",lebedev_num_list[lebnum]," time = ",t2 - t1)
        print("LMAX = ",LMAX_k," time = ",t2 - t1)
    
        fw_umat_vl.write("\n")


    sys.exit()
  



    fw_u = open("umat_nlm.dat",mode="w")
    count = 0

    for n1 in range (node_open + node_close):
        for n2 in range (node_open + node_close):
            for l1 in range (LMAX):
                for l2 in range (LMAX):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                        sys.exit()
                    # to avoid nan in region where it can not interpolate ie: dis > rofi
                    g_V_g = np.where(dis < rofi[-1], g_ln_mat[n1][l1] * V_ang_i * g_ln_mat[n2][l2], 0.)
                    for m1 in range (-l1,l1+1):
                        for m2 in range (-l2,l2+1):
                            print("n1 = {} n2 = {} l1 = {} l2 = {} m1 = {} m2 = {}".format(n1,n2,l1,l2,m1,m2))
                            umat[n1][n2][l1][l2][m1][m2] = np.sum(np.where(dis2 != 0., sph_harm_mat[l1][m1] * g_V_g * sph_harm_mat[l2][m2].conjugate() / dis2, 0.)) * (2 * region[0] * 2 * region[1] * 2 * region[2]) / (igridpx * igridpy * igridpz)
                            fw_u.write(str(count))
                            fw_u.write("{:>15.8f}{:>15.8f}\n".format(umat[n1][n2][l1][l2][m1][m2].real,umat[n1][n2][l1][l2][m1][m2].imag))
                            count += 1

                            print(umat[n1][n2][l1][l2][m1][m2])
                                

    fw_u.close()
    t2 = time.time()
    print("time = ",t2 - t1)


def grid (nx,ny,nz,region):
    x = np.linspace(-region[0],region[0],nx)
    y = np.linspace(-region[1],region[1],ny)
    z = np.linspace(-region[2],region[2],nz)
    return x,y,z


def fitting(r,V,d):
    w = np.polyfit(r,V,d)
    print(w)
    xs = np.linspace(r[0],r[-1],300)
    ys = np.polyval(w,xs)
    return xs,ys
    
if __name__ =="__main__":
    main()
