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
from scipy.integrate import tplquad
from scipy.integrate import simps,cumtrapz
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator
from scipy.special import sph_harm

EPSVAL = 1.e-20

def main():
    # make environment
    t1 = time.time()
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

    #mlab.points3d(np.where(np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] , - my_radial_interfunc(np.sqrt(xx * xx + yy * yy + zz * zz)),0. ),scale_factor=0.4)
    #mlab.show()
    #mlab.points3d(V,scale_factor=0.4)
    #mlab.show()
    #sys.exit()
    V_ang = np.where(np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] , V - my_radial_interfunc(np.sqrt(xx * xx + yy * yy + zz * zz)),0. )
    """
    for i in range(gridpx):
        for j in range(gridpy):
            for k in range(gridpz):
                print(V_ang[i][j][k],end="  ")
            print("")
        print("\n")
    sys.exit()
    """

    #WARING!!!!!!!!!!!!!!!!!!!!!!
    """
    Fujikata rewrote ~/.local/lib/python3.6/site-packages/scipy/interpolate/interpolate.py line 690~702
    To avoid exit with error "A value in x_new is below the interpolation range."
    """
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!

    """
    with open ("V_ang.dat",mode = "w") as fw_a:
        for i in range(len(V_ang)):
            fw_a.write("{:>13.8f}".format(xx[50][i][50]))
            fw_a.write("{:>13.8f}\n".format(V_ang[i][50][50]))
    """

    #mayavi.mlab.plot3d(xx,yy,V_ang)
#    mlab.contour3d(V_ang,color = (1,1,1),opacity = 0.1)
#    obj = mlab.volume_slice(V_ang)
#    mlab.show()

    umat = np.zeros((node_open + node_close,node_open + node_close,LMAX,LMAX,2 * LMAX + 1,2 * LMAX + 1), dtype = np.complex64)
    #umat_t = np.zeros((LMAX,LMAX,2 * LMAX + 1,2 * LMAX + 1,node_open + node_close,node_open + node_close), dtype = np.complex64)

    my_V_ang_inter_func = RegularGridInterpolator((x, y, z), V_ang)

    igridpx = 51
    igridpy = 51
    igridpz = 51

    ix,iy,iz = grid(igridpx,igridpy,igridpz,region)
    ixx,iyy,izz = np.meshgrid(ix,iy,iz)

    V_ang_i = my_V_ang_inter_func((ixx,iyy,izz))

    mlab.points3d(V_ang_i,scale_factor=0.4)
    mlab.show()
#    mlab.contour3d(V_ang_i,color = (1,1,1),opacity = 0.1)
#    obj = mlab.volume_slice(V_ang_i)
#    mlab.show()

    dis = np.sqrt(ixx **2 + iyy **2 + izz **2)
    dis2 = ixx **2 + iyy **2 + izz **2
    theta = np.where( dis != 0., np.arccos(izz / dis), 0.)
    phi = np.where( ixx**2 + iyy **2 != 0. , np.arccos(ixx / np.sqrt(ixx **2 + iyy **2)), 0.)
#    region_t = np.where(dis < rofi[-1],1,0)
    sph_harm_mat = np.zeros((LMAX,2 * LMAX + 1, igridpx,igridpy,igridpz),dtype = np.complex64)
    for l1 in range (LMAX):
        for m1 in range (-l1,l1 + 1):
            sph_harm_mat[l1][m1] = np.where(dis != 0., sph_harm(m1,l1,theta,phi),0.)
            mlab.points3d(sph_harm_mat[l1][m1].real,scale_factor=0.4)
            mlab.show()
    g_ln_mat = np.zeros((node_open + node_close,LMAX,igridpx,igridpy,igridpz),dtype = np.float64)
    for n1 in range (node_open + node_close):
        for l1 in range (LMAX):
            my_radial_g_inter_func = interpolate.interp1d(rofi,all_basis[l1][n1].g[:nr])
            g_ln_mat[n1][l1] = my_radial_g_inter_func(np.sqrt(ixx **2 + iyy **2 + izz **2))
            #mlab.points3d(g_ln_mat[n1][l1],scale_factor=0.4)
            #mlab.show()

    for n1 in range (node_open + node_close):
        for n2 in range (node_open + node_close):
            for l1 in range (LMAX):
                for l2 in range (LMAX):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                        sys.exit()
                    # to avoid nan in region where it can not interpolate ie: dis > rofi
                    g_V_g = np.where(dis < rofi[-1], g_ln_mat[n1][l1] * V_ang_i * g_ln_mat[n2][l2], 0.)
                    #mlab.points3d(g_V_g,scale_factor=0.4)
                    #mlab.show()
                    for m1 in range (-l1,l1+1):
                        for m2 in range (-l2,l2+1):
                            print("n1 = {} n2 = {} l1 = {} l2 = {} m1 = {} m2 = {}".format(n1,n2,l1,l2,m1,m2))
                            #umat[n1][n2][l1][l2][m1][m2] = np.sum(np.where(dis2 != 0., sph_harm_mat[l1][m1] * g_V_g * sph_harm_mat[l2][m2].conjugate() / dis2, 0.)) * (2 * region[0] * 2 * region[1] * 2 * region[2]) / (igridpx * igridpy * igridpz)
                            mlab.points3d(np.where(dis2 != 0., sph_harm_mat[l1][m1] * g_V_g * sph_harm_mat[l2][m2].conjugate() / dis2,0.).real,scale_factor=6.)
                            mlab.show()
                            #umat[l1][l2][m1][m2][n1][n2] = simps(simps(simps(g_V_g * sph_harm_mat[l1][m1] * sph_harm_mat[l2][m2],x = z,even="first"),x = y,even="first"),x = x,even="first")
                            #umat[l1][l2][m1][m2][n1][n2] = simps(simps(simps(g_V_g * sph_harm_mat[l1][m1] * sph_harm_mat[l2][m2],x = z),x = y),x = x)
                            #umat[l1][l2][m1][m2][n1][n2] = cumtrapz(cumtrapz(cumtrapz(g_V_g * sph_harm_mat[l1][m1] * sph_harm_mat[l2][m2],x = z),x = y),x = x)
                            #umat[l1][l2][m1][m2][n1][n2] = np.sum( np.where( dis < rofi[-1] ,sph_harm(m1,l1,theta,phi) * g1 * V_ang * sph_harm(m2,l2,theta,phi) * g2 ,0. )) / (gridpx * gridpy * gridpz)
#                            umat_t[l1][l2][m1][m2][n1][n2] = np.sum( np.where( np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] ,sph_harm(m1,l1,np.arccos(zz / np.sqrt(xx **2 + yy **2 + zz **2)),np.arccos(xx / np.sqrt(xx **2 + yy **2))) * my_radial_g1_inter_func(np.sqrt(xx**2 + yy **2 + zz **2)) * my_V_ang_inter_func((xx,yy,zz)) * sph_harm(m2,l2,np.arccos(zz / np.sqrt(xx **2 + yy **2 + zz **2)),np.arccos(xx / np.sqrt(xx **2 + yy **2))) * my_radial_g2_inter_func(np.sqrt(xx **2 + yy **2 + zz **2)),0. ))
#                            print(umat_t[l1][l2][m1][m2][n1][n2])
                            """
                            for xi in x:
                                for yi in y:
                                    for zi in z:
                                        if np.sqrt(xi**2 + yi **2 + zi **2) < rofi[-1]:
                                            umat_t[l1][l2][m1][m2][n1][n2] += sph_harm(m1,l1,np.arccos(zi / np.sqrt(xi **2 + yi **2 + zi **2)),np.arccos(xi / np.sqrt(xi **2 + yi **2))) * my_radial_g1_inter_func(np.sqrt(xi**2 + yi **2 + zi **2)) * my_V_ang_inter_func((xi,yi,zi)) * sph_harm(m2,l2,np.arccos(zi / np.sqrt(xi **2 + yi **2 + zi **2)),np.arccos(xi / np.sqrt(xi **2 + yi **2))) * my_radial_g2_inter_func(np.sqrt(xi **2 + yi **2 + zi **2))
                            print(umat[l1][l2][m1][m2][n1][n2])
                            print(umat_t[l1][l2][m1][m2][n1][n2])
                            print(umat_t[l1][l2][m1][m2][n1][n2] - umat[l1][l2][m1][m2][n1][n2])
                            """

                            #umat[l1][l2][m1][m2][n1][n2] = tplquad( lambda x, y, z : sph_harm(m1,l1,np.arccos(z / np.sqrt(x **2 + y **2 + z **2)),np.arccos(x / np.sqrt(x **2 + y **2))).real * my_radial_g1_inter_func(np.sqrt(x**2 + y **2 + z **2)) * my_V_ang_inter_func((x,y,z)) * sph_harm(m2,l2,np.arccos(z / np.sqrt(x **2 + y **2 + z **2)),np.arccos(x / np.sqrt(x **2 + y **2))).real * my_radial_g2_inter_func(np.sqrt(x **2 + y **2 + z **2)) if np.sqrt(x **2 + y **2 + z **2) < rofi[-1] else 0 , -region, region , -region, region, -region, region)
                            print(umat[n1][n2][l1][l2][m1][m2])
                                

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
