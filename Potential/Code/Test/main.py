import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from surfaceintegral import surfaceintegral
from call_rseq import rseq
from makepotential import makepotential
from basis import make_basis
from integrate import integrate
from scipy import interpolate
from scipy.integrate import tplquad,trapz
from scipy.integrate import simps,cumtrapz
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator
from scipy.special import sph_harm
from call_genev import solve_genev
from make_V_radial_new import make_V_radial_new
from make_environment import make_environment

EPSVAL = 1.e-20
BETAZERO = 1.e-8

def main():
    # make environment
    fw_t = open("time_log",mode="w")
    t1 = time.time()
  
    pot_region,bound_rad,radius,region,nr,gridpx,gridpy,gridpz,x,y,z,xx,yy,zz,a,b,rofi,pot_type,pot_bottom,pot_show_f,si_method,radial_pot_show_f,new_radial_pot_show_f,node_open,node_close,LMAX = make_environment()

    # make potential
    V = makepotential(xx,yy,zz,pot_region,pot_type=pot_type,pot_bottom=pot_bottom,pot_show_f=pot_show_f)

    # surface integral
    V_radial = surfaceintegral(x,y,z,rofi,V,si_method=si_method,radial_pot_show_f=radial_pot_show_f)

    V_radial_new = make_V_radial_new(V_radial,rofi,pot_region,bound_rad,new_radial_pot_show_f=new_radial_pot_show_f)

    vofi = np.array (V_radial_new)  # select method of surface integral

    # make basis
    node = node_open + node_close
    nstates = node * LMAX**2

    all_basis = make_basis(LMAX,node_open,node_close,nr,a,b,rofi,vofi)
 
    #make spherical matrix element
    Smat = np.zeros((LMAX,LMAX,node,node),dtype = np.float64)
    hsmat = np.zeros((LMAX,LMAX,node,node),dtype = np.float64)
    lmat = np.zeros((LMAX,LMAX,node,node), dtype = np.float64)
    qmat = np.zeros((LMAX,LMAX,node,node), dtype = np.float64)

    for l1 in range (LMAX):
        for l2 in range (LMAX):
            if l1 != l2 :
                continue
            for n1 in range (node):
                for n2 in range (node):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                        sys.exit()
                    Smat[l1][l2][n1][n2] = integrate(all_basis[l1][n1].g[:nr] * all_basis[l2][n2].g[:nr],rofi,nr)
                    hsmat[l1][l2][n1][n2] = Smat[l1][l2][n1][n2] * all_basis[l2][n2].e
                    lmat[l1][l2][n1][n2] = all_basis[l1][n1].val * all_basis[l2][n2].slo
                    qmat[l1][l2][n1][n2] = all_basis[l1][n1].val * all_basis[l2][n2].val
    print("\nSmat")
    print(Smat)
    print ("\nhsmat")
    print (hsmat)
    print ("\nlmat")
    print (lmat)
    print ("\nqmat")
    print (qmat)

    hs_L = np.zeros((LMAX,LMAX,node,node))
    """
    for l1 in range(LMAX):
        for n1 in range(node):
            for l2 in range(LMAX):
                for n2 in range(node):
                    print("{:>8.4f}".format(lmat[l1][l2][n1][n2]),end="")
            print("")
    print("")
    """
    print("hs_L")
    for l1 in range(LMAX):
        for n1 in range(node):
            for l2 in range(LMAX):
                for n2 in range(node):
                    hs_L[l1][l2][n1][n2] = hsmat[l1][l2][n1][n2] + lmat[l1][l2][n1][n2]
                    print("{:>8.4f}".format(hs_L[l1][l2][n1][n2]),end="")
            print("")
    

    t2 = time.time()
    fw_t.write("make environment, make basis and calculate spherical matrix element\n")
    fw_t.write("time = {:>11.8}s\n".format(t2-t1))
    t1 = time.time()

    #make not spherical potential
    my_radial_interfunc = interpolate.interp1d(rofi, V_radial_new)

    V_ang = np.where(np.sqrt(xx * xx + yy * yy + zz * zz) < rofi[-1] , V - my_radial_interfunc(np.sqrt(xx **2 + yy **2 + zz **2)),0. )
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

    ngrid = 10
    fw_u = open("umat_gauss.dat",mode="w")

    umat = np.zeros((LMAX,2 * LMAX + 1,node,LMAX,2 * LMAX + 1,node), dtype = np.complex64)


    igridpx = ngrid
    igridpy = ngrid
    igridpz = ngrid

    gauss_px,gauss_wx = np.polynomial.legendre.leggauss(igridpx)
    gauss_py,gauss_wy = np.polynomial.legendre.leggauss(igridpy)
    gauss_pz,gauss_wz = np.polynomial.legendre.leggauss(igridpz)

    ix,iy,iz = gauss_px * pot_region[0],gauss_py * pot_region[1],gauss_pz * pot_region[2]
    ixx,iyy,izz = np.meshgrid(ix,iy,iz)

    gauss_weight = np.zeros((igridpx,igridpy,igridpz))
    for i in range(igridpx):
        for j in range(igridpy):
            for k in range(igridpz):
                gauss_weight[i][j][k] = gauss_wx[i] * gauss_wy[j] * gauss_wz[k]

    #my_V_ang_inter_func = RegularGridInterpolator((x, y, z), V_ang)
    #V_ang_i = my_V_ang_inter_func((ixx,iyy,izz))

    V_ang_i = np.zeros((igridpx,igridpy,igridpz))
    V_ang_i = -1. - my_radial_interfunc(np.sqrt(ixx * ixx + iyy * iyy + izz * izz))



    #mlab.contour3d(V_ang_i,color = (1,1,1),opacity = 0.1)
    #obj = mlab.volume_slice(V_ang_i)
    #mlab.show()

    dis = np.sqrt(ixx **2 + iyy **2 + izz **2)
    dis2 = ixx **2 + iyy **2 + izz **2
    theta = np.where( dis != 0., np.arccos(izz / dis), 0.)
    phi = np.where( iyy**2 + ixx **2 != 0 , np.where(iyy >= 0, np.arccos(ixx / np.sqrt(ixx **2 + iyy **2)), np.pi + np.arccos(ixx / np.sqrt(ixx **2 + iyy **2))), 0.)
    #phi = np.where( iyy != 0. , phi, 0.)
    #phi = np.where(iyy > 0, phi,-phi)
#    region_t = np.where(dis < rofi[-1],1,0)
    sph_harm_mat = np.zeros((LMAX,2 * LMAX + 1, igridpx,igridpy,igridpz),dtype = np.complex64)
    g_V_g = np.zeros((LMAX,node,LMAX,node,igridpx,igridpy,igridpz))
    for l1 in range (LMAX):
        for m1 in range (-l1,l1 + 1):
            sph_harm_mat[l1][m1] = np.where(dis != 0., sph_harm(m1,l1,phi,theta),0.)
    g_ln_mat = np.zeros((LMAX,node,igridpx,igridpy,igridpz),dtype = np.float64)
    for l1 in range (LMAX):
        for n1 in range (node):
            my_radial_g_inter_func = interpolate.interp1d(rofi,all_basis[l1][n1].g[:nr])
            g_ln_mat[l1][n1] = my_radial_g_inter_func(np.sqrt(ixx **2 + iyy **2 + izz **2))

    for l1 in range (LMAX):
        for n1 in range (node):
            for l2 in range (LMAX):
                for n2 in range (node):
                    # to avoid nan in region where it can not interpolate ie: dis > rofi
                    g_V_g[l1][n1][l2][n2] = np.where(dis < rofi[-1], g_ln_mat[l1][n1] * V_ang_i * g_ln_mat[l2][n2], 0.)

    for l1 in range (LMAX):
        for m1 in range (-l1,l1+1):
            for l2 in range (l1+1):
                for m2 in range (-l2,l2+1):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                        sys.exit()
                    for n1 in range (node):
                        for n2 in range (node):
                            #print("n1 = {} n2 = {} l1 = {} l2 = {} m1 = {} m2 = {}".format(n1,n2,l1,l2,m1,m2))
                            umat[l1][m1][n1][l2][m2][n2] = np.sum(np.where( dis2 != 0., sph_harm_mat[l1][m1].conjugate() * g_V_g[l1][n1][l2][n2] * sph_harm_mat[l2][m2] / dis2 * gauss_weight, 0.)) * pot_region[0] * pot_region[1] * pot_region[2] 
                            #print(umat[n1][n2][l1][l2][m1][m2])

    for l1 in range (LMAX):
        for m1 in range (-l1,l1+1):
            for l2 in range (l1+1,LMAX):
                for m2 in range (-l2,l2+1):
                    if all_basis[l1][n1].l != l1 or all_basis[l2][n2].l != l2:
                        print("error: L is differnt")
                        sys.exit()
                    for n1 in range (node):
                        for n2 in range (node):
                            #print("n1 = {} n2 = {} l1 = {} l2 = {} m1 = {} m2 = {}".format(n1,n2,l1,l2,m1,m2))
                            umat[l1][m1][n1][l2][m2][n2] = umat[l2][m2][n2][l1][m1][n1].conjugate()
                            #print(umat[n1][n2][l1][l2][m1][m2])
 


    count = 0
    for l1 in range (LMAX):
        for l2 in range (LMAX):
            for m1 in range (-l1,l1+1):
                for m2 in range (-l2,l2+1):
                    for n1 in range (node):
                        for n2 in range (node):
                            fw_u.write(str(count))
                            fw_u.write("{:>15.8f}{:>15.8f}\n".format(umat[l1][m1][n1][l2][m2][n2].real,umat[l1][m1][n1][l2][m2][n2].imag))
                            count += 1


    fw_u.close()

    t2 = time.time()
    fw_t.write("calculate U-matrix\n")
    fw_t.write("grid = {:<5} time ={:>11.8}s\n".format(ngrid,t2 - t1))

    
    ham_mat = np.zeros((nstates * nstates),dtype = np.float64)
    qmetric_mat = np.zeros((nstates * nstates),dtype = np.float64)
    for l1 in range(LMAX):
        for m1 in range (-l1,l1+1):
            for n1 in range(node):
                for l2 in range(LMAX):
                    for m2 in range(-l2,l2+1):
                        for n2 in range(node):
                            if l1 == l2 and m1 == m2 :
                                ham_mat[l1 **2 * node * LMAX**2 * node + (l1 + m1) * node * LMAX**2 * node + n1 * LMAX**2 * node + l2 **2 * node + (l2 + m2) * node + n2] += hs_L[l1][l2][n1][n2]
                                qmetric_mat[l1 **2 * node * LMAX**2 * node + (l1 + m1) * node * LMAX**2 * node + n1 * LMAX**2 * node + l2 **2 * node + (l2 + m2) * node + n2] = qmat[l1][l2][n1][n2]
                            ham_mat[l1 **2 * node * LMAX**2 * node + (l1 + m1) * node * LMAX**2 * node + n1 * LMAX**2 * node + l2 **2 * node + (l2 + m2) * node + n2] += umat[l1][m1][n1][l2][m2][n2].real

    lambda_mat = np.zeros(nstates * nstates,dtype = np.float64)
    alphalong = np.zeros(nstates)
    betalong = np.zeros(nstates)
    alpha = np.zeros(LMAX**2)
    beta = np.zeros(LMAX**2)
    reveclong = np.zeros(nstates * nstates)
    revec = np.zeros((LMAX**2,nstates))
    normfac = np.zeros(LMAX**2)
    Wlk = np.zeros((LMAX**2,LMAX**2))
    R_matrix = np.zeros((LMAX**2,LMAX**2))
    t1 = time.time()

    for e_num in range(1,2):
        E = e_num * 0
        lambda_mat = np.zeros(nstates * nstates,dtype = np.float64)
        lambda_mat -= ham_mat
        """
        for i in range(nstates):
            lambda_mat[i + i * nstates] += E
        """
        count = 0
        fw_lam = open("lambda.dat",mode="w")

        for l1 in range(LMAX):
            for m1 in range (-l1,l1+1):
                for n1 in range(node):
                    for l2 in range(LMAX):
                        for m2 in range(-l2,l2+1):
                            for n2 in range(node):
                                if l1 == l2 and m1 == m2 :
                                    lambda_mat[l1 **2 * node * LMAX**2 * node + (l1 + m1) * node * LMAX**2 * node + n1 * LMAX**2 * node + l2 **2 * node + (l2 + m2) * node + n2] += Smat[l1][l2][n1][n2] * E
                                fw_lam.write(str(count))
                                fw_lam.write("{:>15.8f}\n".format(lambda_mat[l1 **2 * node * LMAX**2 * node + (l1 + m1) * node * LMAX**2 * node + n1 * LMAX**2 * node + l2 **2 * node + (l2 + m2) * node + n2]))
                                count += 1
        fw_lam.close()
 
        
        info = solve_genev(nstates,lambda_mat,qmetric_mat,alphalong,betalong,reveclong)

        print("info = ",info)
        print("alphalong")
        print(alphalong)
        print("betalong")
        print(betalong)
        #print(revec)
        
        k = 0
        jk = np.zeros(nstates,dtype=np.int32) 
        for i in range(nstates):
            if abs(betalong[i]) > BETAZERO :
                jk[k] = i
                k += 1

        if k != LMAX**2:
            print(" main PANIC: solution of genev has rank k != nchannels ")
            sys.exit()

        for k in range(LMAX**2):
            j = jk[k]
            alpha[k] = alphalong[j]
            beta[k]  = betalong[j]
            revec[k] = reveclong[j * nstates : (j + 1) * nstates]

        for k in range(LMAX**2):
            dum = 0
            for l in range(LMAX):
                for m in range(-l,l+1):
                    for n in range(node):
                        Wlk[l**2 + (l+m)][k] += revec[k][l**2*node + (l+m)*node + n] * all_basis[l][n].g[nr-1]
                    dum += Wlk[l**2 + (l+m)][k] **2
            normfac[k] = 1 / np.sqrt(dum)
            Wlk[:,k] *= normfac[k]
            """
            for l in range(LMAX):
                for m in range(-l,l+1):
                    Wlk[l**2 + (l+m)][k] *= normfac[k]
            """

        for l in range(LMAX):
            for m in range(-l,l+1):
                for k in range(LMAX**2):
                    print("{:>8.4f}".format(Wlk[l**2 + (l + m)][k]),end=" ")
                print("")

        print("R-matrix")
        for l1 in range(LMAX):
            for m1 in range(-l1,l1+1):
                for l2 in range(LMAX):
                    for m2 in range(-l2,l2+1):
                        for k in range(LMAX**2):
                            R_matrix[l1**2+(l1+m1)][l2**2+(l2+m2)] -= Wlk[l1**2 + (l1+m1)][k] * beta[k] / alpha[k] * Wlk[l2**2 + (l2+m2)][k]
                        print ("{:>8.4f}".format(R_matrix[l1**2+(l1+m1)][l2**2+(l2+m2)]),end="")
                print("")

        sys.exit()



        fw_revec = open("revec.dat",mode="w")
        print("revec")
        for j in range(nstates):
            fw_revec.write("{:>4}".format(j))
            for i in range(LMAX**2):
                fw_revec.write("{:>13.8f}".format(revec[i][j]))
                print("{:>11.6f}".format(revec[i][j]),end="")
            print("")
            fw_revec.write("\n")
        print("")
        fw_revec.close()

        t2 = time.time()
        fw_t.write("solve eigenvalue progrem\n")
        fw_t.write("time = {:>11.8}s\n".format(t2-t1))
        fw_t.close()






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
