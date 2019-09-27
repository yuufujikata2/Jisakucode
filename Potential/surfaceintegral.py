import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import rasen
import call_lebedev
import matplotlib.pyplot as plt

lebedev_num_list = (6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810)

def surfaceintegral(x,y,z,r,V,method="lebedev_f",potshow_f=False):
    fw = open("V_radial.dat",mode="w")
    my_inter_func = RegularGridInterpolator((x, y, z), V)
    V_radial = []
    if method == "lebedev_f" :
        surface_int_method = SI_lebedev_f(-1)
    elif method == "lebedev_py" :
        surface_int_method = SI_lebedev_py()
    elif method == "rasen" :
        surface_int_method = SI_rasen(5000)
    else :
        print ("error: There is no method entered")
        sys.exit()

    for dr in r :
        integral_result = surface_int_method.integral(my_inter_func, dr)
        V_radial.append(integral_result)
        fw.write("{:>11.6f}".format(dr))
        fw.write("{:>11.6f}".format(integral_result))
        fw.write("\n")

    fw.close()

    if potshow_f :
        potentialshow(r,V_radial)

    return V_radial

class SI_lebedev_f():
    def __init__(self,num):
        self.lebedev_num = lebedev_num_list[num]
        self.lebedev_x = np.zeros(self.lebedev_num)
        self.lebedev_y = np.zeros(self.lebedev_num)
        self.lebedev_z = np.zeros(self.lebedev_num)
        self.lebedev_w = np.zeros(self.lebedev_num)
        call_lebedev.lebedev(self.lebedev_num,self.lebedev_x,self.lebedev_y,self.lebedev_z,self.lebedev_w)

    def integral(self,my_inter_func,dr):
        return  sum(my_inter_func(np.array([self.lebedev_x,self.lebedev_y,self.lebedev_z]).T * dr) * self.lebedev_w)

class SI_lebedev_py():
    def __init__ (self):
        #self.scheme = quadpy.sphere.lebedev_019()
        self.scheme = quadpy.sphere.lebedev_131()

    def integral(self,my_inter_func,dr):
        if dr == 0. :
            return my_inter_func([0.0,0.0,0.0])[0]
        else:
            return self.scheme.integrate( my_inter_func,[0.0,0.0,0.0],dr) / (4 * np.pi * dr**2)

class SI_rasen():
    def __init__(self, npoints):
        self.npoints = npoints
        self.gss = rasen.generate_gss(npoints)
        self.relocated_gss = rasen.relocate(self.gss)
        self.coord = np.array([rasen.spherical2cartesian (c) for c in self.relocated_gss])

    def integral(self,my_inter_func,dr):
        newV = my_inter_func(self.coord * dr)
        V_average = np.average(newV)
        return V_average

def potentialshow(r,V):
    plt.plot(r,V,marker=".")
    plt.show()
