import numpy as np
import sys
from scipy.interpolate import RegularGridInterpolator
import call_lebedev

lebedev_num_list = (6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810)


def jump_index(V_ang,x,y,z,num,nr,rofi):
    my_inter_func = RegularGridInterpolator((x, y, z), V_ang)
    lebn = lebedev_num_list[num]
    lebx =np.zeros(lebn)
    leby =np.zeros(lebn)
    lebz =np.zeros(lebn)
    lebw =np.zeros(lebn)
    call_lebedev.lebedev(lebn,lebx,leby,lebz,lebw)
    Vr_omega = np.zeros((nr,lebn))
    Vomega_r_dr = np.zeros((lebn,nr))
    Vomega_r_dr2 = np.zeros((lebn,nr))
    for r in range(nr):
        Vr_omega[r] = my_inter_func(np.array([lebx,leby,lebz]).T * rofi[r])
    Vomega_r = Vr_omega.T

    for omega in range(lebn):
        Vomega_r_dr[omega] = np.gradient(Vomega_r[omega],rofi)
        Vomega_r_dr2[omega] = np.gradient(Vomega_r_dr[omega],rofi)

    for r in range(nr):
        print(Vomega_r[1730][r])
    print("dr")
    for r in range(nr):
        print(Vomega_r_dr[1730][r])
    for r in range(nr):
        print(Vomega_r_dr2[1730][r])
