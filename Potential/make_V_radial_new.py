import numpy as np
import matplotlib.pyplot as plt

def make_V_radial_new(V_radial,rofi,pot_region,bound_rad,pot_show_f=False):
    r1 = pot_region[0] 
    r2 = bound_rad[0]
    r1_i = 0
    r2_i = 0
    while(True):
        if rofi[r2_i] < bound_rad[0] :
            r2_i += 1
        if rofi[r1_i] > pot_region[0] :
            break
        r1_i += 1
    a = - (V_radial[r1_i-2] + (r1 - rofi[r1_i - 2]) / (rofi[r1_i-1] -rofi[r1_i - 2]) * (V_radial[r1_i-1] - V_radial[r1_i -2]))
    #a = - V_radial[r1_i-1] 
    b = - (V_radial[r1_i-1] - V_radial[r1_i - 2]) / (rofi[r1_i -1] - rofi[r1_i-2 ]) 
    a = 0.98 
    b = -0.5

    V_radial_new = np.zeros_like(V_radial)

    for i in range(len(rofi)):
        r = rofi[i]
        if r < r2 :
            V_radial_new[i] = V_radial[i]
        elif r < r1 :
            f = (3 * a - b * (r1 - r2)) * ((r - r2) / (r1 - r2))**2 + (-2 * a + b * (r1 -r2)) * ((r - r2) / (r1 - r2))**3
            V_radial_new[i] = V_radial[i] + f
        else:
            V_radial_new[i] = 0.

    if pot_show_f:
        plt.plot(rofi,V_radial_new,marker=".")
        plt.plot(rofi,V_radial,marker=".")
        plt.show()

    return V_radial_new
