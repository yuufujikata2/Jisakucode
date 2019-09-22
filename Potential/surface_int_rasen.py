import numpy as np
import rasen

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

if __name__ =="__main__":
    main()
