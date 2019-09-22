import numpy as np
import quadpy

class SI_lebedev_py():
    def __init__ (self):
        #self.scheme = quadpy.sphere.lebedev_019()
        self.scheme = quadpy.sphere.lebedev_131()

    def integral(self,my_inter_func,dr):
        if dr == 0. :
            return my_inter_func([0.0,0.0,0.0])[0]
        else:
            return self.scheme.integrate( my_inter_func,[0.0,0.0,0.0],dr) / (4 * np.pi * dr**2)



if __name__ =="__main__":
    main()
