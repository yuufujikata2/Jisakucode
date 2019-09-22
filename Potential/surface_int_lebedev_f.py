import numpy as np
import call_lebedev

lebedev_num_list = (6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810)

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


if __name__ =="__main__":
    main()
