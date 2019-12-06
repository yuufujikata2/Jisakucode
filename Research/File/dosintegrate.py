import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import copy
from scipy.integrate import simps

def main():
    args = sys.argv
    argc = len(args)

    print("Please input start row (default = 1 (\"/\"))")
    input_startrow = input(">>")
    if input_startrow == "/":
        skiprows = None
    else:
        skiprows = int(input_startrow)-1
    print("Please input last row (default = all (\"/\"))")
    input_lastrow = input(">>")
    if input_lastrow == "/":
        maxrows = None
    else:
        maxrows = int (input_lastrow) - skiprows
    print("Please input using culmns. First one is x and other are y (ex: 1 2 3   \"/\")")
    input_column = input(">>")
    if input_column == "/":
        column = None
    else:
        column = input_column.split()
        for i in range(len(column)) :
            column[i] = int(column[i]) -1

    print("Please input xmin and xmax")
    input_xrange = input(">>").split()
    xmin = float(input_xrange[0])
    xmax = float(input_xrange[1])

    """
    print("Please input ratio of sum data (dfalut = 1:1:1...)")
    input_sumratio = input (">>").split()
    if input_sumratio == "/":
        sumratio = [1 for i in  range(argc -1)]
    else :
        sumratio = [int(i) for i in input_sumratio]
    """

    all_data = []
    for i in range (1,argc):
        data = Data(args[i],skiprows,maxrows,column)
        all_data.append(data)

    Data.newx = copy.deepcopy(all_data[0].x)
    for i in range(len(all_data)):
        count1,count2 = all_data[i].integraterange(xmin,xmax)

    for i in range(len(all_data)) :
        all_data[i].interpolate()
    
    for i in range(len(all_data)):
        print(all_data[i].integrate())
        


    fw = open("newdata.dat",mode="w")
    fw.write("#")
    fw.write(" y(file,data)\n")
    fw.write("#")
    fw.write("      x       ")
    for i in range(len(all_data[0].y)):
        for j in range(len(all_data)):
            fw.write("    y({0},{1})   ".format(j+1,i+1))
        fw.write("  sum_data{}  ".format(i+1))
    fw.write("\n")
    for i in range(len(Data.newx)):
        fw.write("{:>13.6f}".format(Data.newx[i]))
        for j in range(len(all_data[0].y)):
            data_sum = 0
            for k in range(len(all_data)):
                fw.write("{:>13.6f}".format(all_data[k].newy[j][i]))
                data_sum += all_data[k].newy[j][i] #* sumratio[k]
            fw.write("{:>13.6f}".format(data_sum))
        fw.write("\n")
    fw.close()

class Data():
    newx = None
    def __init__(self,arg,skiprows,maxrows,column):
        if skiprows == None and maxrows == None and column == None:
            self.file = np.loadtxt(arg,unpack=True)
        elif skiprows == None and maxrows == None:
            self.file = np.loadtxt(arg,usecols=column,unpack=True)
        elif maxrows ==None and column == None:
            self.file = np.loadtxt(arg,skiprows=skiprows,unpack=True)
        elif column == None:
            self.file = np.loadtxt(arg,skiprows=skiprows,max_rows=maxrows,unpack=True)
        else:
            self.file = np.loadtxt(arg,skiprows=skiprows,usecols=column,max_rows=maxrows,unpack=True)
        self.x = self.file[0]
        self.y = self.file[1:]
        self.newy = []

    def integraterange(self,xmin,xmax):
        count1 = 0
        count2 = 0
        while(True):
            if Data.newx[count1] >= xmin:
                break
            count1 += 1
        while(True):
            if Data.newx[count2-1] <= xmax:
                break
            count2 -= 1
        if count2 == 0:
            Data.newx = Data.newx[count1:]
        else:
            Data.newx = Data.newx[count1:count2]
        
        return count1,count2

    def interpolate(self):
        for i in range(len(self.y)):
            #my_inter = RegularGridInterpolator((self.x,), self.y[i])
            #print ("self.newx = ",self.newx)
            #new_oney = my_inter(self.newx)
            f = interpolate.interp1d(self.x,self.y[i])
            new_oney = f(Data.newx)
            self.newy.append(new_oney)

    def integrate(self):
        integ_result=[]
        for i in range(len(self.y)):
            integ_resulti = simps(self.newy[i],self.newx)
            integ_result.append(integ_resulti)
        return integ_result






if __name__ == "__main__":
    main()

