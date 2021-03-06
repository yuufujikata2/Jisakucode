import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import copy

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
    print("Please input using culmns (ex:1 2 3 or \"/\")")
    input_column = input(">>")
    if input_column == "/":
        column = None
    else:
        column = input_column.split()
        for i in range(len(column)) :
            column[i] = int(column[i]) -1

    print("sum of files a culmns or sum of columns in a file? (T or F)")
    input_columnsum = input(">>")
    if input_columnsum == "T":
        columnsum = True
    else:
        columnsum = False

    print("Please input ratio of sum data (dfalut = 1:1:1...)")
    input_sumratio = input (">>")
    if input_sumratio == "/":
        if columnsum :
            sumratio = [1 for i in  range(argc -1)]
        else:
            print("please input number of y columns")
            numofcolumn = int(input(">>"))
            sumratio = [1 for i in  range(numofcolumn)]
    else :
        sumratio_list = input_sumratio.split()
        sumratio = [float(i) for i in sumratio_list]

    all_data = []
    for i in range (1,argc):
        data = Data(args[i],skiprows,maxrows,column)
        all_data.append(data)

    Data.newx = copy.deepcopy(all_data[0].x)
    for i in range(len(all_data)):
        all_data[i].hairetusoroe()
    print("x_min = ",Data.newx[0])
    print("x_max = ",Data.newx[-1])

    for i in range(len(all_data)) :
        all_data[i].interpolate()

    fw = open("newdata.dat",mode="w")
    fw.write("#")
    fw.write(" y(file,data)\n")
    fw.write("#")
    fw.write("      x       ")
    if columnsum:
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
                    data_sum += all_data[k].newy[j][i] * sumratio[k]
                fw.write("{:>13.6f}".format(data_sum))
            fw.write("\n")
    else:
        fw.write("\n")
        for i in range(len(Data.newx)):
            fw.write("{:>13.6f}".format(Data.newx[i]))
            data_sum = 0
            for j in range(len(all_data[0].y)):
                #fw.write("{:>13.6f}".format(all_data[0].newy[j][i]))
                data_sum += all_data[0].newy[j][i] * sumratio[j]
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

    def hairetusoroe(self):
        count1 = 0
        count2 = 0
        while(True):
            if Data.newx[count1] >= self.x[0]:
                break
            count1 += 1
        while(True):
            if Data.newx[count2-1] <= self.x[-1]:
                break
            count2 -= 1
        if count2 == 0:
            Data.newx = Data.newx[count1:]
        else:
            Data.newx = Data.newx[count1:count2]
            
                
            

    def interpolate(self):
        for i in range(len(self.y)):
            #my_inter = RegularGridInterpolator((self.x,), self.y[i])
            #print ("self.newx = ",self.newx)
            #new_oney = my_inter(self.newx)
            f = interpolate.interp1d(self.x,self.y[i])
            new_oney = f(Data.newx)
            self.newy.append(new_oney)





    
if __name__ == "__main__":
    main()

        
