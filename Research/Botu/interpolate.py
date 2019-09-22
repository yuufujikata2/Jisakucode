import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator

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

    all_data = []
    for i in range (1,argc):
        data = Data(args[i],skiprows,maxrows,column)
        print("i = ",i,"x = ",data.x)
        print("i = ", i, "y = ",data.y)
        all_data.append(data)

    for i in range(len(all_data)) :
        all_data[i].newx = all_data[0].x
        all_data[i].interpolate()



class Data():
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
        self.newx = None
        self.newy = []

    def interpolate(self):
        for i in range(len(self.y)):
            my_inter = RegularGridInterpolator((self.x,), self.y[i])
            new_oney = my_inter(self.newx)
            self.newy.append(new_oney)
        print (self.newy)





    
if __name__ == "__main__":
    main()

        
