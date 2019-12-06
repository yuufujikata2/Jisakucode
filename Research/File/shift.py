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

    print("Please input shift volume. (All volumes for each culmns is necessary. noshift is 0)")
    input_shift = input (">>")
    shift_list = input_shift.split()
    shift_volume = [float(i) for i in shift_list]

    all_data = []
    for i in range (1,argc):
        data = Data(args[i],skiprows,maxrows,column)
        all_data.append(data)

    fw = open("newdata.dat",mode="w")
    for i in range(len(all_data[0].x)):
        fw.write("{:>13.6f}".format(all_data[0].x[i] + shift_volume[0]))
        for j in range(len(all_data[0].y)):
            fw.write("{:>13.6f}".format(all_data[0].y[j][i] + shift_volume[j+1]))
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
                
            




    
if __name__ == "__main__":
    main()

        
