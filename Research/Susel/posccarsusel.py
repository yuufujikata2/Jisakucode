import sys
import numpy as np

def main():
    args = sys.argv
    argc = len(args)
    if argc != 2:
        print("error: please input a POSCAR!")
        sys.exit()
    with open(args[1],mode="r") as fp:
        firstline = fp.readline()
        secondline = fp.readline()
        plat = np.zeros((3,3))
        for i in range(3):
            line = fp.readline().split()
            for j in range(3):
                plat[i][j] = float(line[j])
        atomline = fp.readline()
        atom = atomline.split()
        line = fp.readline().split()
        atom_number = []
        for i in range(len(atom)):
            atom_number.append(int(line[i])) 
        sum_atom = sum(atom_number)
        coord_type = fp.readline().strip()
        coord = np.zeros((sum_atom,3))
        for i in range(sum_atom):
            line = fp.readline().split()
            for j in range(3):
                coord[i][j] = float(line[j])

    print(plat)
    print(atom)
    print(atom_number)
    print(coord)

    print ("\n")
    print("please input supercell size")
    su = input().split()
    if len(su) != 3 :
        print("error")
        sys.exit()
    sucel_size = []
    for i in range (3):
        sucel_size.append(int(su[i]))
    times_of_sucel =1
    for i in range(3):
        times_of_sucel = times_of_sucel * sucel_size[i]

    with open ("POSCAR_SU",mode="w") as fw:
        fw.write(firstline)
        fw.write(secondline)
        for i in range(3):
            plat[i] = plat[i]* sucel_size[i]
            for j in range(3):
                fw.write("{:>11.6f}".format(plat[i][j]))
            fw.write("\n")
        fw.write(atomline)
        for i in atom_number:
            fw.write("{:>4}".format(i*times_of_sucel))
        fw.write("\n")
        fw.write(coord_type,"\n")
        for i in range (len(coord)):
            for j in range(sucel_size[0]):
                for k in range(sucel_size[1]):
                    for l in range(sucel_size[2]):
                        fw.write("{:>11.6f}".format(coord[i][0]/sucel_size[0] + j/sucel_size[0]))
                        fw.write("{:>11.6f}".format(coord[i][1]/sucel_size[1] + k/sucel_size[1]))
                        fw.write("{:>11.6f}".format(coord[i][2]/sucel_size[2] + l/sucel_size[2]))
                        fw.write("\n")


        

                
                














if __name__ == "__main__":
    main()
