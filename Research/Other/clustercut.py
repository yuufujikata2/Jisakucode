import sys

BOHR = 0.529177

def main():
    args=sys.argv
    argc=len(args)
    if argc != 2:
        print("error Please input a xyzfile!")
        sys.exit()

    with open(args[1],mode="r") as fr:
        lines=fr.readlines()
        natoms = int(lines[0].split()[0])
        allatom=[]
        for i in range (2, 2 + natoms):
            atom = Atom(lines[i].split()[0],float(lines[i].split()[1],),float(lines[i].split()[2]),float(lines[i].split()[3]))
            allatom.append(atom)
    
    print("Please input cut direction ( x or y or z or -x ...)")
    direction = input()
    print("Plase input cutting height")
    cuthigh = float(input())

    for i in range(natoms-1,-1,-1):
        if direction == "x":
            if allatom[i].x > cuthigh:
                allatom.pop(i)
        elif direction == "y":
            if allatom[i].y > cuthigh:
                allatom.pop(i)
        elif direction == "z":
            if allatom[i].z > cuthigh:
                allatom.pop(i)
        elif direction == "-x":
            if allatom[i].x < cuthigh:
                allatom.pop(i)
            else:
                allatom[i].x = (-1) * allatom[i].x
        elif direction == "-y":
            if allatom[i].y < cuthigh:
                allatom.pop(i)
            else:
                allatom[i].y = (-1) * allatom[i].y
        elif direction == "-z":
            if allatom[i].z < cuthigh:
                allatom.pop(i)
            else:
                allatom[i].z = (-1) * allatom[i].z

    with open ("cutcluster.xyz",mode="w") as fw:
        fw.write(" " + str(len(allatom)) + "\n")
        fw.write("\n")
        for i in range(len(allatom)):
            fw.write(str("{:<4}".format(allatom[i].name)))
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].x))))
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].y))))
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].z))))
            fw.write("\n")
   
    lmaxatoms = []
    for i in allatom:
        count = 0
        for j in lmaxatoms:
            if j.name == i.name:
                count += 1
        if count == 0:
            lmaxatoms.append(i)

    print("please input iz and lmax")
    for i in range(len(lmaxatoms)):
        print(lmaxatoms[i].name)
        lmaxatoms[i].iz = input("iz:")
        lmaxatoms[i].lmax = input("lmax:")

    for i in range(len(allatom)):
        for j in lmaxatoms:
            if allatom[i].name == j.name:
                allatom[i].lmax = j.lmax
                allatom[i].iz = j.iz
        
    with open ("instr",mode="w") as fw:
        fw.write(str("{:>5}".format(len(allatom))))
        fw.write(str("{:>5}".format(len(allatom))))
        fw.write(str("{:>5}".format(2)))
        fw.write("\n")
        for i in range(len(allatom)):
            fw.write(str("{:<4}".format(allatom[i].name))) 
            fw.write(str("{:>2}".format(allatom[i].lmax))) 
            fw.write(str("{:>3}".format(allatom[i].iz))) 
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].x / BOHR))))
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].y / BOHR))))
            fw.write(str("{:>11}".format("{:6f}".format(allatom[i].z / BOHR))))
            fw.write("\n")

    print("succeeded!")

class Atom():
    def __init__(self,name,x,y,z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.lmax = None
        self.iz = None

if __name__ == "__main__":
    main()

