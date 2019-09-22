import sys

def main():
    args=sys.argv
    argc=len(args)
    if argc !=3:
        print("error: please input 2 files!")
        sys.exit()
    with open(args[1],mode="r") as fi:
        lines=fi.readlines()
        subject=lines[0].strip()
        atomchar=lines[5].split()
        atomnumber=lines[6].split()
        allnumber=0
        for i in atomnumber:
            allnumber += int(i)
        allcood=[]
        if lines[7].split()[0] == "Selective":
            for line in lines[9:8+allnumber]:
                cood=[0,0,0]
                for i in range(3):
                    cood[i]=float(line.split()[i])
                allcood.append(cood)
        else:
            for line in lines[8:7+allnumber]:
                cood=[0,0,0]
                for i in range(3):
                    cood[i]=float(line.split()[i])
                allcood.append(cood)

    for i in range (len(atomchar)):
        for j in range(int(atomnumber[i])):
            print(str(atomchar[i]) + " " + str(allcood[j][0]) + " " + str(allcood[j][1]) + " " + str(allcood[j][2]))
    print()
    print("please input calculate condition")
    calccon=input(">>")
    print("Plase input whether coodinate moves. 0 or 1")
    xmove=input("x:")
    ymove=input("y:")
    zmove=input("z:")
    with open(args[2],mode="w") as fw:
        fw.write(calccon)
        fw.write("\n")
        fw.write(subject)
        fw.write("\n")
        fw.write("\n")
        for i in range (len(atomchar)):
            for j in range(int(atomnumber[i])):
                fw.write(str("{:>2}".format(atomchar[i])))
                fw.write(str("{:>12}".format("{:6f}".format(allcood[j][0]))))
                fw.write(str("{:>2}".format(xmove)))
                fw.write(str("{:>12}".format("{:6f}".format(allcood[j][1]))))
                fw.write(str("{:>2}".format(ymove)))
                fw.write(str("{:>12}".format("{:6f}".format(allcood[j][2]))))
                fw.write(str("{:>2}".format(zmove)))
                fw.write("\n")
    print("succeeded!")
                

if __name__ == "__main__" :
    main()

            


