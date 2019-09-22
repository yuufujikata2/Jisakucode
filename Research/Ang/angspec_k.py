import sys
import math
from matplotlib import pyplot

def main():
    args = sys.argv
    argc = len(args)
    if argc != 2:
        print("error: please input a file!")
        sys.exit()
    allspec = []
    with open(args[1],mode="r") as fp:
        lines =fp.readlines()
        for i in range(1,len(lines)):
            spec = Spec(float(lines[i].split()[0]),
                        float(lines[i].split()[1]),
                        float(lines[i].split()[2]),
                        float(lines[i].split()[3]),
                        float(lines[i].split()[4]),
                        float(lines[i].split()[5]),
                        float(lines[i].split()[6]),
                        float(lines[i].split()[7]),
                        float(lines[i].split()[8]),
                        float(lines[i].split()[9]))
            allspec.append(spec)
    for i in range(len(allspec)):
        print("energy =",allspec[i].energy,
              "I_(z,z) =",allspec[i].zz,
              "I_(z,x) = ",allspec[i].zx,"...")

    print("")
    print("Please input range of angular(degree) (ex:0 90)")
    angular = input().split()
    startang = int(angular[0])
    endang = int(angular[1])
    print("please input step of angular point")
    angstep = int(input())
    
    fw_x = open("angspec_k_x.dat",mode = "w")
    fw_x.write("#")
    fw_x.write(" omega, I_shita")
    fw_x.write(" ( sita(degree) = ")
    fw_y = open("angspec_k_y.dat",mode = "w")
    fw_y.write("#")
    fw_y.write(" omega, ,I_shita")
    fw_y.write(" ( sita(degree) = ")
    fw_av = open("angspec_k_av.dat",mode = "w")
    fw_av.write("#")
    fw_av.write(" omega, ,I_shita")
    fw_av.write(" ( sita(degree) = ")
    fw_u = open("spec_k_unpol.dat",mode = "w")
    fw_u.write("#")
    fw_u.write(" omega I Izz Ixx Iyy")
    fw_u.write("\n")
    for ang in range(startang, endang + 1, angstep):
        fw_x.write(str(ang)+ " ")   
        fw_y.write(str(ang)+ " ")   
        fw_av.write(str(ang)+ " ")   
    fw_x.write(")\n")
    fw_y.write(")\n")
    fw_av.write(")\n")
    for ispec in allspec:
        fw_x.write(str("{:>11}".format("{:.6f}".format(ispec.energy))))
        fw_y.write(str("{:>11}".format("{:.6f}".format(ispec.energy))))
        fw_av.write(str("{:>11}".format("{:.6f}".format(ispec.energy))))
        fw_u.write(str("{:>11}".format("{:.6f}".format(ispec.energy))))
        fw_u.write(str("{:>11.6f}".format(ispec.xx+ispec.yy+ispec.zz)))
        fw_u.write(str("{:>11.6f}".format(ispec.zz)))
        fw_u.write(str("{:>11.6f}".format(ispec.xx)))
        fw_u.write(str("{:>11.6f}".format(ispec.yy)))
        fw_u.write("\n")
        for ang in range(startang, endang + 1, angstep):
            newspec_x = ispec.specchange_x(ang,fw_x)
            newspec_y = ispec.specchange_y(ang,fw_y)
            fw_av.write(str("{:>11}".format("{:.6f}".format((newspec_x + newspec_y) / 2))))
        fw_x.write("\n")
        fw_y.write("\n")
        fw_av.write("\n")
    fw_x.close()
    fw_y.close()
    fw_av.close()
    fw_u.close()

    pyplot.figure()

    for i in range((endang - startang ) // angstep + 1):
        angallspec_x = []
        angallspec_y = []
        angallspec_av = []
        angallenergy = []
        for ispec in allspec:
            angallspec_x.append(ispec.angspec_x[i])
            angallspec_y.append(ispec.angspec_y[i])
            angallspec_av.append((ispec.angspec_x[i] + ispec.angspec_y[i]) / 2)
            angallenergy.append(ispec.energy)
#        pyplot.plot(angallenergy,angallspec_x)
#        pyplot.plot(angallenergy,angallspec_y)
        pyplot.plot(angallenergy,angallspec_av)

    pyplot.savefig("angspec.png")
    pyplot.show()
            


class Spec():
    def __init__(self,energy,zz,zx,zy,xz,xx,xy,yz,yx,yy):
        self.energy = energy
        self.zz = zz
        self.zx = zx
        self.zy = zy
        self.xz = xz
        self.xx = xx
        self.xy = xy
        self.yz = yz
        self.yx = yx
        self.yy = yy
        self.angspec_x = []
        self.angspec_y = []
    def specchange_x(self,ang,fw):
        rad = math.radians(ang)
        newspec = math.cos(rad)**2 * self.zz + math.sin(rad)**2 * self.xx + (self.xz + self.zx)* math.sin(rad) * math.cos(rad)
        fw.write(str("{:>11}".format("{:.6f}".format(newspec))))
        self.angspec_x.append(newspec)
        return newspec
    def specchange_y(self,ang,fw):
        rad = math.radians(ang)
        newspec = math.cos(rad)**2 * self.zz + math.sin(rad)**2 * self.yy + (self.yz + self.zy)* math.sin(rad) * math.cos(rad)
        fw.write(str("{:>11}".format("{:.6f}".format(newspec))))
        self.angspec_y.append(newspec)
        return newspec

if __name__ == "__main__":
    main()
