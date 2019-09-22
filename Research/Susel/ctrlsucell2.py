"""
To make super cell

call fortran sub routine

"""
import sys
import numpy as np
import copy
from cProfile import Profile
import time

BOHR=0.529177

def main():
    fp = open("CTRL", mode ="r")
    lines = fp.readlines()
    site_t=None
    class_t=None
    plat_n=0
    for line in lines:
      word = line.split()
      if word[0]=="STRUC":
        worde = line.split("=")
        alat=float(worde[1].strip())
        plat_n=4
      if word[0]=="DIM":
        for w in word:
          worde = w.split("=")
          if worde[0]=="NBAS":
            nbas=int(worde[1])
          if worde[0]=="NCLASS":
            nclass=int(worde[1])
      if plat_n==3:
        Plat=np.zeros((3,3))
        worde = line.split("=")
        worde2 = worde[1].split()
        Plat[0][0]=float(worde2[0])
        Plat[0][1]=float(worde2[1])
        Plat[0][2]=float(worde2[2])
      elif plat_n==2:
        Plat[1][0]=float(word[0])
        Plat[1][1]=float(word[1])
        Plat[1][2]=float(word[2])
      elif plat_n==1:
        Plat[2][0]=float(word[0])
        Plat[2][1]=float(word[1])
        Plat[2][2]=float(word[2])
      plat_n-=1
      if word[0]=="CLASS":
        class_t=True
        class_n=nclass
        allclass_ctrl=[]
      if class_t:
        if len(word)<5:
          continue
        worde=line.split("=")
        class_ctrl=Class_ctrl(worde[1].split()[0],worde[2].split()[0],worde[3].split()[0])
        allclass_ctrl.append(class_ctrl)
        class_n-=1
        if class_n==0:
          class_t=False
      if word[0]=="SITE":
        site_t=True
        site_n=nbas
        allatom=[]
      if site_t:
        atom=Atom()
        worde = line.split("=")
        atom.name=worde[1].split()[0]
        atom.pos=np.array([float(worde[2].split()[0]),float(worde[2].split()[1]),float(worde[2].split()[2])])
        allatom.append(atom)
        site_n-=1
        if site_n==0:
          site_t=False
      if site_t==False:
        break
  
  
  
    for i in allatom:
      i.ctrltuika(allclass_ctrl)  
  
    print ("\nallatom\n")
  
    for i in range (len(allatom)):
      print (allatom[i].name,"iz=",allatom[i].ze,"rad=",allatom[i].rad,"pos=",allatom[i].pos)
  
    print("")
  
    nbas=len(allatom)
    
    
     
    print("How is a supercell size? (plese input a b c)")
    sucellsize=input().split()
    if len(sucellsize)!=3:
      print("error : plese input 3 numbers")
      sys.exit()
  
    sa=int(sucellsize[0])
    sb=int(sucellsize[1])
    sc=int(sucellsize[2])
  
  
    sinPlat=np.zeros((3,3))
  
    sinPlat[0]=Plat[0]*sa
    sinPlat[1]=Plat[1]*sb
    sinPlat[2]=Plat[2]*sc
       
    
    cart=np.zeros((nbas,3))
    allze = [[None] * 2 for i in range(nbas)]
    for i in range(nbas):
      allze[i][0] = allatom[i].ze
      allze[i][1] = allatom[i].name.strip("0123456789")
  
    allze = tuple(set(map(tuple,allze))) 
    zedict = {}
    for i in allze:
        zedict[i] = 0
  
  
    sinallatom_ctrl = []
    sinallatom_xyz = []
  
   
    for i0 in range(sa):
      for i1 in range(sb):
        for i2 in range(sc):
          for k in range(nbas):
            atom_ctrl = Atom()
            atom_xyz  = Atom()
            atom_ctrl.pos = i0*Plat[0] + i1*Plat[1] + i2*Plat[2] + allatom[k].pos
            atom_xyz.pos = atom_ctrl.pos * alat * BOHR
            atom_ctrl.ze = allatom[k].ze
            atom_ctrl.rad = allatom[k].rad
            for zename, zenumber in zedict.items():
                if atom_ctrl.ze == int(zename[0]):
                    atom_ctrl.name = zename[1] + str(zenumber)
                    zedict[zename] += 1
            atom_xyz.name = allatom[k].name
  
            sinallatom_ctrl.append(atom_ctrl)
            sinallatom_xyz.append(atom_xyz)
  
  
    for i in sinallatom_ctrl:
        print(i.name, i.ze, i.rad, i.pos)
  
    with open("susel.xsf", mode = "w") as fw:
        fw.write("CRYSTAL\n")
        fw.write("\n")
        fw.write("PRIMVEC\n")
        for i in range (3):
            for j in range(3):
                fw.write("{:>11.6f}".format(sinPlat[i][j] * alat * BOHR))
            fw.write("\n")
        fw.write("\nPRIMCOORD\n")
        fw.write(str(len(sinallatom_xyz))+" 1\n")
        for i in sinallatom_xyz:
            fw.write("{:>5}".format(i.name))
            for j in i.pos:
                fw.write("{:>11.6f}".format(j))
            fw.write("\n")
            

    with open ("CTRL_su", mode = "w") as fw:
        for i in range (3):
            fw.write(lines[i])
        fw.write("SYMGRP    NGEN=1 GENGRP=E\n")
        fw.write("          SPCGRP=undefined USESYM=F\n")
        fw.write("STRUC     ALAT=" + str(alat)+"\n")
        fw.write("          PLAT=")
        for i in range (3):
            fw.write("{:<12.8f}".format(sinPlat[0][i]))
        fw.write("\n")
        fw.write("               ")
        for i in range (3):
            fw.write("{:<12.8f}".format(sinPlat[1][i]))
        fw.write("\n")
        fw.write("               ")
        for i in range (3):
            fw.write("{:<12.8f}".format(sinPlat[2][i]))
        fw.write(" FIXLAT=T\n")

        for i in range(len(sinallatom_ctrl)):
            if i == 0:
                fw.write("CLASS     ")
            else:
                fw.write("          ")
            fw.write("ATOM=" + "{:<6}".format(sinallatom_ctrl[i].name))
            fw.write(" Z=" + "{:<3}".format(sinallatom_ctrl[i].ze))
            fw.write(" R=" + "{:<13.8f}".format(sinallatom_ctrl[i].rad))
            fw.write("\n")
        for i in range(len(sinallatom_ctrl)):

            if i == 0:
                fw.write("SITE      ")
            else:
                fw.write("          ")
            fw.write("ATOM=" + "{:<6}".format(sinallatom_ctrl[i].name))
            fw.write("POS=")
            for j in sinallatom_ctrl[i].pos:
                fw.write("{:<13.8f}".format(j))
            fw.write("\n")
        for line in lines:
            if line.split()[0] == "SCALE" or line.split()[0] == "STR":
                fw.write(line)




    fp.close()


  
    print("Succeeded!")



class Class_ctrl():
  def __init__(self,name,ze,rad):
    self.name=name
    self.ze=int(ze)
    self.rad=float(rad)


class Atom():
  def __init__(self):
    self.name=None
    self.pos=None
    self.ze=None
    self.rad=None

  def ctrltuika(self,allclass_ctrl):
    for i in allclass_ctrl:
      if i.name==self.name:
        self.ze=i.ze
        self.rad=i.rad
        break



if __name__=="__main__":
  main()
