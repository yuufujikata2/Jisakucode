"""
To check symmetry

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
    
    
    print("Please input rotation axis (input 2 atoms)")
    centeratom=input().split()
    c=-1
    origincart=np.array([0.,0.,0.])
    if len (centeratom) ==1:
      for i in range(nbas):
        if centeratom[0]==allatom[i].name:
          c=i
          break
      origincart=allatom[c].pos
      if c==-1:
        print("error : There is no atom of that name")
        sys.exit()
    elif len(centeratom)==3:
      for i in range(3):
        origincart[i]=float(centeratom[i])
    else:
      print("error")
      sys.exit()

    headatom = input().split()
    headcart=np.array([0.,0.,0.])
    if len (headatom) ==1:
      for i in range(nbas):
        if headatom[0]==allatom[i].name:
          c=i
          break
      headcart=allatom[c].pos
      if c==-1:
        print("error : There is no atom of that name")
        sys.exit()
    elif len(headatom)==3:
      for i in range(3):
        headcart[i]=float(headatom[i])
    else:
      print("error")
      sys.exit()

    print("origincart = ",origincart)
    print("headcart = ",headcart)


    print("Please input times of rotation")
    times = int(input())
    gamma = 2 * np.pi / times    

    #atom position shift
    headcart_2 = headcart - origincart

    allatom_2 = copy.deepcopy(allatom)
    for i in allatom_2 :
      i.pos = i.pos - origincart

    #alpha shita
    if np.sqrt(headcart_2[1]**2 +headcart_2[0]**2) == 0:
      alpha = 0
    elif headcart_2[1] >= 0:
        alpha = np.arccos(headcart_2[0] / np.sqrt(headcart_2[1]**2 +headcart_2[0]**2))
    else:
        alpha = -np.arccos(headcart_2[0] / np.sqrt(headcart_2[1]**2 +headcart_2[0]**2))
    if np.sqrt(headcart_2[1]**2 + headcart_2[0]**2 + headcart_2[2]**2) == 0:
      shita = 0
    else:
      shita = np.arccos(headcart_2[2] / np.sqrt(headcart_2[1]**2 + headcart_2[0]**2 + headcart_2[2]**2))

    print("alpha = ",alpha)
    print("shita = ",shita)

    #rotation matrix
    matrix_a = np.array([[np.cos(-alpha), -np.sin(-alpha), 0.],
                         [np.sin(-alpha),  np.cos(-alpha), 0.],
                         [0.,             0.,            1.]])

    matrix_s = np.array([[np.cos(-shita) , 0., np.sin(-shita)],
                         [0.,              1., 0.            ],
                         [-np.sin(-shita), 0., np.cos(-shita)]])
     
    matrix_g = np.array([[np.cos(gamma), -np.sin(gamma), 0.],
                         [np.sin(gamma),  np.cos(gamma), 0.],
                         [0.,             0.,            1.]])
    matrix_ai = np.array([[np.cos(alpha), -np.sin(alpha), 0.],
                         [np.sin(alpha),  np.cos(alpha), 0.],
                         [0.,             0.,            1.]])

    matrix_si = np.array([[np.cos(shita) , 0., np.sin(shita)],
                         [0.,              1., 0.            ],
                         [-np.sin(shita), 0., np.cos(shita)]])
     
    matrix_gi = np.array([[np.cos(-gamma), -np.sin(-gamma), 0.],
                         [np.sin(-gamma),  np.cos(-gamma), 0.],
                         [0.,             0.,            1.]])
    """
    print(matrix_a)
    print(matrix_s)
    print(matrix_g)
    """

    #rotation
    allatom_3 = copy.deepcopy(allatom_2)

    headcart_2 = np.dot(matrix_a, headcart_2)
    headcart_2 = np.dot(matrix_s, headcart_2)
    
    for i in allatom_3:
      i.pos = np.dot(matrix_a, i.pos)
      i.pos = np.dot(matrix_s, i.pos)

    allatom_4 = copy.deepcopy(allatom_3)


    for t in range(times - 1):
      for i in allatom_4:
        i.pos = np.dot(matrix_g, i.pos)
  
      for i in allatom_4:
        print(i.pos)
  
      for i in range(len(allatom_3)):
        for j in range(1,len(allatom_4)):
          diff_x = abs(allatom_3[i].pos[0] - allatom_4[j].pos[0])
          diff_y = abs(allatom_3[i].pos[1] - allatom_4[j].pos[1])
          diff_z = abs(allatom_3[i].pos[2] - allatom_4[j].pos[2])
          if diff_x < 0.0001 and diff_y < 0.0001 and diff_z < 0.0001:
            allatom_3[i].eqatomname.append(allatom_4[j].name)

    #bunnrui
    allatom_5 = []

    for i in allatom_3:
      for j in i.eqatomname:
        for k in allatom_5:
          if k.name == j:
            i.name = k.name
      allatom_5.append(i)

    #modoshi
    for i in range(len(allatom_5)):
      allatom_5[i].pos = np.dot(matrix_si,allatom_5[i].pos)
      allatom_5[i].pos = np.dot(matrix_ai,allatom_5[i].pos)
      allatom_5[i].pos = allatom_5[i].pos + origincart


    allclass = []
    for i in allatom_5:
        a = 0
        for j in allclass:
            if i.name == j.name:
                a += 1
        if a == 0:
            allclass.append(i)

    for i in allatom_5:
      print(i.name,i.pos)

    #for i in allatom_3:
    #  print(i.eqatomname)
    

    with open("susel.xsf", mode = "w") as fw:
        fw.write("CRYSTAL\n")
        fw.write("\n")
        fw.write("PRIMVEC\n")
        for i in range (3):
            for j in range(3):
                fw.write("{:>11.6f}".format(Plat[i][j] * alat * BOHR))
            fw.write("\n")
        fw.write("\nPRIMCOORD\n")
        fw.write(str(len(allatom_5))+" 1\n")
        for i in allatom_5:
            fw.write("{:>5}".format(i.name))
            for j in i.pos:
                fw.write("{:>11.6f}".format(j))
            fw.write("\n")
            

    with open ("CTRL_sym", mode = "w") as fw:
        for i in range (3):
            fw.write(lines[i])
        fw.write("SYMGRP    NGEN=1 GENGRP=E\n")
        fw.write("          SPCGRP=undefined USESYM=F\n")
        fw.write("STRUC     ALAT=" + str(alat)+"\n")
        fw.write("          PLAT=")
        for i in range (3):
            fw.write("{:<12.8f}".format(Plat[0][i]))
        fw.write("\n")
        fw.write("               ")
        for i in range (3):
            fw.write("{:<12.8f}".format(Plat[1][i]))
        fw.write("\n")
        fw.write("               ")
        for i in range (3):
            fw.write("{:<12.8f}".format(Plat[2][i]))
        fw.write(" FIXLAT=T\n")

        for i in range(len(allclass)):
            if i == 0:
                fw.write("CLASS     ")
            else:
                fw.write("          ")
            fw.write("ATOM=" + "{:<6}".format(allclass[i].name))
            fw.write(" Z=" + "{:<3}".format(allclass[i].ze))
            fw.write(" R=" + "{:<13.8f}".format(allclass[i].rad))
            fw.write("\n")
        for i in range(len(allatom_5)):

            if i == 0:
                fw.write("SITE      ")
            else:
                fw.write("          ")
            fw.write("ATOM=" + "{:<6}".format(allatom_5[i].name))
            fw.write("POS=")
            for j in allatom_5[i].pos:
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
    self.eqatomname = []

  def ctrltuika(self,allclass_ctrl):
    for i in allclass_ctrl:
      if i.name==self.name:
        self.ze=i.ze
        self.rad=i.rad
        break



if __name__=="__main__":
  main()
