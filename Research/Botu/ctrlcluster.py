"""
To create cluster.in from CTRL(lmtofile)

"""
import sys
import numpy as np
import copy
from cProfile import Profile
import time

BOHR=0.529177
FACTOR=5
def main():
  with open ("CTRL",mode="r") as fp:
    site_t=None
    class_t=None
    plat_n=0
    for line in fp:
      word = line.split()
      if word[0]=="STRUC":  
        for w in word:
          worde = w.split("=")
          if worde[0]=="ALAT":
            alat=float(worde[1])
        plat_n=4
      if word[0]=="DIM":
        for w in word:
          worde = w.split("=")
          if worde[0]=="NBAS":
            nbas=int(worde[1])
          if worde[0]=="NCLASS":
            nclass=int(worde[1])
#  print(alat,nbas)
      if plat_n==3:
        Plat=np.zeros((3,3))
        if len(word)==3:
          worde = word[0].split("=")
          Plat[0][0]=float(worde[1])
          Plat[0][1]=float(word[1])
          Plat[0][2]=float(word[2])
        else:
          Plat[0][0]=float(word[1])
          Plat[0][1]=float(word[2])
          Plat[0][2]=float(word[3])
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
        if class_n==nclass:
          if word[2]=="Z=":
            class_ctrl=Class_ctrl(word[1].split("=")[1],word[3],word[4].split("=")[1])
          else: 
            class_ctrl=Class_ctrl(word[1].split("=")[1],word[2].split("=")[1],word[3].split("=")[1])
          allclass_ctrl.append(class_ctrl)
        else:
          if word[1]=="Z=":
            class_ctrl=Class_ctrl(word[0].split("=")[1],word[2],word[3].split("=")[1])
          else: 
            class_ctrl=Class_ctrl(word[0].split("=")[1],word[1].split("=")[1],word[2].split("=")[1])
          allclass_ctrl.append(class_ctrl)
        class_n-=1  
        if class_n==0:
          class_t=False  
      if word[0]=="SITE":
        site_t=True
        site_n=nbas 
        allatom=[]
      if site_t:
        if site_n==nbas:
          atom=Atom()
          for w in word:
            worde = w.split("=")
            if worde[0]=="ATOM":
              atom.name=worde[1]
          if len(word)==5:
            worde = word[2].split("=")
            atom.pos=np.array([float(worde[1]),float(word[3]),float(word[4])])
          elif len(word)==6:
            atom.pos=np.array([float(word[3]),float(word[4]),float(word[5])])
          else:
            print("error")
            sys.exit()
          allatom.append(atom)
        else:
          atom=Atom()
          atom.name=(word[0].split("=")[1])
          if len(word)==4:
            worde = word[1].split("=")
            atom.pos=np.array([float(worde[1]),float(word[2]),float(word[3])])
          elif len(word)==5:
            atom.pos=np.array([float(word[2]),float(word[3]),float(word[4])])
          else:
            print("error")
            sys.exit()
          allatom.append(atom)
        site_n-=1
        if site_n==0:
          site_t=False
      if site_t==False:
        break

  for i in allatom:
    i.ctrltuika(allclass_ctrl)  
  for i in range (len(allatom)):
    print (allatom[i].name,allatom[i].ze,allatom[i].rad,allatom[i].pos)

#  print(Plat)
#  sys.exit()
#  Pos_np=np.array(Pos)
#  print(Atom)
#  print(Pos_np)   
#  Plat=Plat*alat*BOHR
#  Pos_np=Pos_np*alat*BOHR
#  print(Pos_np)   
  
  print("which is center ?(Atomname or Coordinate)")
  centeratom=input().split()
  c=0
  origincart=np.array([0.,0.,0.])
  if len (centeratom) ==1:
    for i in range(nbas):
      if centeratom==allatom[i].name:
        c=i
        break
    origincart=allatom[c].pos
  elif len(centeratom)==3:
    for i in range(3):
      origincart[i]=float(centeratom[i])
  else:
    print("error")
    sys.exit()
    
  
   
  print("How is clustersize ? (please input cluster radius)")
  clradius=float(input()) 
  
  with open ("cluster.in",mode="w") as clin:
    for i in range(3):
      clin.write(str("{:>11}".format("{:.8f}".format(Plat[i][0])))) 
      clin.write(str("{:>11}".format("{:.8f}".format(Plat[i][1])))) 
      clin.write(str("{:>11}".format("{:.8f}".format(Plat[i][2])))) 
      clin.write("\n")

    clin.write(" 0 0 0\n")
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[c].pos[0]))))
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[c].pos[1]))))
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[c].pos[2]))))
    clin.write("\n")
    clin.write(" "+str(nbas)+"\n")
    for i in range(nbas):   
      clin.write(str("{:>3}".format(allatom[i].name)))
      clin.write(str("{:>12}".format("{:.8f}".format(allatom[i].pos[0]))))
      clin.write(str("{:>12}".format("{:.8f}".format(allatom[i].pos[1]))))
      clin.write(str("{:>12}".format("{:.8f}".format(allatom[i].pos[2]))))
      clin.write("\n")
    clin.write(str("{:.4f}".format(alat*BOHR))+"\n")
    for i in range(3):
      clin.write(str(clradius)+"\n")

  sys.exit()

  cutoff=clradius
  i0max=int(cutoff * FACTOR /np.sqrt(np.sum(Plat[0]**2)))
  i1max=int(cutoff * FACTOR /np.sqrt(np.sum(Plat[1]**2)))
  i2max=int(cutoff * FACTOR /np.sqrt(np.sum(Plat[2]**2)))

  clusteratom=[]
  y=np.array([0.,0.,0.])
  n=0  
  pr=Profile()
  pr.enable()
  t1=time.time()
  for i0 in range(-i0max,i0max+1):
    for i1 in range(-i1max,i1max+1):
      for i2 in range(-i2max,i2max+1):
        for k in range(nbas):
          y=i0*Plat[0] + i1*Plat[1] + i2*Plat[2] + allatom[k].pos - origincart
          y*=alat*BOHR
          if np.sqrt(sum(y**2))<clradius:
            atom=Atom()
            atom.name=allatom[k].name
            atom.ze=allatom[k].ze
            atom.pos=copy.deepcopy(y)
            clusteratom.append(atom)
            n+=1
  print(n)
  t2=time.time()
  print(t2-t1)
  sys.exit()
  pr.disable()
  pr.print_stats()
  sys.exit()
  for i in range(len(clusteratom)):
    print(clusteratom[i].name,clusteratom[i].ze,clusteratom[i].pos) 

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
