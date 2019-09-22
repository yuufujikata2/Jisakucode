"""
To calculate cluster.xyz from CTRL(lmtofile)

call fortran sub routine

"""
import sys
import numpy as np
import copy
from cProfile import Profile
import time
from fortcall_xpd import Fort

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
#  print(alat,nbas)
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
  
  print("which is surface ?(Please input Miller index)")
  millerin=input().split()
  miller=np.array([0.,0.,0.])
  if len(millerin)==3:
    for i in range(3):
      miller[i]=float(millerin[i])
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

    clin.write(str(" {:f}".format(miller[0]))+str(" {:f}".format(miller[1]))+str(" {:f}".format(miller[2]))+"\n")
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[0].pos[0]))))
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[0].pos[1]))))
    clin.write(str("{:>12}".format("{:.8f}".format(allatom[0].pos[2]))))
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

#  sys.exit()
  cart=np.zeros((nbas,3))
  alliz=np.zeros((nbas),dtype="int32")
  for i in range(nbas):
    alliz[i]=allatom[i].ze
    for j in range(3):
      cart[i][j]=allatom[i].pos[j]
  
  fort=Fort()
  fort.call_fortran(clradius,Plat,nbas,cart,alat,alliz,miller)

  print("Succeeded!")
 
  """
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
  """
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
