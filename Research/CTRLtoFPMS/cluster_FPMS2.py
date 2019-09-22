"""
To create FPMS format for some material from CTRL(lmtofile) and cluster.xyz 

"""
import sys
import numpy as np
import copy
BOHR=0.529177

def main():
  
  args=sys.argv
  argc=len(args)
  if argc >=5:
    print("error")
    sys.exit()
  if argc==4:
    fpxyz = args[1]
    fpctrl = args[2]
    scalefac = float(args[3])
  elif argc==3:
    fpxyz = args[1]
    fpctrl = args[2]
    print("Please input scalefactor \n(radius of output =　Real radius(Angs) * scalefactor)")
    scalefac = float(input())  
  elif argc==2:
    fpxyz = args[1]
    fpctrl="CTRL" 
    print("Is there CTRL ? (y/n)")
    if input()!="y":
      sys.exit()
    print("Please input scalefactor \n(radius of output =　Real radius(Angs) * scalefactor)")
    scalefac = float(input())  
  elif argc==1:
    print("Please input xyzfile")
    fpxyz = input()
    print("Please input CTRL")
    fpctrl = input()
    print("Please input scalefactor \n(radius of output =　Real radius(Angs) * scalefactor)")
    scalefac = float(input())  

  fw1="clus.fp"
  fw2="clus2.fp"
  
  with open (fpxyz,mode="r") as fp:
    line = fp.readlines()
    number = int(line[0]) 
    allatom_xyz=[]
    for i in range(2,number+2):
      word=line[i].split()
      atom_xyz=Atom_xyz(word[0],word[1],word[2],word[3])
      allatom_xyz.append(atom_xyz)
#      Atom.append(word[0])
#      Posatom=[]
#      Posatom.append(float(word[1]))
#      Posatom.append(float(word[2]))
#      Posatom.append(float(word[3]))
#      Pos.append(Posatom)
      
  with open (fpctrl,mode="r") as fp:
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
        Atom_ctrl=[]
        Pos_ctrl=[]
      if site_t:
        if site_n==nbas:
          Pos_1=[]
          for w in word:
            worde = w.split("=")
            if worde[0]=="ATOM":
              Atom_ctrl.append(worde[1])
          if len(word)==5:
            worde = word[2].split("=")
            Pos_1.append(float(worde[1]))
            Pos_1.append(float(word[3]))
            Pos_1.append(float(word[4]))
          elif len(word)==6:
            Pos_1.append(float(word[3]))
            Pos_1.append(float(word[4]))
            Pos_1.append(float(word[5]))
          else:
            print("error")
            sys.exit()
          Pos_ctrl.append(Pos_1)
        else:
          Pos_1=[]    
          Atom_ctrl.append(word[0].split("=")[1])
          if len(word)==4:
            worde = word[1].split("=")
            Pos_1.append(float(worde[1]))
            Pos_1.append(float(word[2]))
            Pos_1.append(float(word[3]))
          elif len(word)==5:
            Pos_1.append(float(word[2]))
            Pos_1.append(float(word[3]))
            Pos_1.append(float(word[4]))
          else:
            print("error")
            sys.exit()
          Pos_ctrl.append(Pos_1)
        site_n-=1
        if site_n==0:
          site_t=False
      if site_t==False:
        break
  for i in range (len(allclass_ctrl)):
    print (allclass_ctrl[i].name,allclass_ctrl[i].ze,allclass_ctrl[i].rad)
  for i in range(len(allatom_xyz)):
    print(allatom_xyz[i].name,allatom_xyz[i].pos)

  for i in allatom_xyz:
    i.ctrltuika(allclass_ctrl)  
  
  for i in allatom_xyz:
    print(i.name,i.ze,i.pos,i.rad) 

  print("What is shape ? (NM or MT)")
  shape=input() 

  if shape=="NM":
    with open(fw1,mode='w') as fw:
      linenumber=1
      for line in allatom_xyz:
        fw.write(str("{:>4}".format(line.ze)))
        fw.write(str('{:>12}'.format("{:6f}".format(line.pos[0]))))
        fw.write(str('{:>11}'.format("{:6f}".format(line.pos[1]))))
        fw.write(str('{:>11}'.format("{:6f}".format(line.pos[2]))))
        fw.write(str('{:>12}'.format("{:8f}".format(line.rad*BOHR*scalefac))))
        fw.write(str('{:>7}'.format(str(linenumber))))
        fw.write(str('{:>11.6f}'.format(np.sqrt((line.pos[0]-allatom_xyz[0].pos[0])**2+(line.pos[1]-allatom_xyz[0].pos[1])**2+(line.pos[2]-allatom_xyz[0].pos[2])**2))))
        fw.write(str('{:>5}'.format(line.name)))
        fw.write("\n")
        linenumber+=1
    
    with open(fw2,mode="w") as fw:
      allatom_xyz2=copy.deepcopy(allatom_xyz)
      linenumber=1
      sortnumber=0
      sortze=[]
      sortname=[]
      baburu_t=True
      while(baburu_t):
        baburu=0
        baburu_t=False
        for i in range(0,len(allatom_xyz2)-1-baburu):
          if allatom_xyz2[i].ze<allatom_xyz2[i+1].ze:
            kari=allatom_xyz2[i]
            allatom_xyz2[i]=allatom_xyz2[i+1]
            allatom_xyz2[i+1]=kari 
            baburu_t=True
      while(sortnumber<len(allatom_xyz2)):
        if not (allatom_xyz2[sortnumber].name in sortname):
          sortname.append(allatom_xyz2[sortnumber].name)
          for line in allatom_xyz2:
            if line.name==sortname[-1]:
              print(line.name)
              fw.write(str("{:>4}".format(line.ze)))
              fw.write(str('{:>12}'.format("{:6f}".format(line.pos[0]))))
              fw.write(str('{:>11}'.format("{:6f}".format(line.pos[1]))))
              fw.write(str('{:>11}'.format("{:6f}".format(line.pos[2]))))
              fw.write(str('{:>12}'.format("{:8f}".format(line.rad*BOHR*scalefac))))
              fw.write(str('{:>7}'.format(str(linenumber))))
              fw.write(str('{:>11.6f}'.format(np.sqrt((line.pos[0]-allatom_xyz2[0].pos[0])**2+(line.pos[1]-allatom_xyz2[0].pos[1])**2+(line.pos[2]-allatom_xyz2[0].pos[2])**2))))
              fw.write(str('{:>5}'.format(line.name)))
              fw.write("\n")
              linenumber+=1
        sortnumber+=1
      
      
      
  elif shape=="MT":
    with open(fw1,mode='w') as fw:
      linenumber=1
      for line in allatom_xyz:
        if line.ze==0:
          continue
        fw.write(str("{:>4}".format(line.ze)))
        fw.write(str('{:>12}'.format("{:6f}".format(line.pos[0]))))
        fw.write(str('{:>11}'.format("{:6f}".format(line.pos[1]))))
        fw.write(str('{:>11}'.format("{:6f}".format(line.pos[2]))))
        fw.write(str('{:>7}'.format(str(linenumber))))
        fw.write(str('{:>11.6f}'.format(np.sqrt((line.pos[0]-allatom_xyz[0].pos[0])**2+(line.pos[1]-allatom_xyz[0].pos[1])**2+(line.pos[2]-allatom_xyz[0].pos[2])**2))))
        fw.write(str('{:>5}'.format(line.name)))
        fw.write("\n")
        linenumber+=1
        number+=1
    
    with open(fw2,mode="w") as fw:
      allatom_xyz2=copy.deepcopy(allatom_xyz)
      linenumber=1
      sortnumber=0
      sortze=[]
      sortname=[]
      baburu_t=True
      while(baburu_t):
        baburu=0
        baburu_t=False
        for i in range(0,len(allatom_xyz2)-1-baburu):
          if allatom_xyz2[i].ze<allatom_xyz2[i+1].ze:
            kari=allatom_xyz2[i]
            allatom_xyz2[i]=allatom_xyz2[i+1]
            allatom_xyz2[i+1]=kari 
            baburu_t=True
      while(sortnumber<len(allatom_xyz2)):
        if not (allatom_xyz2[sortnumber].name in sortname):
          sortname.append(allatom_xyz2[sortnumber].name)
          for line in allatom_xyz2:
            if line.ze==0:
             continue
            if line.name==sortname[-1]:
              print(line.name)
              fw.write(str("{:>4}".format(line.ze)))
              fw.write(str('{:>12}'.format("{:6f}".format(line.pos[0]))))
              fw.write(str('{:>11}'.format("{:6f}".format(line.pos[1]))))
              fw.write(str('{:>11}'.format("{:6f}".format(line.pos[2]))))
              fw.write(str('{:>7}'.format(str(linenumber))))
              fw.write(str('{:>11.6f}'.format(np.sqrt((line.pos[0]-allatom_xyz2[0].pos[0])**2+(line.pos[1]-allatom_xyz2[0].pos[1])**2+(line.pos[2]-allatom_xyz2[0].pos[2])**2))))
              fw.write(str('{:>5}'.format(line.name)))
              fw.write("\n")
              linenumber+=1
        sortnumber+=1

  print("Succeeded!")

class Atom_xyz():
  def __init__(self,name,x,y,z):
    self.name=name
    self.pos=[float(x),float(y),float(z)]
    self.ze=None
    self.rad=None

  def ctrltuika(self,allclass_ctrl):
    for i in allclass_ctrl:
      if i.name==self.name:
        self.ze=i.ze
        self.rad=i.rad
        break

class Class_ctrl():
  def __init__(self,name,ze,rad):
    self.name=name
    self.ze=int(ze)
    self.rad=float(rad)

if __name__=="__main__":
  main()
