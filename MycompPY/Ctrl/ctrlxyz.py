import sys
import numpy as np
BOHR=0.529177

def main():
  with open ("CTRL",mode="r") as ctrl:
    site_t=None
    plat_n=0
    for line in ctrl:
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
      if word[0]=="SITE":
        site_t=True
        site_n=nbas 
        Atom=[]
        Pos=[]
      if site_t:
        if site_n==nbas:
          Pos_1=[]
          for w in word:
            worde = w.split("=")
            if worde[0]=="ATOM":
              Atom.append(worde[1])
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
          Pos.append(Pos_1)
        else:
          Pos_1=[]    
          Atom.append(word[0].split("=")[1])
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
          Pos.append(Pos_1)
        site_n-=1
        if site_n==0:
          site_t=False
      if site_t==False:
        break
#  print(Plat)
#  sys.exit()
  Pos_np=np.array(Pos)
#  print(Atom)
#  print(Pos_np)   
  Plat=Plat*alat*BOHR
  Pos_np=Pos_np*alat*BOHR
#  print(Pos_np)   
  with open ("ctrl.xyz",mode="w") as fp:
    fp.write(" "+str(nbas)+"\n")
    fp.write("\n")  
    for i in range(nbas):   
      fp.write(str("{:>4}".format(Atom[i])))
      fp.write(str("{:>11}".format("{:.6f}".format(Pos_np[i][0]))))
      fp.write(str("{:>11}".format("{:.6f}".format(Pos_np[i][1]))))
      fp.write(str("{:>11}".format("{:.6f}".format(Pos_np[i][2]))))
      fp.write("\n")
  with open ("ctrl.xsf", mode="w") as fp:
    fp.write("CRYSTAL\n")
    fp.write("\n")
    fp.write("PRIMVEC\n")
    for i in range(3):
      fp.write(str("{:>6}".format("{:.3f}".format(Plat[i][0])))) 
      fp.write(str("{:>7}".format("{:.3f}".format(Plat[i][1])))) 
      fp.write(str("{:>7}".format("{:.3f}".format(Plat[i][2])))) 
      fp.write("\n")
    fp.write("\n")
    fp.write("PRIMCOORD\n")
    fp.write(str(nbas)+" 1\n")
    for i in range(nbas):
      fp.write(str("{:>4}".format(Atom[i])))
      fp.write(str("{:>12}".format("{:.8f}".format(Pos_np[i][0])))) 
      fp.write(str("{:>12}".format("{:.8f}".format(Pos_np[i][1])))) 
      fp.write(str("{:>12}".format("{:.8f}".format(Pos_np[i][2])))) 
      fp.write("\n")

if __name__=="__main__":
  main()
