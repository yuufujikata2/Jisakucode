"""
To create FPMS format for V2O5 from cluster.xyz file

"""
import sys

args=sys.argv
argc=len(args)
shape="NM"
if argc==1:
  print ("error : please input xyzfile!")
  sys.exit()
elif argc==3:
  shape=args[2]
elif argc>3:
  print("error : ")
fw="clus.fp"
fw2="clus2.fp"

r_V="1.47874454"
r_O="1.29733269"
r_O1="0.98117380"
r_O2="1.45194048"
r_E="2.68331799"
r_E1="1.37471562"
r_E2="0.96616869"
r_E3="1.00534521"
r_E4="0.77305703"
r_E5="0.69294226"
r_E6="0.65198121"


f=open(args[1])

number=0
gyou=0
if shape=="NM":
  with open(fw,mode='w') as fw:
    for i in f:
      gyou+=1
      if gyou<=2:
        continue
      line=i.split()
      number+=1
      if line[0]=="V":
        fw.write("  "+str(23)+str('{:>12}'.format(line[1]))+str('{:>11}'.format(line[2]))+str('{:>11}'.format(line[3]))+str('{:>12}'.format(r_V))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("V"))+"\n")
      elif line[0]=="O":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format(line[3]))+str('{:>12}'.format(r_O))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O"))+"\n")
      elif line[0]=="O1":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_O1))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O1"))+"\n")
      elif line[0]=="O2":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_O2))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O2"))+"\n")
  
      elif line[0]=="E":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E"))+"\n")
      elif line[0]=="E1":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E1))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E1"))+"\n")
      elif line[0]=="E2":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E2))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E2"))+"\n")
      elif line[0]=="E3":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E3))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E3"))+"\n")
      elif line[0]=="E4":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E4))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E4"))+"\n")
      elif line[0]=="E5":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E5))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E5"))+"\n")
      elif line[0]=="E6":
        fw.write("   "+str(0)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>12}'.format(r_E6))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("E6"))+"\n")
  
  f.close()
  
  with open(fw2,mode="w") as fw2:
    fw=open("clus.fp")
    array_V=[]
    array_O=[]
    array_O1=[]
    array_O2=[]
    array_E=[]
    array_E1=[]
    array_E2=[]
    array_E3=[]
    array_E4=[]
    array_E5=[]
    array_E6=[]
    
    for fw in fw:
      line2=fw.split()
      if line2[7]=="V":
        array_V.append(fw)
      elif line2[7]=="O":
        array_O.append(fw)
      elif line2[7]=="O1":
        array_O1.append(fw)
      elif line2[7]=="O2":
        array_O2.append(fw)
      elif line2[7]=="E":
        array_E.append(fw)
      elif line2[7]=="E1":
        array_E1.append(fw)
      elif line2[7]=="E2":
        array_E2.append(fw)
      elif line2[7]=="E3":
        array_E3.append(fw)
      elif line2[7]=="E4":
        array_E4.append(fw)
      elif line2[7]=="E5":
        array_E5.append(fw)
      elif line2[7]=="E6":
        array_E6.append(fw)
  
    for i in array_V:
      fw2.write(i)
    for i in array_O:
      fw2.write(i) 
    for i in array_O1:
      fw2.write(i) 
    for i in array_O2:
      fw2.write(i) 
    for i in array_E1:
      fw2.write(i) 
    for i in array_E2:
      fw2.write(i) 
    for i in array_E3:
      fw2.write(i) 
    for i in array_E4:
      fw2.write(i) 
    for i in array_E5:
      fw2.write(i) 
    for i in array_E6:
      fw2.write(i) 
    for i in array_E:
      fw2.write(i) 
    
elif shape=="MT":
  with open(fw,mode='w') as fw:
    for i in f:
      gyou+=1
      if gyou<=2:
        continue
      line=i.split()
      number+=1
      if line[0]=="V":
        fw.write("  "+str(23)+str('{:>12}'.format(line[1]))+str('{:>11}'.format(line[2]))+str('{:>11}'.format(line[3]))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("V"))+"\n")
      elif line[0]=="O":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format(line[3]))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O"))+"\n")
      elif line[0]=="O1":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O1"))+"\n")
      elif line[0]=="O2":
        fw.write("   "+str(8)+str('{:>12}'.format((line[1])))+str('{:>11}'.format((line[2])))+str('{:>11}'.format((line[3])))+str('{:>7}'.format(str(number)))+str('{:>10}'.format("1.00000"))+str('{:>4}'.format("O2"))+"\n")
  
  
  f.close()
  
  with open(fw2,mode="w") as fw2:
    fw=open("clus.fp")
    array_V=[]
    array_O=[]
    array_O1=[]
    array_O2=[]
    
    for fw in fw:
      line2=fw.split()
      if line2[6]=="V":
        array_V.append(fw)
      elif line2[6]=="O":
        array_O.append(fw)
      elif line2[6]=="O1":
        array_O1.append(fw)
      elif line2[6]=="O2":
        array_O2.append(fw)
  
    for i in array_V:
      fw2.write(i)
    for i in array_O:
      fw2.write(i) 
    for i in array_O1:
      fw2.write(i) 
    for i in array_O2:
      fw2.write(i) 

