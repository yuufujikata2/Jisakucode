"""
To chenge oneline times BOHR or etc

"""

path='indata'
path_w='outdata'
f=open(path)
with open(path_w,mode='w') as f2:
  for i in f:
    f2.write(str('{:.8f}'.format(float(i)*13.6  ))+'\n')
f.close()
    
