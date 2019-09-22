import numpy as np
import sys

def main():
  n=input('Please input dimension: ').split()
  if len(n)==0:
    print('error:no dimension')
    sys.exit()
  elif len(n)==1:
    n1=1
    n2=int(n[0])
    A=np.zeros(n2)
  elif len(n)==2:
    n1=int(n[0])
    n2=int(n[1])
    A=np.zeros((n1,n2))
  elif len(n)>=3:
    print('error:It is not matrix')
    sys.exit()

  for i in range(n1):
    while True:
      ni=input('please input matrix raw'+str(i)+':').split()
      if len(ni)!=n2:
        print('error:column is bad')
      else:
        if len(n)==1: 
          for j in range(n2):
            A[j]=float(ni[j])
        else:
          for j in range(n2):
            A[i][j]=float(ni[j])
        break  
  print(A)
  A_inv=np.linalg.inv(A)
  print(A_inv)
  m=input('Please input dimension2: ').split()
  if len(m)==0:
    print('error:no dimension')
    sys.exit()
  elif len(m)==1:
    m1=1
    m2=int(m[0])
    B=np.zeros(m2)
  elif len(m)==2:
    m1=int(m[0])
    m2=int(m[1])
    B=np.zeros((m1,m2))
  elif len(m)>=3:
    print('error:It is not matrix')
    sys.exit()

  for i in range(m1):
    while True:
      mi=input('please input matrix2 raw'+str(i)+':').split()
      if len(mi)!=m2:
        print('error:column is bad')
      else:
        if len(m)==1: 
          for j in range(m2):
            B[j]=float(mi[j])
        else:
          for j in range(m2):
            B[i][j]=float(mi[j])
        break  

  if len(m)==1:
    B=B.T
    m1=m2
  print(B)
  ca=input('calculation? y or n:')
  if ca!='y':
    sys.exit()
  if n2!=m1:
    print('error: column of matrix1 is not raw of matrix2')
    sys.exit()
  C=np.dot(A_inv,B)
  print('result')
  print(C)

if __name__=='__main__':
  main()
