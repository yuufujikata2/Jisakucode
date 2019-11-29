import numpy as np

x,y=np.polynomial.legendre.leggauss(7)
print(x*2/np.sqrt(3))
print(y)
