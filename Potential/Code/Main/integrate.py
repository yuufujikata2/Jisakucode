def integrate(f,x,Nx):
    sfin=0.0
    sgro=0.0
    sum_=0.0

    for i in range(2,Nx,2):
        sgro=sgro+ (f[i-2]+ f[i])*(x[i]-x[i-2])
  
        sfin=sfin+(f[i-2] + f[i-1])*(x[i-1]-x[i-2])
        sfin=sfin+(f[i-1]+f[i])*(x[i]-x[i-1])
  
        sum_=(4*sfin-sgro)/6.
      
    i=((Nx-1)>>1)*2


    if i != (Nx-1):
        i += 1
        h=x[i-1]-x[i-2]
        d=x[i]-x[i-1]
        p=-d*d*d/6./h/(d+h)
        q=d*(d+3*h)/6./h
        r=d*(2*d+3*h)/6./(d+h)

        fp=f[i-2]
        fq=f[i-1]
        fr=f[i]

        sum_=sum_ + p*fp+q*fq+r*fr

    return sum_

