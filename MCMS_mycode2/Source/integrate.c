/***************************************************************************
                          integrate.cc  -  description
                             -------------------
    begin                : Sat Jan 1 2000
    copyright            : (C) 2000 by Alessandro MIRONE
    email                : mirone@lure.u-psud.fr
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
void integrate(double *f, double *fint, double *x, int Nx )
{
  int i;
  double sfin,sgro,d,h,p,q,r,fp,fq,fr;
  
  if(Nx<3)
    {
      printf(" Errore: Nx<3 (Nx = %d) nella routine integrate \n", Nx) ;
      exit(0);
    }
  sfin=0.0;
  sgro=0.0;
  fint[0]=0.0;
  for( i=2;i< (Nx) ;i+=2)
    {
      sgro=sgro+ (f[i-2]+ f[i])*(x[i]-x[i-2]);       

      sfin=sfin+(f[i-2] + f[i-1])*(x[i-1]-x[i-2]);
      sfin=sfin+(f[i-1]+f[i])*(x[i]-x[i-1]);

      fint[i]=(4*sfin-sgro)/6.;
    }  
  
  for( i=1;i< (Nx) ;i+=2)
    {
      if(i<Nx-1)
	{
	  d=x[i]-x[i-1];
	  h=x[i+1]-x[i];

	  p=d*(2*d+3*h)/6./(d+h);
	  q=d*(d+3*h)/6./h;
	  r=-d*d*d/6./h/(d+h);

	  fp=f[i-1];
	  fq=f[i];
	  fr=f[i+1];

	  fint[i]=fint[i-1]+p*fp+q*fq+r*fr;
	}
      else
	{
	  h=x[i-1]-x[i-2];
	  d=x[i]-x[i-1];
	  p=-d*d*d/6./h/(d+h);
	  q=d*(d+3*h)/6./h;
	  r=d*(2*d+3*h)/6./(d+h);

	  fp=f[i-2];
	  fq=f[i-1];
	  fr=f[i];

	  fint[i]=fint[i-1]+p*fp+q*fq+r*fr;
	}
    }
}


double  defintegr(double *f,  double *x, int Nx )
{
  int i;
  double sfin,sgro,sum,d,h,p,q,r,fp,fq,fr;
  
  if(Nx<3)
    {
      printf(" Errore: Nx<3 (Nx = %d) nella routine defintegr \n", Nx) ;
      exit(0);
    }
  sfin=0.0;
  sgro=0.0;
  sum=0.0; 

  for( i=2;i< (Nx) ;i+=2)
    {
      sgro=sgro+ (f[i-2]+ f[i])*(x[i]-x[i-2]);

      sfin=sfin+(f[i-2] + f[i-1])*(x[i-1]-x[i-2]);
      sfin=sfin+(f[i-1]+f[i])*(x[i]-x[i-1]);

      sum=(4*sfin-sgro)/6.;
    }  

  {
    i=((Nx-1)>>1)*2;

    if(i!=(Nx-1))
      {
        i++;
	h=x[i-1]-x[i-2];
	d=x[i]-x[i-1];
	p=-d*d*d/6./h/(d+h);
	q=d*(d+3*h)/6./h;
	r=d*(2*d+3*h)/6./(d+h);
	
	fp=f[i-2];
	fq=f[i-1];
	fr=f[i];
	
	sum=sum + p*fp+q*fq+r*fr;
      }
  }
  return sum;

}
