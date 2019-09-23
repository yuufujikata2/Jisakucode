#include <math.h>
void gintsr(double a, double b, double* ga, double* gb, double* gfac,
	    int nr, double* rofi, double* sum);
void norm1g(double a, double b, double* g, double* gfac, int nr, double* rofi, double* slo, double* val){
/*
C- Integrate inner product of a wave function \int g*g dr + rel. corrections
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Cio  g     :wave function (times r)
Ci   gfac  :normalization factor for first component of g
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Cio  slo   :Derivative of wave function at sphere radius
Cio  val   :value of wave function at sphere radius
Co Outputs: for normalized wave function
Cio  g     :wave function (times r)
Cio  slo   :Derivative of wave function at sphere radius
Cio  val   :value of wave function at sphere radius
Cr Remarks:
Cr   fac   :factor to normalize g
C ----------------------------------------------------------------------
*/
/*  Local variables: */
    int ir,jr;
    double fac,sum;

    gintsr(a,b,g,g,gfac,nr,rofi,&sum);
    fac = 1./sqrt(sum);
/*
      if (g[1]<0.) fac = -fac;      
        fac = dsign(1.d0,g(2))/sqrt(sum)
*/
    for (ir=0; ir<2*nr; ir++) g[ir] *= fac; 
    *slo *= fac;
    *val *= fac;
}
