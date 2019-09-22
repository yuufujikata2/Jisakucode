#include <math.h>
#define CSL 274.071979 
void rsq1(double a, double b, double e, double* g, int l, int nc,
	  int* nod, int nr, double* rofi, double* slonc, double* v, 
	  double* valnc, int iz) {
/*
C- integrate the scalar relativistic eqn outward from 0 to rofi(nc)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a*i) -1] 
          i=0,...,nr-1  before in Fortran :   b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   nc    :integration from 0.0 to rofi(nc)
Ci   nr    :number of mesh points   
Ci   rofi  :radial mesh points
Ci   v     :spherical potential = v_true
Ci   iz    :nuclear charge
Co Outputs:
Co   g     :wave function times r
Co   nod   :number of nodes
Co   slonc :slope of g at rofi(nc)
Co   valnc :value of g at rofi(nc)
Cr Remarks:
Cr   Boundary condition does not fix value of wave function near the
Cr   origin; integration can be scaled by an arbritrary factor.
Cr   Scalar relativistic version
C ----------------------------------------------------------------------
*/
/*  Local variables: */
    int ir,irm1,jr,jrm1;
    double aa,b1,b2,d[2][3],det,df1,df2,df3,dg1,dg2,dg3;
    double drdi,f0,fllp1,g0,h83,phi,r,r1,r2,r3,r83sq,s,sf,u,x,x2,y,zz;
    *nod  = 0;
    zz    = (double) (2*iz);
    fllp1 = l*(l + 1);
    r83sq = 64./9.;
    r1    =  1./9.;
    r2    = -5./9.;
    r3    = 19./9.;
    h83   =  8./3.;

/* --- Approximate g,f by leading term near zero */
    g0   = 1.;
      if (iz<0.9) {
	  s  = l+1;
	  sf = l;
	  f0 = l/CSL;
      }
      else {
	  aa = zz/CSL;
	  s  = sqrt(fllp1 + 1. - aa*aa);
	  sf = s;
	  f0 = g0*(s - 1.)/aa;
      }
      *g    = 0.;
      *(g+nr) = 0.;
      for (ir = 1; ir < 4; ir++){
	  jr    = ir + nr;
	  r     = *(rofi+ir);
	  drdi  = a*(r + b);
	  *(g+ir) = g0*pow(r,s);
	  *(g+jr) = f0*pow(r,sf);
	  d[0][ir-1] = drdi* *(g+ir) * s/r;
	  d[1][ir-1] = drdi* *(g+jr) * sf/r;
      }
/* --- integrate over rest of points ------ */
    dg1 = d[0][0];
    dg2 = d[0][1];
    dg3 = d[0][2];
    df1 = d[1][0];
    df2 = d[1][1];
    df3 = d[1][2];
    for (ir = 4; ir <= nc; ir++){
        jr    = ir + nr;
        irm1  = ir - 1;
        jrm1  = irm1 + nr;
        r     = *(rofi+ir);
        drdi  = a*(r + b);
        phi   = (e + zz/r - *(v+ir))*drdi/CSL;
        u     = drdi*CSL + phi;
        x     = -drdi/r;
        x2    = x*x;
        y     = -fllp1*x2/u + phi;
        det   = r83sq-x2+u*y;
        b1    = *(g+irm1) * h83 + r1*dg1 + r2*dg2 + r3*dg3;
        b2    = *(g+jrm1) * h83 + r1*df1 + r2*df2 + r3*df3;
        *(g+ir) = (b1*(h83-x) + b2*u)/det;
        *(g+jr) = (b2*(h83+x) - b1*y)/det;
        if  (*(g+ir) * *(g+irm1) < 0.) (*nod)++ ;
        dg1   = dg2;
        dg2   = dg3;
        dg3   = u * *(g+jr) - x * *(g+ir);
        df1   = df2;
        df2   = df3;
        df3   = x * *(g+jr) - y * *(g+ir);
    }
    *valnc = *(g+nc);
    drdi  = a*(*(rofi+nc) + b);
    *slonc = dg3/drdi;
}
