#include <math.h>
#define CSL 274.071979 
void rsq2(double a, double b, double e, double* g, int l, int *nc,
	  int ncmax, int ncmin, int *nod, int nr, int nre, double *rofi,
	  double *slonc, double slore, double *v, double *valnc,
	  double valre, int iz){
/* Integrate the scalar relativistic eqn inward from nre to nc
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   ncmax :upper bound to cutoff nc
Ci   ncmin :lower bound to cutoff nc
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   slore :slope of g at rofi(nre)
Ci   v     :spherical potential = v_true
Ci   valre :value of g(nre) at rofi(nre)
Ci   iz    :nuclear charge
Ci   nre   :rofi(nre) = radius from which RSQ2 integrates
Co Outputs:
Co   g     :wave function times r
Co   nc    :radius to which RSQ2 integrates
Co   nod   :number of nodes between nc and nre
Co   slonc :slope of g(nc)
Co   valnc :value of g(nc)
Cr Remarks:
Cr   Integrates inward from nre to nc
Cr   Cutoff nc is chosen at first maximum, but nc.ge.ncmin
Cr                                         and nc.lt.ncmax
C ----------------------------------------------------------------------*/
/* Local variables: */
    int i,ir,irp1,jr,jrp1;
    double af1,af2,af3,ag1,ag2,ag3,b1,b2,d[2][3],det,df1;
    double df2,df3,dg1,dg2,dg3,dr,drdi,ff,fllp1,gg,h83;
    double phi,q,r,r1,r2,r3,r83sq,rpb,u,vb,x,y,zz;

    *nod  = 0;
    zz    = (double) (2*iz);
    fllp1 = l*(l + 1);
    r83sq = 64./9.;
    r1    =  1./9.;
    r2    = -5./9.;
    r3    = 19./9.;
    h83   = -8./3.;

/* --- FIRST POINT ------ */
    r      = rofi[nre];
    drdi   = a*(r+b);
    phi    = (e + zz/r - v[nre])*drdi/CSL;
    u      = drdi*CSL + phi;
    x      = -drdi/r;
    y      = -fllp1*x*x/u + phi;
    g[nre] = valre;
    g[nre+nr] = (slore*drdi + x*valre)/u;
    ag1    = slore*drdi;
    af1    = x*g[nre+nr] - y*g[nre];
    ir     = nre;
    dg3    = ag1;

/* --- RUNGE-KUTTA FOR NEXT THREE POINTS ----- */
    rpb  = b*exp(a*nre);      /* Fortran: rpb  = b*dexp(a*nre - a) */
    q    = 1./sqrt(exp(a));
    for ( i = 0 ; i < 3; i++ ) {
        irp1 = ir;
        jrp1 = ir + nr;
        ir   = ir - 1;
        jr   = ir + nr;

        rpb   = rpb*q;
        drdi  = rpb*a;
        r     = rpb - b;
        gg    = g[irp1]-0.5*ag1;
        ff    = g[jrp1]-0.5*af1;
        vb    = (3.*v[irp1] + 6.*v[ir] - v[ir-1])*0.125;
        phi   = (e + zz/r - vb)*drdi/CSL;
        u     = drdi*CSL + phi;
        x     = -drdi/r;
        y     = -fllp1*x*x/u + phi;
        ag2   = u*ff - x*gg;
        af2   = x*ff - y*gg;
        gg    = g[irp1]-0.5*ag2;
        ff    = g[jrp1]-0.5*af2;
        ag3   = u*ff - x*gg;
        af3   = x*ff - y*gg;

        rpb   = rpb*q;
        drdi  = a*rpb;
        r     = rpb - b;
        phi   = (e + zz/r - v[ir])*drdi/CSL;
        u     =  drdi*CSL + phi;
        x     = -drdi/r;
        y     = -fllp1*x*x/u + phi;
        gg    = g[irp1] - ag3;
        ff    = g[jrp1] - af3;
        g[ir] = g[irp1] - (ag1+ 2.*(ag2 + ag3) + u*ff - x*gg)/6.;
        g[jr] = g[jrp1] - (af1+ 2.*(af2 + af3) + x*ff - y*gg)/6.;
        if (g[ir]*g[irp1] < 0.) *nod = *nod + 1 ;
        ag1   = u*g[jr] - x*g[ir];
        af1   = x*g[jr] - y*g[ir];
        if (ir == ncmin) goto three ;
	d[0][i] = ag1;
	d[1][i] = af1;
    }

/* --- ALL REMAINING POINTS ----- */
    dg1 = d[0][0];
    dg2 = d[0][1];
    dg3 = d[0][2];
    df1 = d[1][0];
    df2 = d[1][1];
    df3 = d[1][2];
/*      do ir = nre-4,0,-1 */
    for ( ir = nre-4 ; ir > -2 ; ir-- ){
        jr   = ir + nr;
        irp1 = ir + 1;
        jrp1 = irp1 + nr;

        r     = rofi[ir];
        drdi  = a*(r + b);
        phi   = (e + zz/r - v[ir])*drdi/CSL;
        u     = drdi*CSL + phi;
        x     = -drdi/r;
        y     = -fllp1*x*x/u + phi;
        det   = r83sq-x*x+u*y;
        b1    = g[irp1]*h83 + r1*dg1 + r2*dg2 + r3*dg3;
        b2    = g[jrp1]*h83 + r1*df1 + r2*df2 + r3*df3;
        g[ir] = (b1*(h83-x) + b2*u)/det;
        g[jr] = (b2*(h83+x) - b1*y)/det;
        if (g[ir]*g[irp1] < 0.) *nod = *nod + 1;
	dg1   = dg2;
        dg2   = dg3;
        dg3   = u*g[jr] - x*g[ir];
        df1   = df2;
        df2   = df3;
        df3   = x*g[jr] - y*g[ir];
	if ( ir%2 == 0 && ir <= ncmax && 
	    (ir <= ncmin || g[ir]*dg3 >= 0.) ) goto three ;
/*  Fortran->C  changed ir=odd -> ir=even since ir, nr, nre are shifted by -1 
       if (mod(ir,2).ne.0.and.ir.le.ncmax.and.
      .     (ir.le.ncmin.or.g(ir)*dg3.ge.0.d0)) goto 3 */
    }

/*  ----Integration done */
 three:
    *nc    = ir;
    *valnc = g[ir];
    drdi   = a*(rofi[ir] + b);
    *slonc = dg3/drdi;
    return ;
}
