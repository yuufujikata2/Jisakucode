#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NRMAX 1000
#define NITMAX 40
#define TOLRSQ 1.e-12
/*       parameter (nrmax=1000,nitmax=40,tolrsq=1.d-12) */
void fctp0(int l, int *nctp0, int nr, double *rofi, double *v, int iz) ;
void fctp(double a, double b, double e, int l, int *nctp, int nctp0,
	  int nr, double *rofi,double *v, int iz) ;
void mkgfac(double e, double* gfac, int l, int nr, double* rofi, double* v, int iz);
void rsq1(double a, double b, double e, double* g, int l, int nc,
          int* nod, int nr, double* rofi, double* slonc, double* v,
          double* valnc, int iz);
void rsq2(double a, double b, double e, double* g, int l, int *nc,
	  int ncmax, int ncmin, int *nod, int nr, int nre, double *rofi,
	  double *slonc, double slore, double *v, double *valnc,
	  double valre, int iz);
void norm1g(double a, double b, double* g, double* gfac, int nr, double* rofi, double* slo, double* val);
void rseq(double a, double b, double *e, double eb1, double eb2,
	  double *g, double *gfac, int l, int nod, int nr, int *nre,
	  double *rofi, double *slo, double *v, double *val, int iz) {
/* Solves radial Schroedinger equation for given BCs and number of nodes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   clabl :name of the different inequivalent atoms
Ci   eb1   :lower bound for energy
Ci   eb2   :upper bound for energy
Ci   l     :angular momentum
Ci   nod   :number of nodes
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential = v_true
Ci   z     :nuclear charge
Ci Inputs/Outputs:
Cio  e     :energy
Cio  slo   :slope of g(nre)
Cio  val   :value of g(nre)
Co Outputs:
Co   g     :normalized wave function times r
Co   gfac  :normalization factor for first component of g
Co   nre   :rofi(nre)=radius to which rseq integrates
Cr Remarks:
Cr   tolrsq:maximum relative error in energy; absolute error if |e| < 1
Cr   Nonrelativistic calculation is obtained by replacing the
Cr   light velocity by a very large number (1.d10) in ATSCPP.
Cr
Cr   Output wavefunction normalized to 1: so g(nr).ne. val*wsr
Cr
Cr   Note: if r=b(exp(a*z) -1] then dr=A(r + B) dz
C ---------------------------------------------------------------------- */
/*      double g(2*nr),gfac(nr), rofi(nr),v(nr) */
    int ir,nc,ncmax,ncmin,nctp,nctp0,ierr,nit,nod1,nod2,node,i;
    double ein,de,dphi,drdi[NRMAX],e1,e2,sum,fi,phi,r,ratio,re,slonc1;
    double slonc2,slore,tsrme,valnc1,valnc2,valre;
/*      character*144 messg
      character*6 name
      character*1 cl(0:9)
C External calls:
      external dinit,dscal,errmsg,fctp0,fctp,iprint,lunit,norm1g,rsq1,
     .         rsq2
C Intrinsic functions:
      intrinsic dabs,dlog,dmax1,dsqrt,idnint,max0,min0
C Data statements:
      data cl/'s','p','d','f','g','h','i','j','k','l'/
*/
    
/*  if (nod < 0 || nod > 7 ) { */
    if (nod < 0 ) {	
	printf(" RSEQ  :bad number of nodes\n");
	exit(1);
    }
    fctp0(l,&nctp0,nr,rofi,v,iz);
    nctp=(nctp0+nr-1)/2;

    for (ir = 1; ir < nr; ir++)
        drdi[ir] = a*(rofi[ir]+b);
    ein = *e;
    
/* --- START ITERATIONS TO FIND ENERGY ------ */
    ierr = 0;
one:     
    e1   = eb1;
    e2   = eb2;
    *e    = ein;
    for ( nit = 1 ; nit <= NITMAX ; nit++ ) {
        if ( *e <= e1 || *e >= e2 ) *e = (e1+e2)*0.5;

/* ----- get nre */
        fctp(a,b,*e,l,&nctp,nctp0,nr,rofi,v,iz);
        re  = 15. * rofi[nctp];
        *nre = log(re/b + 1.)/a ;
        *nre = ((*nre+1)/2)*2 ;
        *nre = *nre > 34 ? *nre : 34 ;     /* max0(nre,35) */
        *nre = *nre < nr-1 ? *nre : nr-1 ; /* min(nre,nr)  */
        if (ierr != 0) *nre = nr-1;

/* ----- get valre,slore */
        if (*nre == nr-1) {
	    valre=*val;
	    slore=*slo;
	}
        else {
	    valre =  1.e-30;
	    slore = -1.e-30;
	}

/* ----- get ncmin */
        if (nod != 0) 
	    ncmin = 30 ;
        else if (valre*slore > 0.)
	    ncmin = *nre-10 ;
        else
	    ncmin = *nre/3 ;

/* ----- get ncmax */
        if (ierr == 0 || nit < 5) ncmax = *nre ;

/* ----- Integrate the scalar relativistic equation inward from rofi(nre)
   ----- to rofi(nc) and then outward from rofi(0) to rofi(nc)
*/
        rsq2(a,b,*e,g,l,&nc,ncmax,ncmin,&nod2,nr,*nre,rofi,&slonc2,
		  slore,v,&valnc2,valre,iz);
        rsq1(a,b,*e,g,l,nc,&nod1,nr,rofi,&slonc1,v,&valnc1,iz);
        node  = nod1 + nod2;
        ratio = valnc2/valnc1;
	
        if (node != nod) {
	    if (nit >= NITMAX-5) 
		printf("WARNING: node != nod && nit >= NITMAX-5 in rseq\n");
	    if (node > nod) e2 = *e ;
	    if (node < nod) e1 = *e ;
	    *e = ( e1 + e2 ) * 0.5 ;
	}
        else { 
/* ------  Calculate sum=norm of wave function from trapezoidal rule */
	    sum = 0.;
	    for ( ir = 1 ; ir <= nc ; ir++ )         /* do ir = 2, nc */
		sum += drdi[ir]*g[ir]*g[ir] ;
	    sum *= ratio*ratio ;
	    for ( ir = nc + 1 ; ir <= *nre ; ir++ )  /* do ir = nc + 1, nre */
		sum += drdi[ir]*g[ir]*g[ir] ;
	    sum -= drdi[*nre]*g[*nre]*g[*nre]*0.5 ;   /* nre kept */

/* ------- DE=estimated correction to eigenvalue */
	    de = -valnc2 * (slonc2 - ratio*slonc1)/sum ;
	    if ( de > 0. ) 
		e1 = *e ;
	    else
		e2 = *e ;
	    *e = *e + de ;
	    if (fabs(de/(fabs(*e)>1.?fabs(*e):1.)) < TOLRSQ) goto two ;
        }
    }
/* --- END OF ITERATION TO FIND ENERGY ------
   --- Search for eigenvalue failed   */
    printf("Search for eigenvalue failed\n");
    if (nod != node) 
	printf("nod != node: NITMAX=%d, nod=%d, node=%d, e=%lf, de=%lf\n",
		 NITMAX,nod,node,*e,de);
    else
	printf("NITMAX=%d, e=%lf, de=%lf\n", NITMAX, *e, de);
    ierr++;
/*    call errmsg(messg,ierr) */
/* c     if (dabs(de/dmax1(dabs(e),1.d0)) .lt. TOLRSQ*10000**ierr) goto 2*/
    goto one ;

two:
/* --- Normalize G */
    *e = *e - de;
    for (i=0;i<=nc;i++) { g[i] *= ratio ; g[nr+i] *= ratio ; }
    if (*nre < nr-1) for (i=*nre+1;i<nr;i++) { g[i] = 0.; g[nr+i] = 0.; }
    
    mkgfac(*e,gfac,l,nr,rofi,v,iz);

    *val  = valre;
    *slo  = slore;
    norm1g(a,b,g,gfac,nr,rofi,slo,val);
    phi  = *val/rofi[nr-1];
    dphi = (*slo-phi)/rofi[nr-1];
}

/*
      if (iprint().ge.120) write(lunit(1),305)
      if (iprint().ge.100) then
        write(lunit(1),303)name,e,sum,nr,nre,nc,nit
        write(lunit(1),304)name,phi,dphi
      endif

299   format(a4,i1,a1)
300   format(/' RSEQ  : ',a6,'  BC=',d8.2,d9.2,'  E1,E2=',f8.1,f6.1)
301   format(' NIT L NODE NOD NRE NC',9x,'E1',14x,'E',13x,'E2',11x,'DE',
     .       /80('-'))
302   format(2i3,1x,i1,'+',i1,3i4,3f15.7,1p,d13.4)
303   format(' RSEQ  : ',a6,'  E=',f12.5,'  SUM=',1pd9.2,
     .       '  NR/NRE/NC=',3i4,'  NIT=',i2)
304   format(' RSEQ  : ',a6,'  PHI ,DPHI= ',2f11.6)
305   format(80('-'))

400   format(' RSEQ  : ',a6,'  NIT >',i3,' and bad nodes,  NOD=',i1,
     .       ' NODE=',i2,
     .       '|                E=',1pd12.5,' DE=',1pd12.5,'$')
401   format(' RSEQ  : ',a6,'  NIT >',i3,' and E+DE outside [E1,E2]',
     .       '|                E=',1pd12.5,' DE=',1pd12.5,'$')
402   format(14x,'-------> RETRY')
      end

*/




