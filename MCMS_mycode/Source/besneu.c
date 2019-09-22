#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define SMALLARG 1.e-2
void bessel( int l, double x, double *pjl, double *pjlp ) ;
void neuman( int l, double x, double *pnl, double *pnlp ) ;
void besselmod( int l, double x, double *pil, double *pilp ) ;
void hankelmod( int l, double x, double *phl, double *phlp ) ;
/*
  besneu() computes spherical Bessel (j_l) and Neuman (n_l) functions
  or modified versions (i_l and -m_l) - ATTENTION SIGN nl == -m_l 
  value of "modified"    *pjl   *pnl
               0          j_l    n_l
               1          i_l   -m_l 
  as well as their first derivatives *pnlp == d/dx n_l or -d/dx m_l, etc

      ATTENTION FOR NON STANDARD SIGN OF MODIFIED NEUMAN FUNCTION
*/
void besneu( int modified, double x, int l, double* pjl, double* pnl, double* pjlp, double* pnlp) 
{
  if ( modified ) {
    double hl, hlp ;
    besselmod( l, x, pjl, pjlp ) ;
    hankelmod( l, x, &hl, &hlp ) ;
    *pnl  = hl - *pjl ;
    *pnlp = hlp - *pjlp ;
  } 
  else {
    bessel( l, x, pjl, pjlp ) ;
    neuman( l, x, pnl, pnlp ) ;
  }
}   

double dfactfact( int n ) 
{
  double result = 1. ;
  while ( n > 0) {
    result *= n ;
    n -= 2 ;
  }
  return result ;
}
void bessel( int l, double x, double *pjl, double *pjlp ) 
{
  if ( l < 0 || x < 0. ) {
    printf("l < 0 || x < 0. in bessel\n") ; exit(1) ; }
  if ( x < SMALLARG ) {       /* series expansion up to x^4 [A.S. 10.1.2] */
    double leading, halfx2, expansion ;
    leading = pow(x,l) / dfactfact( 2*l + 1 ) ;
    halfx2 = 0.5*x*x ;
    expansion = 1. - halfx2/(2*l+3) + halfx2*halfx2/(2*(2*l+3)*(2*l+5)) ;
    *pjl = leading*expansion ;
    if ( l == 0 ) {
      leading = -x/3. ;
      expansion = 1 - halfx2/5. + halfx2*halfx2/70. ; 
    } 
    else {
      leading = pow(x,l-1) / dfactfact( 2*l + 1 ) ;
      expansion = l - (l+2)*halfx2/(2*l+3) + (l+4)*halfx2*halfx2/(2*(2*l+3)*(2*l+5)) ;
    }
    *pjlp = leading*expansion ;
  }
  else { /* downward recursion with Miller's method [A.S. 10.1.19,22; 10.5] */
    int istart, i ;
    double Fi, Fiplus1, Fiplus2, Fl, dFl, FloverF0, dFloverF0, j0 ;
    istart = l + (int) x + 10 ;
    Fiplus2 = 0. ;
    Fiplus1 = 1. ;
    for ( i = istart ; i >= 0 ; i-- ) {
      Fi = (2*i+3)*Fiplus1/x -Fiplus2 ;
      if ( i == l ) { 
	Fl = Fi ;
	dFl = i*Fi/x - Fiplus1 ;
      }
      Fiplus2 = Fiplus1 ;
      Fiplus1 = Fi ;
    }
    FloverF0 = Fl/Fi ;
    dFloverF0 = dFl/Fi ;
    j0 = sin(x)/x ;
    *pjl = FloverF0 * j0 ;
    *pjlp = dFloverF0 * j0 ;
  }
}
void neuman( int l, double x, double *pnl, double *pnlp ) 
{
  if ( l < 0 || x <= 0. ) {
    printf("l < 0 || x <= 0. in neuman\n") ; exit(1) ; }
  if ( x < SMALLARG ) {       /* series expansion up to x^4 [A.S. 10.1.3] */
    double halfx2, leading, expansion ;
    halfx2 = 0.5*x*x ;
    leading = -dfactfact( 2*l - 1 ) / pow(x,l+1) ; 
    expansion = 1. - halfx2/(1-2*l) + halfx2*halfx2/(2*(1-2*l)*(3-2*l)) ;
    *pnl = leading*expansion ;
    leading =  dfactfact( 2*l - 1 ) / pow(x,l+2) ; 
    expansion = l+1 - (l-1)*halfx2/(1-2*l) + (l-3)*halfx2*halfx2/(2*(1-2*l)*(3-2*l)) ;
    *pnlp = leading*expansion ;
  }
  else {   /* upward recurrence relation [A.S. 10.1.19 and 10.1.22] */
    int i ;
    double ni, niplus1, niminus1 ;
    ni = -cos(x)/x ;
    niplus1 = -cos(x)/x/x - sin(x)/x ;
    for ( i = 1 ; i <= l ; i++ ) {
      niminus1 = ni ;
      ni = niplus1 ;
      niplus1 = (2*i+1)*ni/x - niminus1 ;
    }
    *pnl = ni ;
    *pnlp = l*ni/x - niplus1 ;
  }
}
/* 
   Modified spherical Bessel functions of the first kind [A.S. 10.2.2]
   i_l(x) = sqrt(pi/2/x) I_{l+1/2}(x)    
   output: *pil == i_l, *pilp == (d/dx)i_l 
*/
void besselmod( int l, double x, double *pil, double *pilp ) 
{
  if ( l < 0 || x < 0. ) {
    printf("l < 0 || x < 0. in bessel\n") ; exit(1) ; }
  if ( x < SMALLARG ) {       /* series expansion up to x^4 [A.S. 10.2.5] */
    double leading, halfx2, expansion ;
    leading = pow(x,l) / dfactfact( 2*l + 1 ) ;
    halfx2 = 0.5*x*x ;
    expansion = 1. + halfx2/(2*l+3) + halfx2*halfx2/(2*(2*l+3)*(2*l+5)) ;
    *pil = leading*expansion ;
    if ( l == 0 ) {
      leading = x/3. ;
      expansion = 1 + halfx2/5. + halfx2*halfx2/70. ; 
    } 
    else {
      leading = pow(x,l-1) / dfactfact( 2*l + 1 ) ;
      expansion = l + (l+2)*halfx2/(2*l+3) + (l+4)*halfx2*halfx2/(2*(2*l+3)*(2*l+5)) ;
    }
    *pilp = leading*expansion ;
  }
  else { /* downward recursion with Miller's method [A.S. 10.2.18,21; 10.5] */
    int istart, i ;
    double Fi, Fiplus1, Fiplus2, Fl, dFl, FloverF0, dFloverF0, j0 ;
    istart = l + (int) x + 10 ;
    Fiplus2 = 0. ;
    Fiplus1 = 1. ;
    for ( i = istart ; i >= 0 ; i-- ) {
      Fi = (2*i+3)*Fiplus1/x + Fiplus2 ;
      if ( i == l ) { 
	Fl = Fi ;
	dFl = i*Fi/x + Fiplus1 ;
      }
      Fiplus2 = Fiplus1 ;
      Fiplus1 = Fi ;
    }
    FloverF0 = Fl/Fi ;
    dFloverF0 = dFl/Fi ;
    j0 = 0.5*(exp(x)-exp(-x))/x ;
    *pil = FloverF0 * j0 ;
    *pilp = dFloverF0 * j0 ;
  }
}
/*       modified spherical Bessel functions of the first kind (i_l)
   MINUS modified spherical Bessel functions of the second kind (m_l)
   h_l(x) = i_l(x) - m_l(x) = sqrt(pi/2/x) ( I_{l+1/2}(x) - I_{-l-1/2}(x) )   
   Also: h_l = 2/PI * (-)^{n+1} * k_l, where k_l(x) = sqrt(pi/2/x) K_{l+1/2}(x)
   output: *phl == h_l, *phlp == (d/dx)h_l                  [A.S. 10.2.2,3,4]
*/
void hankelmod( int l, double x, double *phl, double *phlp ) 
{
  if ( l < 0 || x <= 0. ) {
    printf("l < 0 || x <= 0. in hankel\n") ; exit(1) ; }
  if ( x < SMALLARG ) {       /* series expansion up to x^5 [A.S. 10.2.6] */
    double halfx2, leading, expansion ;
    if ( l == 0 ) {
      *phl  = -1./x + 1. - x/2. + x*x/6. - x*x*x/24. + x*x*x*x/120. ;
      *phlp =  1./x/x    - 1/2. + x/3.   - x*x/8.    + x*x*x/30. ;
    } 
    else if ( l == 1 ) {
      *phl  =  1./x/x    - 1./2. + x/3.   - x*x/8. + x*x*x/30. ;
      *phlp =  -2./x/x/x         + 1/3.   - x/4.   + x*x/10. ;
    }
    else if ( l == 2 ) {
      *phl  =  -3./x/x/x    + 0.5/x    -  x/8.   +  x*x/15. ;
      *phlp =   9./x/x/x/x  - 0.5/x/x  - 1./8.   + 2.*x/15. ;
    }
    else {
      halfx2 = 0.5*x*x ;
      leading = (l%2 ? 1:-1) * dfactfact( 2*l - 1 ) / pow(x,l+1) ; 
      expansion = 1. + halfx2/(1-2*l) + halfx2*halfx2/(2*(1-2*l)*(3-2*l)) ;
      *phl = leading*expansion ;
      leading = (l%2 ? 1:-1) * dfactfact( 2*l - 1 ) * (-1.) / pow(x,l+2) ; 
      expansion = l+1 + (l-1)*halfx2/(1-2*l) + (l-3)*halfx2*halfx2/(2*(1-2*l)*(3-2*l)) ;
      *phlp = leading*expansion ;
    }
  }
  else {   /* upward recurrence relation [A.S. 10.2.18 and 10.1.21] */
    int i ;
    double hi, hiplus1, himinus1 ;
    hi = -exp(-x)/x ;
    hiplus1 = -hi*(1.+1./x) ;
    for ( i = 1 ; i <= l ; i++ ) {
      himinus1 = hi ;
      hi = hiplus1 ;
      hiplus1 = -(2*i+1)*hi/x + himinus1 ;
    }
    *phl = hi ;
    *phlp = l*hi/x + hiplus1 ;
  }
}
#undef SMALLARG
