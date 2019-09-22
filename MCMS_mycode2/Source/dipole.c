#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "w3j.h"
#include "globals.h"

#define EPSILON 1.e-10

double defintegr(double *f,  double *x, int Nx ) ;


/* calculate \int_0^r_0 P_inish(r) r P_finsh(r) dr */
double radipmatele( int inish, int finsh, int nr, 
		    double *rofi, double **ppofi ) {
  int ir ;
  double *f, result ;

  f = ( double * ) malloc( nr * sizeof( double ) ) ;
  
  for ( ir = 0 ; ir < nr ; ir++ ) {
    f[ir] = ( ppofi[inish][ir] * ppofi[finsh][ir] 
	      + ppofi[inish][nr+ir] * ppofi[finsh][nr+ir] ) * rofi[ir] ;
  }
  
  result =  defintegr(f,rofi,nr) ;

  printf("In radipmatele: inish = %d, finsh = %d, radipmatele = %lf\n",
	 inish, finsh, result );
  free (f) ;
  
  return result ;
}


/* 
   dipmatele2list  calculates the matrix elements of the dipole operator 
   r^(1)_q == r C^(1)_q, using formulae (11.15), (14.58) of Cowan app F-7. 
   for given  nli  and q  (where  l,m <=> lf,mf and l',m' <=> li,mi).
   It writes the non-zero elements to a linked list. 
   v(i1,i2) c^+_i2 c_i1, where i1 is a core and i2 a valence orbital.
*/

int dipole2list( int nshells, int *lsh, int *sorb1sh, int inish, int q,
		 int nr, double *rofi, double **ppofi,
                 struct O1plistitem **ppitem0 )
{
  int finsh, li, lf, lgreater, sign, mi, mf, sigma, iorbi, iorbf ;
  int counter = 0 ;
  double reducedmatele, val ;
  
  /* uses formulae (14.58) (11.15) of Cowan app F.7 */
  /* with l,m <=> lf,mf and l',m' <=> li,mi         */

  li = lsh[inish] ;
  for ( finsh = 0 ; finsh < nshells ; finsh++ ) {
    lf = lsh[finsh] ;
    if ( lf == li + 1 || lf == li - 1 ) {
      /* calculate reduced matrix element using (14.58) */
      lgreater = li > lf ? li : lf ;
      sign =  ( lf + lgreater ) % 2 == 0 ? 1 : -1 ;
      reducedmatele = sign * sqrt( lgreater ) * 
	radipmatele( inish, finsh, nr, rofi, ppofi ) ;
      /* calculate matrix element using Wigner Eckart Theorem (11.15) */
      for (  mi = -li ; mi <= li ; mi++ ) {
	mf = mi + q ;
	if ( -lf <= mf && mf <= lf ) {
	  sign = ( lf - mf ) % 2 == 0 ? 1 : -1 ;
	  val = sign * w3j( 2*lf, 2, 2*li, -2*mf, 2*q, 2*mi ) * reducedmatele ;
	printf("li=%d, mi=%d, lf=%d, mf=%d, si=%d, w3j=%lf, re=%lf, val=%lf\n",
	li,mi,lf,mf,sign,w3j(2*lf,2,2*li,-2*mf,2*q,2*mi),reducedmatele,val ) ;
	  if ( fabs(val) > EPSILON ) {
	    /* write matrix elements in list */
	    for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	      iorbi = sporbindex( sorb1sh, lsh, inish, mi, sigma ) ;
	      iorbf = sporbindex( sorb1sh, lsh, finsh, mf, sigma ) ;
	      o1plistitemadd( ppitem0, iorbi, iorbf, val ) ;
	      counter ++ ;
	    }
	  }
	}
      }
    }
  }
  return counter ;
}

#undef EPSILON
