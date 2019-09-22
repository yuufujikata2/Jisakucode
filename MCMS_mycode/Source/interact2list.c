#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "w3j.h"

int sporbindex( int *sorb1sh, int *lsh, int ish, int m, int sigma )
{
  return sorb1sh[ish] + 2 * ( lsh[ish] + m ) + ( sigma > 0 ? 1 : 0 ) ;
}

double cklmlm( int k, int l1, int m1, int l2, int m2 ) {
  if ( !w3jdefined( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) ||
       !w3jnonzero( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) )
    return 0. ;
  else 
    return ( m1 % 2 == 0 ? 1 : -1 )              /* = (-)^m == (-)^{-m} */
      * sqrt( ( 2 * l1 + 1 ) * ( 2 * l2 + 1 ) )  /* = [l1,l2]^{1/2} */
      * w3j( 2*l1, 2*k, 2*l2, 0, 0, 0 )
      * w3j( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) ; 
}


int coulomb2list( int nshells, int *lsh, int *sorb1sh, 
		  struct CIlistitem *pcilistitem0, 
		  struct O2plistitem **ppo2plistitem0 )
{
  int ish1, ish2, ish3, ish4, l1, l2, l3, l4, m1, m2, m3, m4, s1, s2, s3, s4 ;
  int i1, i2, i3, i4, nk, k, counter = 0 ;
  double *r4312, v0 ;
  struct CIlistitem *ip = pcilistitem0 ;

  *ppo2plistitem0 = 0 ;

  while( ip != 0 ) {
    ish1 = ip->ish1 ; ish2 = ip->ish2 ; ish3 = ip->ish3 ; ish4 = ip->ish4 ; 
    l1 = lsh[ish1] ; l2 = lsh[ish2] ; l3 = lsh[ish3] ; l4 = lsh[ish4] ;
    nk = ip -> nk ; /* nk=min(l4+l1,l2+l3)+1. ck(4,1)ck(2,3)!=0 for k=0..nk-1*/
    r4312 = ip -> rmx ; /* r4312[k] defined for k=0..nk-1 */

    for ( m1 = -l1 ; m1 <= l1 ; m1++ ) 
      for ( m2 = -l2 ; m2 <= l2 ; m2++ )
	for ( m3 = -l3 ; m3 <= l3 ; m3++ ) { 
	  m4 = m1 + m2 - m3 ;
	  if ( m4 >= -l4 && m4 <= l4 ) 
	    for ( s1 = -1 ; s1 <= 1 ; s1 += 2 )
	      for ( s2 = -1 ; s2 <= 1 ; s2 += 2 ) {
		s4 = s1 ; 
		s3 = s2 ;
		i1 = sporbindex( sorb1sh, lsh, ish1, m1, s1 ) ;
		i2 = sporbindex( sorb1sh, lsh, ish2, m2, s2 ) ;
		i3 = sporbindex( sorb1sh, lsh, ish3, m3, s3 ) ;
		i4 = sporbindex( sorb1sh, lsh, ish4, m4, s4 ) ;
		if ( i1 < i2 && i3 != i4 ) { 
		  v0 = 0. ;
		  for ( k = 0 ; k < nk ; k++ )
		    v0 += cklmlm(k,l4,m4,l1,m1) * cklmlm(k,l2,m2,l3,m3) 
		          * r4312[k] ;
		  if ( v0 != 0. ) {
		    o2plistitemadd( ppo2plistitem0, i1, i2, i3, i4, v0 ) ;
		    counter++ ;
		  }
		}

	      }
	}
    ip = ip -> next ;
  }
  return counter ;
}

int spinorbit2list( int nshells, int *lsh, int *sorb1sh, double *ksish, 
		    struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, i2, counter = 0 ;
  double ksi0, v0 ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    ksi0 = ksish[ish] ;
    if ( l0 != 0 && ksi0 != 0. ) {
      /* lz sz term (diagonal) */
      for ( m = -l0 ; m <= l0 ; m++ )
	if ( m != 0 ) {
	  v0 = 0.5 * m * ksi0 ;
	  for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	    i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
	    o1plistitemadd( ppo1plistitem0, i1, i1, sigma * v0 ) ;
            /*	printf("%2d %2d %20.12e\n", i1, i1, sigma * v0 ) ; */
	    counter++ ;
	  }
	}
      /* 0.5 ( l+ s- + l- s+ ) term (off-diagonal) */
      for ( m = -l0 ; m < l0 ; m++ ) {
	v0 = 0.5 * sqrt( ( l0 - m ) * ( l0 + m + 1 ) ) * ksi0 ;
	i1 = sporbindex( sorb1sh, lsh, ish, m,    1 ) ; 
	i2 = sporbindex( sorb1sh, lsh, ish, m+1, -1 ) ; 
	o1plistitemadd( ppo1plistitem0, i1, i2, v0 ) ;
	o1plistitemadd( ppo1plistitem0, i2, i1, v0 ) ;
	/*  printf("%2d %2d %20.12e (l+s-) and ", i1, i2, v0 ) ;
	    printf("%2d %2d %20.12e (l-s+)\n",    i2, i1, v0 ) ;  */
	counter += 2 ;
      }
    }
  }
  return counter ;
}

/* calculates  magnetic field operator = -hz * 2 * sz */
int magneticfield2list( int nshells, int *lsh, int *sorb1sh, double hz,
			struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, counter = 0 ;

  if ( hz != 0. ) { 
    for ( ish = 0 ; ish < nshells ; ish++ ) {
      l0   = lsh[ish] ; 
      for ( m = -l0 ; m <= l0 ; m++ )
	for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	  i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
/* case of magnetic field operator = -hz * 2 * sz */
	  o1plistitemadd( ppo1plistitem0, i1, i1, -hz * sigma ) ;
/* optionnally: case of magnetic field operator = -hz * jz */
/*        o1plistitemadd( ppo1plistitem0, i1, i1, -hz *( m + 0.5*sigma ) ) ;*/
	  counter++ ;
	}
    }
  }
  return counter ;
}

/* applies crystal field in Oh symmetry */
int crystalfieldOh2list( int nshells, int *lsh, int *sorb1sh, 
			 double *tendqsh, struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, i2, counter = 0 ;
  double tendq, shift ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    tendq = tendqsh[ish] ;
    if ( tendq != 0. && l0 > 1 ) {
      if ( l0 == 2 ) {
      /* diagonal elements: 0:eg. 1,-1:t2g,  <2|CF|2>=<-2|CF|-2>=(eg+t2g)/2 */
        for ( m = -l0 ; m <= l0 ; m++ ) {
	  if ( m == 0 ) shift = 0.6 * tendq ;
	  else if ( m == 1 || m == -1 ) shift = -0.4 * tendq ;
	  else if ( m == 2 || m == -2 ) shift =  0.1 * tendq ;
	  else { printf("TROUBLE in crystalfieldOh2list\n") ; exit(1) ; }
	  for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	    i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
	    o1plistitemadd( ppo1plistitem0, i1, i1, shift ) ;
	    counter++ ;
	  }      
	}
	/* off diagonal elements: <2|CF|-2> = <-2|CF|2> = (eg-t2g)/2 */
	shift = 0.5 * tendq ;
	for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	  i1 = sporbindex( sorb1sh, lsh, ish, 2, sigma ) ;
	  i2 = sporbindex( sorb1sh, lsh, ish,-2, sigma ) ;
	  o1plistitemadd( ppo1plistitem0, i1, i2, shift ) ;
	  o1plistitemadd( ppo1plistitem0, i2, i1, shift ) ;
	  counter += 2 ;
	}      
      }
      else { /* l > 2 */
	printf("CF not applied to shell No %d. Don't know how for l = %d\n",
	       ish, l0  ) ;
      }
    }
  }
  return counter ;
}


/*  writes radial h1pq[ish][jsh] == <p_ish|h(1)|p_jsh> to interaction list */
int h1pq2list( int nshells, int *lsh, int *sorb1sh, double **h1pq,
	       struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, jsh, l0, li, lj, m, sigma, i1, i2, counter = 0 ;
  double value ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    for ( jsh = 0 ; jsh < nshells ; jsh++ ) {
      value = h1pq[ish][jsh] ;
      if ( value != 0. ) {
	li = lsh[ish] ; lj = lsh[jsh] ;
	l0 = ( li < lj ? li : lj ) ;       /* l0 = min(li,lj) */
	for ( m = -l0 ; m <= l0 ; m++ )
	  for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	    i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ;
	    i2 = sporbindex( sorb1sh, lsh, jsh, m, sigma ) ; 
	    o1plistitemadd( ppo1plistitem0, i1, i2, value ) ;
	    counter++ ;
	  }
      }
    }
  }
  return counter ;
}


