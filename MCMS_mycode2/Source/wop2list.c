#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#define EPSILON 1.e-10

int wop2list( int nshells, int *lsh, int *sorb1sh, int wsh, double *pr0o,
	      struct O1plistitem **ppitem0 )
{
  int ish, i10, i20, ims ;
  int counter = 0 ;
  double val ;
  
  /*
  printf("wop2list. Attention! Very unportable code.\n-- Only works if ") ;
  printf("coresh is first, then only valence shells with the same l\n") ;
  */

  if ( nshells < 3 ) { printf("wop2list. PANIC: nshells < 3\n"); exit(1); }
  if ( wsh < 2 ) { printf("wop2list. PANIC: wsh < 2\n"); exit(1); }

  i10 = sorb1sh[wsh] ;
  for ( ish = 1 ; ish < nshells ; ish++ ) {
    if ( lsh[ish] != lsh[wsh] ) { 
      printf("wop2list. PANIC: lsh[ish] != lsh[wsh]\n"); exit(1); 
    }
    val = pr0o[ish] ;
    if ( fabs( val ) > EPSILON ) {
      i20 = sorb1sh[ish] ;
      for ( ims = 0 ; ims < 4 * lsh[ish] + 2 ; ims++ ) {
        o1plistitemadd( ppitem0, i10 + ims, i20 + ims, val ) ;
        counter ++ ;
      }
    }
  }
  return counter ;
}

#undef EPSILON
