#include <stdio.h>
#include <math.h>
#include "globals.h"
#define EPSILON 1.e-10

void printstatevector( int nstates, double *vec, struct Fock *state ) 
{
  int i ;
  for ( i = 0 ; i < nstates ; i++ ) {
    /* long version */
    if ( fabs(vec[i]) > EPSILON ) {
      printf("%15.10lf ", vec[i] ) ;
      fockdisplay( &state[i] ) ;
      printf("\n") ;
    }
  }
  printf("\n") ;
}

#undef EPSILON
