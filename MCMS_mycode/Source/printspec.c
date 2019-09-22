#include <stdio.h>
#include <math.h>
#define EPSILON 1.e-10

int printspec( int dim, double *lambda, int *pgstdeg, double *pgstenergy ) 
{
    int nlevels, degeneracy, i, gstdeg ;
    double gstenergy, newlambda, oldlambda ;
    nlevels = 0 ;
    for ( i = 0 ; i < dim ; i++ ) {
      newlambda = lambda[i] ;
    if ( nlevels == 0 ) {
      oldlambda = newlambda ;
      degeneracy = 1 ;
      nlevels = 1 ;
    } 
    else if ( fabs( newlambda - oldlambda ) < EPSILON )
      degeneracy++ ;
    else {
      printf(" %12.6lf  -  degeneracy = %d\n", oldlambda, degeneracy ) ;
      if ( nlevels == 1 ) {
	gstenergy = oldlambda ;
	gstdeg = degeneracy ;
      }
      oldlambda = newlambda ;
      degeneracy = 1 ;
      nlevels++ ;
    }
  }      
  printf(" %12.6lf  -  degeneracy = %d\n", oldlambda, degeneracy ) ;
  if ( nlevels == 1 ) {
    gstenergy = oldlambda ;
    gstdeg = degeneracy ;
  }

  *pgstdeg = gstdeg ;
  *pgstenergy = gstenergy ;
  return nlevels ;
}

#undef EPSILON
