#include <stdio.h>
#include <stdlib.h>

int transfe0pr0dpr0( int nrwf, double **smat, double **tmat, 
		     double *e0p, double *pr0, double *dpr0, 
		     double ***pe0po, double **ppr0o, double **pdpr0o ) 
{
  int i, j, k ;
  double **e0po, *pr0o, *dpr0o, sum, sumd ;
  
/* allocation */
  e0po = ( double ** ) malloc( nrwf * sizeof( double * ) ) ;
  for ( i = 0 ; i < nrwf ; i++ ) 
    e0po[i]  = ( double * ) calloc( nrwf, sizeof( double ) ) ;
  pr0o  = ( double * ) calloc( nrwf, sizeof( double ) ) ;
  dpr0o = ( double * ) calloc( nrwf, sizeof( double ) ) ;


  /* calculate e0po[i][j], i.e. h^(1) matrix in new radial basis */

  for ( i = 0 ; i < nrwf ; i++ ) {
    for ( j = 0 ; j < nrwf ; j++ ) {
      sum = 0. ;
      for ( k = 0 ;  k < nrwf ; k++ )
	sum += smat[i][k] * e0p[k] * tmat[k][j] ;
      e0po[i][j] = sum ;
    }
  }

  for ( i = 0 ; i < nrwf ; i++ ) {   
    sum  = 0. ;
    sumd = 0. ;
    for ( k = 0 ;  k < nrwf ; k++ ) {
	sum  +=  pr0[k] * tmat[k][i] ;
	sumd += dpr0[k] * tmat[k][i] ;
    }
    pr0o[i]  = sum ;
    dpr0o[i] = sumd ;
  }    

  printf("\nh^(1) matrix in orthonormal radial basis: e0po[i][j]\n") ;
  for ( i = 0 ; i < nrwf ; i++ ) {
    for ( j = 0 ; j < nrwf ; j++ )
      printf(" %10.6lf", e0po[i][j] ) ;
    printf("\n") ;
  }

  printf("\n pr0o[i] = ") ;
  for ( i = 0 ; i < nrwf ; i++ ) 
    printf(" %10.6lf", pr0o[i] ) ;
  printf("\ndpr0o[i] = ") ;
  for ( i = 0 ; i < nrwf ; i++ ) 
    printf(" %10.6lf", dpr0o[i] ) ;
  printf("\n") ;
  
  printf("\n Q-matrix in orthonormal radial basis\n") ;
  for ( i = 0 ; i < nrwf ; i++ ) {
    for ( j = 0 ; j < nrwf ; j++ )
      printf(" %10.6lf", pr0o[i]*pr0o[j] ) ;
    printf("\n") ;
  }

  printf("\n L-matrix in orthonormal radial basis\n") ;
  for ( i = 0 ; i < nrwf ; i++ ) {
    for ( j = 0 ; j < nrwf ; j++ )
      printf(" %10.6lf", pr0o[i]*dpr0o[j] ) ;
    printf("\n") ;
  }

  printf("\n h+L\n") ;
  for ( i = 0 ; i < nrwf ; i++ ) {
    for ( j = 0 ; j < nrwf ; j++ )
      printf(" %10.6lf", e0po[i][j] + pr0o[i]*dpr0o[j] ) ;
    printf("\n") ;
  }

  *pe0po = e0po ; *ppr0o = pr0o ; *pdpr0o = dpr0o ;

  return 1 ;
} 
