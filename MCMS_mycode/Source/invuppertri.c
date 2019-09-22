#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* 
   invert "manually" the _regular_ upper triangular square matrix mat[][]
   (conditions:   mat[i][i] != 0   and   mat[i][j] == 0 if i > j )
   allocate memory and return (through pointer) the inverse matrix invmat
*/  
int invuppertri(int dim, double **mat, double ***pinvmat )
{
  int i, j, k ;
  double **invmat, sum ;

  /* allocation */
  invmat = ( double ** ) malloc( dim * sizeof( double * ) ) ;
  for ( i = 0 ; i < dim ; i++ )
    invmat[i] = ( double * ) calloc( dim , sizeof( double ) ) ;

  for ( i = dim - 1 ; i >= 0 ; i-- ) {
    if ( fabs(mat[i][i]) < 1.e-4 ) {
      printf("invuppertri: BIG DANGER: fabs(mat[i][i]) < 1.e-4 \n\n EXIT\n") ;
      exit(1) ;
    }
    for ( j = 0 ; j < dim ; j++ ) {
      sum = 0. ;
      for ( k = i + 1 ; k < dim ; k++ )
	sum += mat[i][k] * invmat[k][j] ;
      invmat[i][j] = ( ( i==j ? 1. : 0. ) - sum ) / mat[i][i] ;
    }
  }

  *pinvmat = invmat ;

  printf("\n invuppertri.c :  invmat[i][j]\n\n") ;
  for ( i = 0 ; i < dim ; i++ ) {
    for ( j = 0 ; j < dim ; j++ ) 
      printf(" %10.6lf", invmat[i][j] ) ;
    printf("\n") ;
  }

  return 1 ;
} 


