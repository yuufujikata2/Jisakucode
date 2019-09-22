#include <stdio.h>
#include <stdlib.h>
#define LWORKFAC 10

int inv_remat( int matdim, double* mat ) {
/*
  inverts general square matrix M_ij = mat[ i + matdim * j ]
  using Lapack routines : dgetrf.f and dgetri.f
      subroutine  dcopy(int n, double dx, int incx=1 , double dy, int incy=1)
         copies the vector dx, to a vector, dy, both of length n
      SUBROUTINE DGETRF(int M, int N, double* A, int LDA, int* IPIV, int INFO)
         LU factorization of M x N matrix A. LDA >= N.
	 IPIV    (output) INTEGER array, dimension (min(M,N)).
         INFO    (output) INTEGER.    = 0: successful exit,  etc
      SUBROUTINE DGETRI(int N, double* A, int LDA, int* IPIV, double* WORK, 
                        int LWORK, int INFO )
         computes inverse of matrix A using the LU factorization from DGETRF
	 WORK is workspace: double array of dimention LWORK	 
*/
    int LWORK, INFO, i, *IPIV ;
    double *WORK ;
 
    if ( matdim <= 0 ) return 0 ;
    
    IPIV  = ( int * ) calloc( matdim , sizeof( int ) ) ;
    LWORK = matdim * LWORKFAC ;
    WORK  = ( double * ) calloc( LWORK, sizeof( double ) ) ;

    dgetrf_( &matdim, &matdim, mat, &matdim, IPIV, &INFO ) ;

    dgetri_( &matdim, mat, &matdim, IPIV, WORK, &LWORK, &INFO ) ;
 
    /*    printf(" dgetri: LWORK: actual = %d, optimal = %d\n", 
	   LWORK, (int)(WORK[0]) ) ;
    */
    free( WORK ) ;
    free( IPIV ) ;

    return INFO ;
}
