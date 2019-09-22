#include <stdio.h>
#include <stdlib.h>

#define LWORKFAC 10
int invdcmplxmat( int matdim , double* mat ) {
/*
  inverts general double complex square matrix M_ij = mat[ i + matdim * j ]
  using Lapack routines : zgetrf.f and zgetri.f
      SUBROUTINE ZGETRF(int M, int N, (complex*16) *A, int LDA, int *IPIV, 
                        int INFO)
         LU factorization of M x N matrix A. LDA >= N.
	 IPIV    (output) INTEGER array, dimension (min(M,N)).
         INFO    (output) INTEGER.    = 0: successful exit,  etc
      SUBROUTINE ZGETRI(int N, (complex*16) *A, int LDA, int *IPIV, 
                           (complex*16) *WORK, int LWORK, int INFO )
         computes inverse of matrix A using the LU factorization from ZGETRF
	 WORK is workspace: complex*16 array of dimention LWORK	 
*/
    int LWORK, INFO, i, *IPIV ;
    double *WORK ;
 
    if ( matdim <= 0 ) return 0 ;
    
    IPIV  = ( int * ) calloc( matdim , sizeof( int ) ) ;
    LWORK = matdim * LWORKFAC ;
    WORK  = ( double * ) malloc( LWORK * 2 * sizeof( double ) );/* 2 for Re,Im*/

    zgetrf_( &matdim, &matdim, mat, &matdim, IPIV, &INFO ) ;

    zgetri_( &matdim, mat, &matdim, IPIV, WORK, &LWORK, &INFO ) ;
 
    /*    printf(" zgetri: LWORK: actual = %d, optimal = %d\n", 
	   LWORK, (int)(WORK[0]) ) ;  
    */
    free( WORK ) ;
    free( IPIV ) ;

    return INFO ;
}

