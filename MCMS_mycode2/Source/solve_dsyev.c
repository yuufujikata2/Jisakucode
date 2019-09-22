#include <stdio.h>
#include <stdlib.h>

#define LWORKFAC 40

int solve_dsyev(int n, double* a, double* lambda ) {
/*
   C-wrapper function for the Lapack routine :

     SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   which solves symmetric eigenvalue problem

   \sum_j a_{ij} v_{jk} = v_{ik} lambda_k 

   where the v_{*k}'s are orthonormal. 
   On output, the array  a  contains the matrix  v.
*/
    char JOBZ = 'V', UPLO = 'L' ;    
    int LWORK, INFO ;
    double *WORK ;
 
    if ( n <= 0 ) return 0 ;

    LWORK = n * LWORKFAC ;
    WORK = ( double * ) calloc( LWORK, sizeof( double ) ) ;

    dsyev_( &JOBZ, &UPLO, &n, a, &n, lambda, WORK, &LWORK, &INFO ) ;

    printf(" LWORK: actual = %d, optimal = %d\n", LWORK, (int)(WORK[0]) ) ;
    free( WORK ) ;

    return INFO ;
}

#undef LWORKFAC
