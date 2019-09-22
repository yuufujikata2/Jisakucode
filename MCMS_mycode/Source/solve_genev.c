#include <stdio.h>
#include <stdlib.h>

#define LWORKFAC 40

int solve_genev(int n, double* a, double* b,
	        double* alpha, double* beta, double* revec) {
/*
  routine solves generalized eigenvalue problem

   \sum_j a_{ij} revec_{jk} = \sum_j b_{ij} revec_{jk} alpha_k/beta_k

   even "solutions" with beta_k = 0 are calculated
*/
/*  uses Lapack routine :
    SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
*/
    char JOBVL = 'N', JOBVR = 'V' ;    
    int LDVL = 1, LWORK, INFO, i ;
    double *ALPHAI, *VL, *WORK ;
 
    if ( n <= 0 ) return 0 ;

    ALPHAI = ( double * ) calloc( n, sizeof( double ) ) ;
    VL = ( double * ) calloc( LDVL * n, sizeof( double ) ) ;
    LWORK = n * LWORKFAC ;
    WORK = ( double * ) calloc( LWORK, sizeof( double ) ) ;

    dggev_( &JOBVL, &JOBVR, &n, a, &n, b, &n, alpha, ALPHAI, 
	    beta, VL, &LDVL, revec, &n, WORK, &LWORK, &INFO ) ;

/*  printf(" LWORK: actual = %d, optimal = %d\n", LWORK, (int)(WORK[0]) ) ;*/
/*  print out nonzero ALPHAI_k's (should not exist) */
    for ( i = 0 ;  i < n ; i++ )
	if ( ALPHAI[i] != 0. ) 
	    printf(" ALPHAI[%d] = %lg != 0.\n", i, ALPHAI[i]);
    free( WORK ) ;
    free( VL ) ;
    free( ALPHAI ) ; 

    return INFO ;
}



