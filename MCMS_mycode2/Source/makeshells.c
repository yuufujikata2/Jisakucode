#include <stdio.h>
#include <stdlib.h>

void makeshells( int nshells, int lcorsh, int lvalsh, double ksicor, 
		 double ksival, double tendqval, int **plsh, int **psorb1sh, 
                 double **pksish, double **ptendqsh) 
{
  int *lsh, *sorb1sh, i ;
  double *ksish, *tendqsh ;

  lsh     = ( int * ) malloc( nshells * sizeof( int ) ) ;
  sorb1sh = ( int * ) malloc( ( nshells + 1 ) * sizeof( int ) ) ;
  ksish = ( double * ) calloc( nshells, sizeof( double ) ) ;
  tendqsh = ( double * ) calloc( nshells, sizeof( double ) ) ;

  sorb1sh[0] = 0 ;
  for ( i = 0 ; i < nshells ; i++ ) {
    lsh[i] = ( i == 0 ? lcorsh : lvalsh ) ;
    sorb1sh[i+1] = sorb1sh[i] + 4 * lsh[i] + 2 ;
  }

  ksish[0] = ksicor ;
  ksish[1] = ksival ;
  
  tendqsh[1] = tendqval ;
  
  printf("nshells = %d, nsorbs = %d\n",	nshells, sorb1sh[nshells] ) ;
  printf(" i    lsh[i]    sorb1sh[i+1]   ksish[i]   tendqsh[i] \n") ;
  for ( i = 0 ; i < nshells ; i++ ) 
    printf("%3d %3d %3d %8.4lf %8.4lf\n", 
           i, lsh[i], sorb1sh[i+1], ksish[i], tendqsh[i] ) ; 

  *plsh = lsh ;
  *psorb1sh = sorb1sh ;
  *pksish = ksish ;
  *ptendqsh = tendqsh ;
}
