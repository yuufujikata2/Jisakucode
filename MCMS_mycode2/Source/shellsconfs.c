#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

void makeshells( int nshells, int lcorsh, int lvalsh, double ksicor, 
		 double ksival, int **plsh, int **psorb1sh, double **pksish) 
{
  int *lsh, *sorb1sh, i ;
  double *ksish ;

  lsh     = ( int * ) malloc( nshells * sizeof( int ) ) ;
  sorb1sh = ( int * ) malloc( ( nshells + 1 ) * sizeof( int ) ) ;
  ksish = ( double * ) malloc( nshells * sizeof( double ) ) ;

  sorb1sh[0] = 0 ;
  for ( i = 0 ; i < nshells ; i++ ) {
    lsh[i] = ( i == 0 ? lcorsh : lvalsh ) ;
    sorb1sh[i+1] = sorb1sh[i] + 4 * lsh[i] + 2 ;
  }

  for ( i = 0 ; i < nshells ; i++ ) 
    ksish[i] = ( i == 0 ? ksicor : ksival ) ;
 
  
  printf("nshells = %d, nsorbs = %d\n",	nshells, sorb1sh[nshells] ) ;
  for ( i = 0 ; i < nshells ; i++ ) 
    printf("%3d %3d %8.4lf\n", lsh[i], sorb1sh[i+1], ksish[i] ) ; 

  *plsh = lsh ;
  *psorb1sh = sorb1sh ;
  *pksish = ksish ;
}

void makegsconfs( int nshells, int *lsh, int nvalel, 
		  int *pnconfs, int ***pocc, int *pnstates ) 
{
  int nconfs, nelectrons, nelectronsbefore, nstatesinconf, i, j;
  int **occ, nstates ;

  nconfs = 1 ;
  occ = ( int ** ) malloc( nconfs * sizeof( int * ) ) ;

  nstates = 0 ;
  occ[0] = ( int * ) malloc( nshells * sizeof( int ) ) ;
  
  occ[0][0] = 4 * lsh[0] + 2 ;   /* shell # 0 = core shell */
  occ[0][1] = nvalel ;           /* shell # 1 = valence shell */
  for ( j = 2 ; j < nshells ; j++ )
    occ[0][j] = 0  ;             /* continuum shells  */

  nstatesinconf = 1 ;
  for ( j = 0 ; j < nshells ; j++ )
    nstatesinconf *= noverk( 4 * lsh[j] + 2, occ[0][j] ) ;

  nstates +=  nstatesinconf ;

  *pnconfs = nconfs ;
  *pocc = occ ;
  *pnstates = nstates ;
}

void makefsconfs( int nshells, int *lsh, int nvalel, 
	      int *pnconfs, int ***pocc, int *pnstates ) 
{
  int nconfs, nelectrons, nelectronsbefore, nstatesinconf, i, j;
  int **occ, nstates ;

  nconfs = nshells - 1 ;
  occ = ( int ** ) malloc( nconfs * sizeof( int * ) ) ;

  nstates = 0 ;
  for ( i = 0 ; i < nconfs ; i++ ) {
    occ[i] = ( int * ) malloc( nshells * sizeof( int ) ) ;

    occ[i][0] = 4 * lsh[0] + 1 ;   /* shell # 0 = core shell with one hole */
    occ[i][1] = nvalel + ( i==0 ? 1 : 0 ) ; /* shell # 1 = valence shell */
    for ( j = 2 ; j < nshells ; j++ )
      occ[i][j] = ( j == i+1 ? 1 : 0 ) ;    /* continuum shells occ = 0 or 1 */

    nstatesinconf = 1 ;
    for ( j = 0 ; j < nshells ; j++ )
      nstatesinconf *= noverk( 4 * lsh[j] + 2, occ[i][j] ) ;

    nstates +=  nstatesinconf ;
  }

  *pnconfs = nconfs ;
  *pocc = occ ;
  *pnstates = nstates ;
}
