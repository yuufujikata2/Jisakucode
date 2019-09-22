/* which == 0, -1, 1  for ground state, channel function, final state, resp */
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
void makeconfs( int nshells, int *lsh, int nvalel, int which,
		int *pnconfs, int ***pocc, int *pnstates ) 
{
  int nconfs, nstatesinconf, i, j;
  int **occ, nstates, corehole, photoel ;


  if ( which == 1 && nvalel > 4 * lsh[1] + 1 ) {
    printf("makeconf PANIC:  nvalel > 4 * lsh[1] + 1\n") ;  exit(1) ; 
  }
  corehole = ( which == 0 ? 0 : 1 ) ;
  photoel  = ( which == 1 ? 1 : 0 ) ;

  nconfs = ( which > 0 ? nshells - 1 : 1 ) ;
  occ = ( int ** ) malloc( nconfs * sizeof( int * ) ) ;

  nstates = 0 ;
  for ( i = 0 ; i < nconfs ; i++ ) {
    occ[i] = ( int * ) malloc( nshells * sizeof( int ) ) ;

    occ[i][0] = 4 * lsh[0] + 2 - corehole;      /* shell # 0 = core shell    */
    
    occ[i][1] = nvalel + ( i==0 ? photoel : 0 ); /* shell # 1 = valence shell */
    for ( j = 2 ; j < nshells ; j++ )
      occ[i][j] = ( j == i+1 ? photoel : 0 ) ;  /* cont. shells occ = 0 or 1 */

    nstatesinconf = 1 ;
    for ( j = 0 ; j < nshells ; j++ )
      nstatesinconf *= noverk( 4 * lsh[j] + 2, occ[i][j] ) ;

    nstates +=  nstatesinconf ;
  }

  *pnconfs = nconfs ;
  *pocc = occ ;
  *pnstates = nstates ;
}
