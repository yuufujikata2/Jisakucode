#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

int w2spama( int nketbasis, struct Fock *pketbasis, int nbrabasis, 
	     struct Fock *pbrabasis, int nrwf, int *lsh, int *sorb1sh, 
	     double *pr0o, struct Spamaline **ppsml0 ) 
{
  int ketj, brai, irwf, i1, sign1, totnumele ;
  double *columnj, val ;
  struct Fock state0, state1 ;
  struct Spamaline *psmcol0 ;
	

  psmcol0 = ( struct Spamaline * )
    malloc( nketbasis * sizeof( struct Spamaline ) ) ;

  columnj = ( double * ) malloc( nbrabasis * sizeof( double ) ) ;
  totnumele = 0 ;

  for ( ketj = 0 ; ketj < nketbasis ; ketj++ ) {
    state0 = pketbasis[ketj] ;
    /* fockdisplay( &state0 ) ; printf("\n") ; */

    /* clear columnj */
    for ( brai = 0 ; brai < nbrabasis ; brai++ ) 
      columnj[brai] = 0. ;

/* apply \sum_{nu,ms} pr0o(nu) c^+_{nu,ms} to state0 */

    for ( irwf = 1 ; irwf < nrwf ; irwf++ ) { /* loop over radial basis */
      for ( i1 = sorb1sh[irwf] ; i1 < sorb1sh[irwf+1] ; i1++ ) { /* ms-loop */
	state1 = state0 ;
	sign1 = fockcre( &state1, i1 ) ;
	if ( sign1 != 0 ) {
	  val = sign1 * pr0o[irwf] ;
	  if( findstate( nbrabasis, pbrabasis, &state1, &brai ) ) {
	    printf("c+%d  -> %lf ", i1, val ); 
	    fockdisplay( &state1 ) ; 
	    printf(" = state No %d\n", brai ) ; 
	    columnj[brai] += val ;
	  }
	  else { 
	    printf("c+%d -> %lf ", i1, val ); 
	    fockdisplay( &state1 ) ;
	    printf(" = new state under No %d\n", brai ) ;
	  } 
	}
      }
    }

    totnumele += darray2spamaline( psmcol0 + ketj, nbrabasis, columnj ) ;
  }


  free( columnj ) ; 

  printf("\n psmcol0[] totnumele = %d.\n\n", totnumele ) ;
  
  for ( ketj = 0 ; ketj < nketbasis ; ketj++ ) {
    printf("ketj =%3d --", ketj ) ;   
    spamalinewrite( psmcol0[ketj] ) ;
  }

  totnumele = spamatranspose( &psmcol0, nketbasis, nbrabasis ) ;


  for ( brai = 0 ; brai < nbrabasis ; brai++ ) {
    printf("brai =%3d --", brai) ;   
    spamalinewrite( psmcol0[brai] ) ;
  }

  *ppsml0 = psmcol0 ;
  return totnumele ;
}
