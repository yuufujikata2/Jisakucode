#include <stdio.h>
#include <stdlib.h>
#include "globals.h"
int h1pq2list( int nshells, int *lsh, int *sorb1sh, double **h1pq,
	       struct O1plistitem **ppo1plistitem0 ) ;
/* calcham1p
   Calculates one-particle (hamiltonian, or other) operator.
   special simple case of the algorithm for calcham, 
   which calculates the hamiltonian matrix, and stores it as a sparse
   matrix named " ham == *pham ". ( see sparsemat.c for details )
   ham is a pointer to (array of) struct sparsematline. 
   Return value : number of stored (i.e. non-zero) matrix elements.
   All memory allocation is done here.
   The matrix elements are stored as follows :
   < state[i] | H | state[j] > = ham[i].v[k], where i = i and j = ham[i].j[k]
*/
int calcham1p( struct Spamaline **pham, int nstates, struct Fock *state,
	     int nconfs, int nshells, int *lsh, int *sorb1sh, int **occ, 
	     double **h1pq, int nr, double *rofi, double **ppofi ) 
{
  int k1, k2, totnumele ;
  int i1, i2, sign1, sign2 ;
  int ist, jst ;
  int newstatecounter = 0 ;
  int  no1plistitems ;
  double val, *hamcol, diagsum ;
  struct Fock state0, state1, state2 ;
  struct O1plistitem *po1plistitem0 ;
  struct O1p h1p ;
  struct O1p2 h1p1 ;
  struct Spamaline *hamtransp ;


  po1plistitem0 = 0 ;

  no1plistitems = 0 ;  

  no1plistitems += 
    h1pq2list( nshells, lsh, sorb1sh, h1pq, &po1plistitem0 ) ;

  printf("no1plistitems = %d\n", no1plistitems ) ;

  h1p = o1pmake( po1plistitem0 ) ;

  o1plistdelete( po1plistitem0 ) ;

  o1pprint( h1p ) ;

  hamtransp = ( struct Spamaline * )
    malloc( nstates * sizeof( struct Spamaline ) ) ;

  hamcol = ( double * ) malloc( nstates * sizeof( double ) ) ;
  totnumele = 0 ;
  diagsum = 0. ;
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    state0 = state[ist] ;
    /* fockdisplay( &state0 ) ; printf("\n") ; */

    /* clear hamcol */
    for ( jst = 0 ; jst < nstates ; jst++ ) 
      hamcol[jst] = 0. ;

/* apply h1p to state */
    for ( k1 = 0 ; k1 < h1p.n ; k1++ ) {
      state1 = state0 ;
      i1 = h1p.i[k1] ;
      sign1 = fockdes( &state1, i1 ) ;
      if ( sign1 != 0 ) {
	h1p1 = h1p.p[k1] ;
	for ( k2 = 0 ; k2 < h1p1.n ; k2++ ) {
	  state2 = state1 ;
	  i2 = h1p1.i[k2] ;
	  sign2 = fockcre( &state2, i2 ) ;
	  if ( sign2 != 0 ) {
	    val = sign2 * sign1 *  h1p1.v[k2] ;
	    /* printf("c+%d c%d  -> %lf ", i2, i1, val ); 
	       fockdisplay( &state2 ) ; */
	    if( findstate(nstates, state, &state2, &jst ) ) {
	      /* printf(" = state No %d\n", jst ) ; */
	      hamcol[jst] += val ;
	    } 
	    else { 
	      newstatecounter++ ;
	      /*
	      printf("c+%d c%d  -> %lf ", i2, i1, val ); 
	      fockdisplay( &state2 ) ;
	      printf(" = new state under No %d\n", jst ) ;
	      */
	    } 
	  }
	}
      }
    }
    totnumele += darray2spamaline( &hamtransp[ist], nstates, hamcol ) ;
    diagsum += hamcol[ist] ;
  }

  printf("Application of the operator generated %d new states\n", 
	 newstatecounter );

  o1pdelete( h1p ) ;

  free( hamcol ) ; 

  /*
  printf("\n qmetrictransp. totnumele = %d. diagaverage = %lf \n\n",  
	 totnumele, diagsum / ( double ) nstates ) ;
  */
  
  totnumele = spamatranspose( &hamtransp, nstates, nstates ) ;
 
  printf("\n after spamatranspose qmetric totnumele = %d\n",  totnumele ) ;
  
  /*
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    printf("i =%3d --", ist ) ;   
    spamalinewrite( hamtransp[ist] ) ;
  }
  */

  *pham = hamtransp ;

  return totnumele ;
}
