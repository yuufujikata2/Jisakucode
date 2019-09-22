#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

int crystalfieldOh2list( int nshells, int *lsh, int *sorb1sh, 
			 double *tendqsh, struct O1plistitem **ppo1plistitem0 ) ;


int isaconf( int nshells, int *lsh, int *occ ) ;
int confinset( int nconfs, int nshells, int **occ, int *occdum ) ;
int cilistitemadd( struct CIlistitem **ip0, int j1, int j2, int j3, int j4,
		   int *lsh, int nr, double *rofi, double **ppofi,
		   double rkscale ) ;
void cilistdelete( struct CIlistitem *ip ) ;
double slaterrk(int nr, double* rofi, int k, double* pi, double* pj, 
                double* pt, double* pu) ;
int h1pq2list( int nshells, int *lsh, int *sorb1sh, double **h1pq,
	       struct O1plistitem **ppo1plistitem0 ) ;


/* calcham 
   calculates the hamiltonian matrix, and stores it as a sparse
   matrix named " ham == *pham ". ( see sparsemat.c for details )
   ham is a pointer to (array of) struct sparsematline. 
   Return value : number of stored (i.e. non-zero) matrix elements.
   All memory allocation is done here.
   The matrix elements are stored as follows :
   < state[i] | H | state[j] > = ham[i].v[k], where i = i and j = ham[i].j[k]
*/
int calcham( struct Spamaline **pham, int nstates, struct Fock *state,
	     int nconfs, int nshells, int *lsh, int *sorb1sh, int **occ, 
	     double **h1pq, double *ksish, double *tendqsh, double hz, 
             int nr, double *rofi, double **ppofi, double rkscale )
{
  int i, ish, iconf, k1, k2, k3, k4, totnumele, k ;
  int i1, i2, i3, i4, sign1, sign2, sign3, sign4 ;
  int ist, jst ;
  int *occdum, ish1, ish2, ish3, ish4 ;
  int  no1plistitems, no2plistitems ;
  int newstatecounter = 0 ;
  double val, *hamcol, diagsum ;
  double confavfpd ;
  struct Fock state0, state1, state2, state3, state4 ;
  struct O1plistitem *po1plistitem0 ;
  struct O1p h1p ;
  struct O1p2 h1p1 ;
  struct O2plistitem *po2plistitem0, *lip ;
  struct O2p vcb ;
  struct O2p2 vcb1 ;
  struct O2p3 vcb2 ;
  struct O2p4 vcb3 ;
  struct CIlistitem *interactlist0, *ip ;
  struct Spamaline *hamtransp ;


  po1plistitem0 = 0 ;

  no1plistitems = 0 ;  
 
  no1plistitems += 
    h1pq2list( nshells, lsh, sorb1sh, h1pq, &po1plistitem0 ) ;

  no1plistitems += 
    spinorbit2list( nshells, lsh, sorb1sh, ksish, &po1plistitem0 ) ;

  no1plistitems += 
    magneticfield2list( nshells, lsh, sorb1sh, hz, &po1plistitem0 ) ;

  printf("A magnetic field hz = %lf ( with E = -hz*2*sz ) is applied.\n", hz) ;
 
  no1plistitems += 
    crystalfieldOh2list( nshells, lsh, sorb1sh, tendqsh, &po1plistitem0 ) ;

  printf("no1plistitems = %d\n", no1plistitems ) ;

  h1p = o1pmake( po1plistitem0 ) ;

  o1plistdelete( po1plistitem0 ) ;

  o1pprint( h1p ) ;


  interactlist0 = 0 ;

  occdum = ( int * ) malloc( nshells * sizeof( int ) ) ;

  for ( iconf = 0 ; iconf < nconfs ; iconf++ ) {

    for ( ish = 0 ; ish < nshells ; ish++ ) 
      occdum[ish] = occ[iconf][ish] ;

    for ( ish1 = 0 ; ish1 < nshells ; ish1++ ) {
      occdum[ish1]-- ;
      if ( isaconf( nshells, lsh, occdum ) )
	for ( ish2 = 0 ; ish2 < nshells ; ish2++ ) {
	  occdum[ish2]-- ;
	  if ( isaconf( nshells, lsh, occdum ) )
	    for ( ish3 = 0 ; ish3 < nshells ; ish3++ ) {
	      occdum[ish3]++ ;
	      if ( isaconf( nshells, lsh, occdum ) )
		for ( ish4 = 0 ; ish4 < nshells ; ish4++ ) {
		  occdum[ish4]++ ;
		  if ( isaconf( nshells, lsh, occdum ) ) 
		    if ( confinset( nconfs, nshells, occ, occdum ) )
		      cilistitemadd( &interactlist0,ish1,ish2,ish3,ish4,lsh,
				     nr,rofi,ppofi, rkscale ) ;
		  occdum[ish4]-- ;
		}
	      occdum[ish3]-- ;
	    }
	  occdum[ish2]++ ;
	}
      occdum[ish1]++ ;
    }
  }
  free( occdum ) ;

  printf("\nList of generalized Slater integrals\n");
  printf("Attention! F^0_ij == <ij|R^0|ij> is set to zero\n");
  printf(" #) i j t u <i,j|R^k|t,u>, k = 0,1,..\n");
  confavfpd = 0. ;
  for( ip = interactlist0, i = 0 ; ip != 0 ; ip = ip -> next ) 
    if ( ip->ish1 == 0 && ip->ish2 == 1 && ip->ish3 == 0 && ip->ish4 == 1 )
      if ( lsh[0] != 1 || lsh[1] != 2 || ip -> nk != 4 ) {
        printf("calcham: lsh[0] != 1 || lsh[1] != 2 || ip -> nk != 4\n") ;
        printf("Set confavfpd to zero TO BE CORRECTED ! \n") ;
	confavfpd = 0. ;
/*        exit(1) ; */
      }  
      else    
        confavfpd = (1./15.)*ip->rmx[1] + (3./70.)*ip->rmx[3] ;
  for( ip = interactlist0, i = 0 ; ip != 0 ; ip = ip -> next ) {
    printf("%2d) %d %d %d %d",++i, ip->ish1, ip->ish2, ip->ish4, ip->ish3);
    if ( ip -> nk  >  0 )
      ip -> rmx[0] = 0. ;
    /* EXTRA FOR EXPERIMENT: put to 0 all Rk except for core-3d */
    /*
    if ( ip->ish1 > 1 || ip->ish2 > 1 || ip->ish3 > 1 || ip->ish4 > 1 )
      for ( k = 0 ; k < ip-> nk ; k++ )
        ip -> rmx[k] = 0. ;
    */
    if ( ip->ish1 == ip->ish4 &&  ip->ish2 == ip->ish3 && ip -> nk > 0 ) 
      if ( ip->ish1 == 0 && ip->ish2 == 1 || ip->ish1 == 1 && ip->ish2 == 0 )
        ip-> rmx[0] = confavfpd ;
      else if ( ip->ish1 == 1 && ip->ish2 == 1 )
        ip-> rmx[0] = (2./63.)*ip->rmx[2] + (2./63.)*ip->rmx[4] ;
      else if ( ip->ish1 == 0 && ip->ish2 == 0 )
        ip-> rmx[0] = (2./25.)*ip->rmx[2] ;
    for ( k = 0 ; k < ip-> nk ; k++ )
      printf(" %lf", ip -> rmx[k] ) ;
    printf("\n") ;
  } 
  

  po2plistitem0 = 0 ;

  no2plistitems = 
    coulomb2list( nshells, lsh, sorb1sh, interactlist0, &po2plistitem0 ) ;

  printf("no2plistitems = %d\n", no2plistitems ) ;

  cilistdelete( interactlist0 ) ; 

  printf("just before o2pmake\n") ;
  vcb = o2pmake( po2plistitem0 ) ;

  o2plistdelete( po2plistitem0 ) ;

  /*
  o2pprint( vcb ) ;
  */  

  hamtransp = ( struct Spamaline * )
    malloc( nstates * sizeof( struct Spamaline ) ) ;

  hamcol = ( double * ) malloc( nstates * sizeof( double ) ) ;
  totnumele = 0 ;
  diagsum = 0. ;
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    state0 = state[ist] ;
    /*
    fockdisplay( &state0 ) ; printf("\n") ;
    */
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
	    /*
	    printf("c+%d c%d  -> %lf ", i2, i1, val ); 
     	    fockdisplay( &state2 ) ;  
	    */
	    if( findstate(nstates, state, &state2, &jst ) ) {
	      /*	      printf(" = state No %d\n", jst ) ; */ 
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
/* apply vcb to state */
    for ( k1 = 0 ; k1 < vcb.n ; k1++ ) {
      state1 = state0 ;
      i1 = vcb.i[k1] ;
      sign1 = fockdes( &state1, i1 ) ;
      if ( sign1 != 0 ) {
	vcb1 = vcb.p[k1] ;
	for ( k2 = 0 ; k2 < vcb1.n ; k2++ ) {
	  state2 = state1 ;
	  i2 = vcb1.i[k2] ;
	  sign2 = fockdes( &state2, i2 ) ;
	  if ( sign2 != 0 ) {
	    vcb2 = vcb1.p[k2] ;
	    for ( k3 = 0 ; k3 < vcb2.n ; k3++ ) {
	      state3 = state2 ;
	      i3 = vcb2.i[k3] ;
	      sign3 = fockcre( &state3, i3 ) ;
	      if ( sign3 != 0 ) {
		vcb3 = vcb2.p[k3] ;
		for ( k4 = 0 ; k4 < vcb3.n ; k4++ ) {
		  state4 = state3 ;
		  i4 = vcb3.i[k4] ;
		  sign4 = fockcre( &state4, i4 ) ;
		  if ( sign4 != 0 ) {
		    val = sign4 * sign3 * sign2 * sign1 *  vcb3.v[k4] ;
		    /* printf("c+%d c+%d c%d c%d  -> %lf ", i4,i3,i2,i1,val );
		       fockdisplay( &state4 ) ; */
                    if( findstate(nstates, state, &state4, &jst ) ) {
		      /* printf(" = state No %d\n", jst ) ; */
		      hamcol[jst] += val ;
		    }
		    else {
		      newstatecounter++ ;
		      /*
		      printf("c+%d c+%d c%d c%d  -> %lf ", i4,i3,i2,i1,val );
		      fockdisplay( &state4 ) ;		      
	              printf(" = new state under No %d\n", jst ) ;
		      */
		    } 
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    totnumele += darray2spamaline( &hamtransp[ist], nstates, hamcol ) ;
    diagsum += hamcol[ist] ;
  }
  printf("application of h produced %d new states\n", newstatecounter ) ;
  newstatecounter = 0 ;


  o1pdelete( h1p ) ;
  o2pdelete( vcb ) ;

  free( hamcol ) ; 

  printf("\n hamtransp. totnumele = %d. diagaverage = %lf \n\n",  
	 totnumele, diagsum / ( double ) nstates ) ;
  
  /*
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    printf("i =%3d --", ist ) ;   
    spamalinewrite( hamtransp[ist] ) ;
  }
  */
  
  totnumele = spamatranspose( &hamtransp, nstates, nstates ) ;
 
  printf("\n after spamatranspose. ham. totnumele = %d\n",  totnumele ) ;

  /*   
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    printf("i =%3d --", ist ) ;   
    spamalinewrite( hamtransp[ist] ) ;
  }
  */ 

  *pham = hamtransp ;

  return totnumele ;
}


/* 
    definition of functions declared at the top of the file 
*/
int isaconf( int nshells, int *lsh, int *occ ) {
  int ish ;
  for ( ish = 0 ; ish < nshells ; ish++ ) 
    if ( occ[ish] < 0 || occ[ish] > 4 * lsh[ish] + 2 ) 
      return 0 ;
  return 1 ;
}

int confinset( int nconfs, int nshells, int **occ, int *occdum ) {
  int iconf, ish, occallequal ;
  for ( iconf = 0 ; iconf < nconfs ; iconf++ ) {
    for ( occallequal = 1, ish = 0 ; occallequal && ish < nshells ; ish++ )
      occallequal = ( occdum[ish] == occ[iconf][ish] ? 1 : 0 ) ;
    if ( occallequal ) 
      return 1 ;    /* "occdum" is in set of confs defined by occ[][] */
  }
  return 0 ;
}

int cilistitemadd( struct CIlistitem **ip0, int j1, int j2, int j3, int j4,
		   int *lsh, int nr, double *rofi, double **ppofi, 
		   double rkscale ) {
  /* algorithm for non ordered list -- could be improved with ordered list */
  int nitems = 0 ;
  int n1, n2, nk, k ;
  int ten2k ;
  struct CIlistitem *ip ;

  /* quick return if all four shells are number 0, i.e. the core-shell -- 
     in this case we put the slater integrals Fk to zero   */
  if ( j1 == 0 && j2 == 0 && j3 == 0 && j4 == 0  )
    return 0 ;

  ip = *ip0 ;
  while( ip != 0 ) 
    if ( j1 == ip->ish1 && j2 == ip->ish2 && j3 == ip->ish3 && j4 == ip->ish4 )
      return 0 ;
    else { 
      ip = ip-> next ;
      ++nitems ;
    }
  ip = ( struct CIlistitem * ) malloc( sizeof( struct CIlistitem ) ) ;
  ip->ish1 = j1 ; ip->ish2 = j2 ; ip->ish3 = j3 ; ip->ish4 = j4 ;

  /* calculate nk = min(l4+l1,l2+l3) + 1.   ck(4,1)ck(2,3)!=0 for k=0..nk-1 */
  n1 = lsh[j4] + lsh[j1] ;
  n2 = lsh[j2] + lsh[j3] ;
  nk = ( n1 < n2 ? n1 : n2 ) + 1 ;
  ip -> nk = nk ;
  ip -> rmx = ( double * ) calloc( nk, sizeof( double ) ) ; 
  /* calculate  ip -> rmx[k] == rk(4,3;1,2)  */
  for ( k = 0 ; k < nk ; k++ ) { 
    ip -> rmx[k] = ( k == 0 ? 0. :
	slaterrk(nr,rofi,k,ppofi[j4],ppofi[j3],ppofi[j1],ppofi[j2]) ) ;
    ip -> rmx[k] *= rkscale ; 
  }
  ip -> next = *ip0 ;
  *ip0 = ip ;
  return ++nitems ;
}

void cilistdelete( struct CIlistitem *ip )
{
  int counter = 0 ;
  struct CIlistitem *ipn;
/* printf("Clearing CI list with cilistdelete(..) Cleared items (No,nk):\n");*/
  while ( ip != 0 ) {
    ipn = ip -> next ;
/*  printf(" (%d,%d)", ++counter, ip->nk ) ; */
    free( ip -> rmx ) ;
    free( ip ) ;
    ip = ipn ;
  }
/*  printf("\n") ; */
}
  
#undef HZ
