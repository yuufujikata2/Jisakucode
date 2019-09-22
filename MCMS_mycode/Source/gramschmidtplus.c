#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1.e-10

double overlap( int nr, double *rofi, double *pofi1, double *pofi2 ) ;

int gramschmidtplus(int nr, double *rofi, int nrwf, double ***pppofi, 
		double ***psmat )
{
  int i, j, ir ;
  double **pqofi, **smat, sum, ovl, normpisq, norm ;
  FILE *fp ;

  /* allocation */
  pqofi = ( double ** ) malloc( nrwf * sizeof( double * ) ) ;
  for ( i = 0 ; i < nrwf ; i++ )
    pqofi[i] = ( double * ) calloc( 2 * nr , sizeof( double ) ) ;
  smat = ( double ** ) malloc( nrwf * sizeof( double * ) ) ;
  for ( i = 0 ; i < nrwf ; i++ )
    smat[i] = ( double * ) calloc( nrwf , sizeof( double ) ) ;

/* the first radial wf is the core wf. We assume it to be of different
   l value of the valence shell wfs, so no orthogonalization is needed */

/* make sure that core orbital is correctly normalized. If already done,
   one could replace the next 4 lines by the statement smat[0][0] = 0.;  */
  norm = sqrt( overlap( nr, rofi, (*pppofi)[0], (*pppofi)[0] ) ) ;
  smat[0][0] = norm ;
  for ( ir = 0 ; ir < 2*nr ; ir++ )
    pqofi[0][ir] = (*pppofi)[0][ir] / norm ;
  
  for ( i = 1 ; i < nrwf ; i++ ) {
    sum = 0. ;
    for ( j = 1 ; j < i ; j++ ) {
      ovl = overlap( nr, rofi, pqofi[j], (*pppofi)[i] );
      smat[j][i] = ovl ;
      sum += ovl * ovl ;
    }
    
    normpisq = overlap( nr, rofi, (*pppofi)[i], (*pppofi)[i] );

    norm = sqrt( normpisq - sum ) ;

    smat[i][i] = norm ;

    for ( ir = 0 ; ir < 2*nr ; ir++ ) {
      sum = 0. ;
      for ( j = 0 ; j < i ; j++ )  
	sum += pqofi[j][ir] * smat[j][i] ;
      pqofi[i][ir] = ( (*pppofi)[i][ir] - sum ) / norm ;
    }
  }
		  

  /* deallocation and moving of pointers */
  printf("deallocation of ppofi[i][j] for i<nrwf. should be changed to") ;
  printf("i<nrwfall (one day)\n") ;
  for ( i = 0 ; i < nrwf ; i++ ) {
    free( (*pppofi)[i] ) ;
    /*    (*pppofi)[i] = pqofi[i] ;
     */
  }
  free( (*pppofi) ) ;
  *pppofi = pqofi ;
  *psmat = smat ;


/* write wavefunctions into file "wavefcts" */
    fp = fopen("wavefcts2","w") ;
    fprintf(fp,"#r | p_i(r)\n");
    for ( ir = 0 ; ir < nr ; ir++ ) {
	fprintf(fp,"%15.8lf ", rofi[ir]) ;
	for ( i = 0 ; i < nrwf ; i++ )
	    fprintf(fp," %15.8e", (*pppofi)[i][ir]) ;
	fprintf(fp,"\n");
    }
    fclose(fp) ;


    printf("\n smat[i][j]\n\n") ;
    for ( i = 0 ; i < nrwf ; i++ ) {
      for ( j = 0 ; j < nrwf ; j++ ) 
	printf(" %10.6lf", smat[i][j] ) ;
      printf("\n") ;
    }

/*    printf("\n Check orthonormality: (pqofi[i]|pqofi[j])\n" ) ; */
    for ( i = 0 ; i < nrwf ; i++ )
      for ( j = 0 ; j < nrwf ; j++ ) {
	ovl = overlap(nr, rofi, (*pppofi)[i], (*pppofi)[j]) ;
	if ( i!=j && fabs(ovl) > EPSILON || i==j && fabs(ovl-1.) > EPSILON ) 
	  printf("i = %d j = %d (qi|qj) = %10.6lf\n", i, j, ovl );
      }

  return 1 ;
} 


