#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSCO 1.e-5


double* makerofi(int nr, double rws, double* b, double a) ;

void interpolar(int nin, double *xin, double *yin,
                int nout, double *xout, double *yout) ;

double* readcofi(char* cowffilename, int iz, int nr, 
		 double a, double rws, int* ncorel, int* lcorel, 
		 double* ecore) {
/* reads corewavefct from external file (cowffilename) into array 
   cofi[i] (say), allocates memory and returns pointer to it
   also reads in ncorel, lcorel, ecore
   this version only for non-spin-polarized core-wave-fct
*/
    int i, c, izco, nrco, nris, nspco, ispco, dointerpolation ;
    char coreclass[10] ;
    double aco, rwsco, cofifac, *cofi, bco, b, *rco, *rpot, *newcofi ;
    FILE * fpcowf ; 

    fpcowf = fopen(cowffilename,"r");

    for ( i = 0 ; i < 3 ; i++ )
	while((c=fgetc(fpcowf))!=EOF&&c!='\n'); /*3 comment lines*/

    fscanf(fpcowf,"%s%d%d", coreclass, ncorel, lcorel);
    printf("%s\t%d\t%d\n", coreclass, *ncorel, *lcorel);

    for ( i = 0 ; i < 2 ; i++ )
	while((c=fgetc(fpcowf))!=EOF&&c!='\n'); /* comment line */

    fscanf(fpcowf,"%d%d%d%lf%lf",&izco,&nrco,&nspco,&aco,&rwsco);
    if ( nspco != 1 ) {
        printf("can't use this readcofi.c written for nspin = 1\n") ;
        exit(1);
    }
    if ( izco != iz || rwsco < rws - EPSCO ) {
        printf("FATAL inconsistency between core-wf-file and pot-file\n");
	printf("izco != iz || rwsco < rws + EPSCO \n") ;
	exit(1) ;
    }
    else if( nrco != nr || /* nspco != nspin || */ 
	     fabs(aco-a) > EPSCO || fabs(rwsco-rws) > EPSCO ) { 
        printf("WARNING: inconsistency between core-wf-mesh and pot-mesh\n");
        printf("Core-wf will be interpolated on pot-mesh\n") ;
/*        printf("    CHECK CODE !!!!\n") ; */
	dointerpolation = 1 ;
    }
    else {
        dointerpolation = 0 ;
    }

    for ( i = 0 ; i < 2 ; i++ )
	while((c=fgetc(fpcowf))!=EOF&&c!='\n'); /* comment line */

    fscanf(fpcowf,"%d%lf%d", &ispco, ecore, &nris );
    if ( ispco != 1 || nrco != nris ) {
	printf("trouble with ispco or nris in core-wf-file\n");
	exit(1) ;
    }
    cofi = ( double * ) malloc( 2 * nrco * sizeof( double ) ) ;
    for ( i = 0 ; i < nrco ; i++ ) {
	fscanf( fpcowf,"%lf%lf%lf",
		&cofi[i],&cofi[nr+i],&cofifac ) ;
	cofi[i] *= sqrt(cofifac) ;
    }
/*  without interpolation we had the return statement here:
    return cofi ;
*/
    
    if (dointerpolation) {
      rco  = makerofi(nrco, rwsco, &bco, aco) ;
      rpot = makerofi(nr, rws, &b, a) ;
      newcofi = ( double * ) malloc( 2 * nr * sizeof( double ) ) ;
      newcofi[0] = cofi[0]; 
      interpolar(nrco, rco, cofi, nr-1, rpot+1, newcofi+1) ;
      newcofi[nr] = cofi[nrco];
      interpolar(nrco, rco, cofi+nrco, nr-1, rpot+1, newcofi+nr+1) ;
      free(rpot); free(rco); free(cofi) ;
      return newcofi ;
    }
    else {
      return cofi ;
    }
}

#undef EPSCO
