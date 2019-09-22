#include <stdlib.h>

void reorderrwfs( double *cofi, double ecore, int *pnrwf, int nopen1, 
		  double ***pppofi, double **pe0p, double **ppr0, 
		  double **pdpr0 ) 
{
  int nrwf, i ;
  double *e0p, *pr0, *dpr0, **ppofi ;
  
  nrwf = (*pnrwf) + 1 ;

  e0p  = ( double * ) malloc( nrwf * sizeof( double ) ) ;
  pr0  = ( double * ) malloc( nrwf * sizeof( double ) ) ;
  dpr0 = ( double * ) malloc( nrwf * sizeof( double ) ) ;
  ppofi = ( double ** ) malloc( nrwf * sizeof( double * ) ) ;

  e0p[0]   = ecore ;
  pr0[0]   = 0. ;
  dpr0[0]  = 0. ;
  ppofi[0] = cofi ;

  e0p[1]   = (*pe0p)[nopen1] ; 
  pr0[1]   = (*ppr0)[nopen1] ; 
  dpr0[1]  = (*pdpr0)[nopen1] ;
  ppofi[1] = (*pppofi)[nopen1] ;

  for ( i = 0 ; i < nopen1 ; i++ ) {
    e0p[i+2]   = (*pe0p)[i] ; 
    pr0[i+2]   = (*ppr0)[i] ; 
    dpr0[i+2]  = (*pdpr0)[i] ;
    ppofi[i+2] = (*pppofi)[i] ;
  }
  for ( i = nopen1 + 1 ;  i < (*pnrwf) ; i++ ) {
    e0p[i+1]   = (*pe0p)[i] ; 
    pr0[i+1]   = (*ppr0)[i] ; 
    dpr0[i+1]  = (*pdpr0)[i] ;
    ppofi[i+1] = (*pppofi)[i] ;
  }

  free(*pe0p) ;
  free(*ppr0) ;
  free(*pdpr0) ;
  free(*pppofi) ;

  *pe0p   = e0p ; 
  *ppr0   = pr0 ; 
  *pdpr0  = dpr0 ;
  *pppofi = ppofi ;

  *pnrwf = nrwf ;
}
