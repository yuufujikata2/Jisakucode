#include <stdlib.h>

int invdcmplxmat( int matdim , double* mat ) ;

void invtau(int nep, int ntas, double ****tau) {
  int matdim, iep, lm1, lm2, ire;
  double *mat;

  mat = ( double * ) malloc( 2 * ntas * ntas * sizeof( double ) );/* 2=Re,Im */
  for ( iep = 0 ; iep < nep ; iep++ ){
    for ( lm1 = 0 ; lm1 < ntas ; lm1++ ) {
      for ( lm2 = 0 ; lm2 < ntas ; lm2++ ) {
	ire = 2 *( lm1 + ntas * lm2 ) ;
	mat[ire  ] = tau[iep][lm1][lm2][0] ;
	mat[ire+1] = tau[iep][lm1][lm2][1] ;
      }
    }
    invdcmplxmat( ntas , mat ) ;
    for ( lm1 = 0 ; lm1 < ntas ; lm1++ ) {
      for ( lm2 = 0 ; lm2 < ntas ; lm2++ ) {
	ire = 2 *( lm1 + ntas * lm2 ) ;
	tau[iep][lm1][lm2][0] = mat[ire  ] ;
	tau[iep][lm1][lm2][1] = mat[ire+1] ;
      }
    }    
  }
  free(mat);
}
    /* another possiblity for innermost loop :
	for ( ic = 0 ; ic < 2 ; ic++ ) {
	  imat = ic + 2 * ( lm1 + ntas * lm2 ) ;
	  mat[imat] = tau[iep][lm1][lm2][ic] ;
        }
    */
