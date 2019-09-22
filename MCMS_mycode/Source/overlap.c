#include <stdlib.h>

double defintegr(double *f,  double *x, int Nx ) ;

double overlap( int nr, double *rofi, double *pofi1, double *pofi2 ) 
{
  int ir ;
  double *f, ovl ;
  
  f = ( double * ) malloc( nr * sizeof( double ) ) ;

  for ( ir = 0 ; ir < nr ; ir++ )
    f[ir] = pofi1[ir] * pofi2[ir]  +  pofi1[nr+ir] * pofi2[nr+ir] ;
 
  ovl = defintegr( f, rofi, nr ) ;
  
  free( f ) ;

  return ovl ;
}
