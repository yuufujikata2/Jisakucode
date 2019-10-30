#include<stdio.h>
#include <math.h>
#include <stdlib.h>
void main(){
    int i,nr=201;
    double *rofi,*v,a=0.03,b,wsr=2.,r;
    double V0 =-1.;
    b = wsr / ( exp(a*(nr-1)) -1. ) ;
    rofi = ( double * ) calloc( nr, sizeof( double ) ) ;
    v    = ( double * ) calloc( nr, sizeof( double ) ) ;

    for ( i = 0 ; i < nr ; i++ ) rofi[i] = b *( exp(a*i) -1. ) ;

    for ( i = 0 ; i < nr ; i++) {
    r = rofi[i];
    if ( r < wsr / sqrt(3)) {
      v[i] = V0;
    }
    else if ( r >= wsr / sqrt(3) && r < sqrt(2) * wsr / sqrt(3) ) {
      v[i] = V0 * ( 1 - 12 * M_PI * r * ( r - 0.5 * wsr )/(4 * M_PI * pow( r , 2. )));
    }
    else if ( r >= sqrt(2) * wsr / sqrt(3) && r <= wsr ){
      v[i] = V0 * ( r - 0.5 * sqrt(3) * wsr ) * ( 2 - 3 * sqrt(2) * ( sqrt(2) - 1 ) ) / ( pow( r , 2. ) * pow( sqrt(2) - sqrt(3) , 2.));
    }
    printf("%9.6lf  %9.6lf\n",rofi[i],v[i]);
  }
}

