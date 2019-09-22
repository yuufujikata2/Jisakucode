#include <stdlib.h>
#include <math.h>
double* makerofi(int nr, double rws, double* b, double a){
/* rofi(i) = b [e^(a*i) -1]; note: rofi[0]=0, nr,rws and a on input,
   b on output such that rofi[nr-1]=rws; i.e. b=rws/(e^(a*(nr-1))-1)
*/
    int i;
    double * rofi;
    *b = rws/(exp(a*(nr-1))-1.);
    rofi = ( double * ) calloc( nr, sizeof( double ) ) ;
    for (i=0; i<nr; i++) *(rofi+i) = (*b)*(exp(a*i)-1.);
    return rofi;
}



