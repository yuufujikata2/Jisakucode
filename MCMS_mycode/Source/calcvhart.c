#include <stdlib.h>
void integrate(double *f, double *fint, double *x, int Nx ) ;

void calcvhart(int nr, double* rofi, double* pofi, double* vhart) {
/*
 *   calculates the bare (=unscreened) HARTREE POTENTIAL DUE TO the charge 
 *   density of ONE ELECTRON in a spherically symmetric orbital pofi 
 *
 *   vhart(r) = 2 (1/r) \int_0^r dr' P(r')^2 + 2 \int_r^{r0} dr' P(r')^2 / r'
 *
 *   with P(rofi[i]) = pofi[i] and r0 = rofi[nr-1] is the atomic radius
 *   which is also the cut off of the second integral
 */ 
    int i;
    double dum, psqr_over_r_defint ;
    double *psqr, *psqr_over_r, *psqr_int, *psqr_over_r_int ;

    psqr            = ( double * ) malloc( nr * sizeof( double ) ) ;
    psqr_over_r     = ( double * ) malloc( nr * sizeof( double ) ) ;
    psqr_int        = ( double * ) malloc( nr * sizeof( double ) ) ;
    psqr_over_r_int = ( double * ) malloc( nr * sizeof( double ) ) ;
    
    psqr[0] = 0.;
    psqr_over_r[0] = 0.;

    for (i=1;i<nr;i++) {
      dum = pofi[i] * pofi[i] + pofi[nr+i] * pofi[nr+i] ;
      psqr[i]        = dum ;
      psqr_over_r[i] = dum / rofi[i] ;
    }

    integrate(psqr,psqr_int,rofi,nr);
    integrate(psqr_over_r,psqr_over_r_int,rofi,nr);
    
    psqr_over_r_defint = psqr_over_r_int[nr-1] ;

    vhart[0] = psqr_over_r_defint ;
    for (i=1;i<nr;i++) {
      vhart[i] = 2. * ( psqr_int[i] / rofi[i] 
	              + psqr_over_r_defint -  psqr_over_r_int[i] ) ;
    }

    free(psqr_over_r_int) ;
    free(psqr_int) ;    
    free(psqr_over_r) ;
    free(psqr) ;
}
