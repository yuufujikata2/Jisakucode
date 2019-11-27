#include <stdio.h>
#include <stdlib.h>
void gintsr(double a, double b, double* ga, double* gb, double* gfac,
	    int nr, double* rofi, double* sum){
/*
C- Integrate inner product of two wave equations
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   ga    :First wave function (times r)
Ci   gb    :Second wave function (times r)
Ci   gfac  :normalization factor for first component of g
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Co Outputs:
Co   sum   :inner product
Cr Remarks:
Cr   Uses Simpson's rule
Cr   Attention: an integration from 0 to nre cannot be done by putting
Cr   nr=nre
C ----------------------------------------------------------------------
*/
/*  Local variables: */
    double fi,wgt,dum;
    int ir,jr;

    if (nr%2 != 1) {
	printf("FATAL ERROR in gintersr: nr must be odd");
	exit(1) ;
    }
    dum = 0.;
    for (ir = 1; ir < nr-1; ir += 2){
/*       do  ir = 2, nr-1, 2   */
       jr = ir + nr;
       fi = gfac[ir]*ga[ir]*gb[ir] + ga[jr]*gb[jr];
       wgt = rofi[ir] + b;
       dum = dum + wgt * fi;
    }
    dum = dum + dum;

    for (ir = 2; ir < nr-2; ir += 2){
/*   do  ir = 3, nr-2, 2   */
        jr = ir + nr;
        fi = gfac[ir]*ga[ir]*gb[ir] + ga[jr]*gb[jr];
        wgt = rofi[ir] + b;
        dum = dum + wgt * fi;
    }
    ir = 0;
    jr = ir + nr;
    fi = gfac[ir]*ga[ir]*gb[ir] + ga[jr]*gb[jr];  
    wgt = (rofi[ir] + b)*0.5;
    dum = dum + wgt * fi;

    ir = nr-1;
    jr = ir + nr;
    fi = gfac[ir]*ga[ir]*gb[ir] + ga[jr]*gb[jr];
    wgt = (rofi[ir] + b)*0.5;
    dum = dum + wgt * fi;

    *sum = dum*a*2./3.;
}
