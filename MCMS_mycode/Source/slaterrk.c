#include <stdlib.h>
void integrate(double *f, double *fint, double *x, int Nx ) ;
double defintegr(double *f,  double *x, int Nx ) ;

double xton(double x, int n){
/* returns x^n for integer n >=0 only */
    double p = 1.; 
    while(n>0){ p *= x; n--; };
    return p;
} 

double slaterrk(int nr, double* rofi, int k, double* pi, double* pj, 
                double* pt, double* pu){
/*   calculates the slater integral R^k(ij,tu), see Cowan (6.23)
 *   \int_r1 \int_r2 2*(rs^k/rg^{k+1}) pi(r1)*pj(r2)*pt(r1)*pu(r2)
 *   where rs = min(r1,r2) and rg = max(r1,r2)                     
 */ 
    int i;
    double r, rtok, rtomkm1, pipt, f2sum, pjpu, result;
    double *f1, *f2, *f1int, *f2int;

    f1 = ( double * ) malloc( 4 * nr * sizeof( double ) ) ;
    f2 = f1+nr;
    f1int = f2+nr;
    f2int = f1int+nr;

    f1[0] = 0.;
    f2[0] = 0.;
    for (i=1;i<nr;i++) {
        r    = rofi[i];
        rtok = xton(r,k);
        rtomkm1 = 1./rtok/r;
        pipt = pi[i]*pt[i]+pi[nr+i]*pt[nr+i];  
        f1[i] = pipt * rtok;
        f2[i] = pipt * rtomkm1;
    }
    integrate(f1,f1int,rofi,nr);
    integrate(f2,f2int,rofi,nr);
    f1[0] = 0.;
    f2sum = f2int[nr-1];
    for (i=1;i<nr;i++) {
        r    = rofi[i];
        rtok = xton(r,k);
        rtomkm1 = 1./rtok/r;
        pjpu = pj[i]*pu[i]+pj[nr+i]*pu[nr+i];
        f1[i] = ( rtomkm1 * f1int[i] + rtok * (f2sum - f2int[i]) ) * pjpu;
    }
    result = 2.*defintegr(f1,rofi,nr);

    free(f1) ;
    
    return result; 
}
