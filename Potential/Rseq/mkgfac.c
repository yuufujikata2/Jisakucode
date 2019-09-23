#define CSL 274.071979 
void mkgfac(double e, double* gfac, int l, int nr, double* rofi, double* v, int iz){
/*
C-Makes normalization factor for first component of g
C ----------------------------------------------------------------------
Ci Inputs:
Ci   e     :Energy
Ci   l     :angular momentum
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Co Outputs:
Co   gfac  :normalization factor for first component of g
*/

    int ir;
    double r,tmc,tmcr,fi,zz,fllp1;

    fllp1 = l*(l+1);
    zz = (double) (2*iz); 
    for (ir = 1; ir < nr; ir++){
/*      do  ir = 2, nr    */
        r = rofi[ir];
	tmc = CSL + (e + zz/r - v[ir])/CSL;
        tmcr= tmc*r;
        gfac[ir]= 1. + fllp1/tmcr/tmcr;
    }
}
