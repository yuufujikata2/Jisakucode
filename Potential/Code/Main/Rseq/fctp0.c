#define IRMIN 10
/*  parameter (irmin=11) */
void fctp0(int l, int *nctp0, int nr, double *rofi, double *v, int iz) {
/* Initialize things for FCTP, which finds classical turning point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :angular momentum
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci   nr    :number of mesh points
Co Outputs:
Co   nctp0 :minimum of effective potential
Co          or nr if v(nr) > vmin + 3 Ry   //pk: sign in contradict. to below!
C ----------------------------------------------------------------------*/
/*Local variables: */
    int ir,irmin;
    double veff,zz,fllp1,vi,vim1,vnctp0,vnrm1 ;

    zz = (double)(2*iz);
    fllp1 = (double)(l*(l+1));
    
    ir=IRMIN;
    vi=fllp1/(rofi[ir]*rofi[ir])-zz/rofi[ir]+v[ir];
    vim1=fllp1/(rofi[ir-1]*rofi[ir-1])-zz/rofi[ir-1]+v[ir-1];
    while (vi <= vim1 && ir < nr-1) {
	ir=ir+1;
	vim1=vi;
	vi=fllp1/(rofi[ir]*rofi[ir])-zz/rofi[ir]+v[ir];
    }
    
    *nctp0 = ir-1;
    vnctp0 = fllp1/(rofi[ir-1]*rofi[ir-1])-zz/rofi[ir-1]+v[ir-1];
    vnrm1  = fllp1/(rofi[nr-1]*rofi[nr-1])-zz/rofi[nr-1]+v[nr-1];
    if (vnctp0 >= vnrm1 - 3.) 
	*nctp0 = nr-1;
    return;
}
