#include <math.h>
#define dmax1(a,b) ((a)>(b)?(a):(b))
void fctp(double a, double b, double e, int l, int *nctp, int nctp0,
	  int nr, double *rofi,double *v, int iz) {
/* finds classical turning point
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   nctp0 :minimum of effective potential
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci Inputs/Outputs:
Cio  nctp  :rofi(nctp) is classical turning point
C  --------------------------------------------------------------------- */
    int ir,irep,n1,n2,nlast,ntry;
    double r,veff,fllp1,dvdr,fofr,dfdr,rtry,zz,veffnrm1,veffnctp0;

    zz    = (double)(2*iz);
    fllp1 = (double)(l*(l+1));
    
    veffnrm1 = fllp1/(rofi[nr-1]*rofi[nr-1])-zz/rofi[nr-1]+v[nr-1];
    veffnctp0 = fllp1/(rofi[nctp0]*rofi[nctp0])-zz/rofi[nctp0]+v[nctp0];
    if (nctp0 == nr || e > veffnrm1)
	*nctp = nr-1;
    else if (e < veffnctp0) 
	*nctp = 1;
    else {
	n1 = nctp0;
	n2 = nr-2;
	nlast = -10; /*not changed*/
	for (irep = 0 ; irep < 20 ; irep++){
	    if (*nctp > n2 || *nctp < n1) *nctp = (n1 + n2)/2;
	    r = rofi[*nctp];
	    fofr=fllp1/(rofi[*nctp]*rofi[*nctp])-zz/rofi[*nctp]+v[*nctp] - e;
	    dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r +
		(v[*nctp+1] - v[*nctp-1])/(2.*a*(r + b));
	    rtry = dmax1(r-fofr/dfdr,rofi[1]);
	    ntry = log(rtry/b+1.)/a + 1.5;
	    if (nlast == *nctp){
		if (*nctp == nctp0+1) *nctp = 1;
		return ;
	    }
	    if (fofr > 0.)
		n2 = *nctp;
	    else
		n1 = *nctp;
	    nlast = *nctp ;
	    *nctp = ntry < (nr-2) ? ntry : (nr-2) ; /* min0(ntry,nr-1) */
	}
	if (*nctp == nctp0+1) *nctp = 1;
    }
}
