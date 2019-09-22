#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ndimarraycalloc.h"

#define LMAX_cont 3
#define D_cont   10
#define NP_cont 200
#define NB_cont   2
#define NSPIN 1
#define iabs(i) ( (i) > 0 ? (i) : -(i) )

struct Lmtopot {
  int z; double wsr; double a; double b; int nr; int nspin; double *v;
};

double  defintegr(double *f,  double *x, int Nx ) ;
double* readvofi(char* potfile, int* iz, int* nr, int* nspin,
		double* a, double* wsr) ;
struct Lmtopot * readrinopot( char *filename, int *nineqat ) ;
double* readcofi(char* cowffilename, int iz, int nr, 
	      double a, double wsr, int* ncorel, int* lcorel, 
		 double* ecore) ;
double* makerofi(int nr, double wsr, double* b, double a);
void rsq1(double a, double b, double e, double* g, int l, int nc,
	  int* nod, int nr, double* rofi, double* slonc, double* v, 
	  double* valnc, int iz) ;
void mkgfac(double e, double* gfac, int l, int nr, double* rofi, double* v, 
	    int iz) ;
void besneu(double x, int l, double* jl, double* nl, double* jlp, double* nlp);
/*   jl = j_l(x), nl = n_l(x), jlp = (d/dx) j_l(x), nlp = (d/dx) n_l(x) */


/*  void continuum_(float *avmtz, int *nep, float energyf[], 
		   float wnwjf[][D_cont][LMAX_cont+1], int *ntasc, 
		   float rhof[][NB_cont][NB_cont][2]; 
		   float chif[][NB_cont][NB_cont][2]; int *ylmisx) ;
 ATTENTION: 
   rhof, chif are actually complex, so in C we need two floats per element: 
   real part{ RHOF(N1,..,NX) } = rhof[n1]..[nx][0]
   imag.part{ RHOF(N1,..,NX) } = rhof[n1]..[nx][1]
*/
    
void invtau(int nep, int ntas, double ****tau) ;


void environrino(int nr, double *cofi, char* rinopotfile,
		 int nshell, double vmtzero, double imvmtzero, 
		 int nep, double *energy, int lsym, 
		 int ntas, double ****rho, double ****chi ) {
    int iz, nspin, ir, iep, iat, l, nod, lm1, lm2, ylmisx, ntasc ; 
    /*    int envfilestemlen, envfilelen ;  */
    int ish, is, i ;
    double a, b, wsr, r ; 
    double energy0, ke, bes, dbes, neu, dneu, pr0, dpr0, radipme0 ;
    double *rofi, *vofi, *pofi, *gfac, *fwork ;
    double zfac, rWjR, rWnR, wnwj0, sigma0, sigmaat, imtauavd ;
    float energyf[NP_cont], avmtz, aimvz ;
    float wnwjf[NP_cont][D_cont][LMAX_cont+1];
    float rhof[NP_cont][NB_cont][NB_cont][2];
    float chif[NP_cont][NB_cont][NB_cont][2];
    float radipmef[NP_cont][LMAX_cont+1];
    double ****tau;
    struct Lmtopot *potp ;

    FILE * fp ;
    
/* checks */
    if ( nep > NP_cont) {
      printf("nep > NP_cont -- increase NP_cont, NP_ to %d\n",nep ) ;
      exit(1) ;
    }
    if ( ntas > NB_cont ) {
      printf("nep > NB_cont -- increase NB_cont, NB_ to %d\n", ntas ) ;
      exit(1) ;
    }

/* clear Fortran array wnwjf   -   should be useless */
    for (iep = 0 ; iep < NP_cont ; iep++ ) 
      for ( ish = 0 ; ish < D_cont ; ish++ )
	for ( l = 0 ; l <= LMAX_cont ; l++ ) 
	  wnwjf[iep][ish][l] = 0. ;


/* read (again!) rinopotfile */
    potp = readrinopot(rinopotfile, &nshell ) ;
 

    fp = fopen("veevenofr","w");

/* loop over shells */
    for ( ish = 0 ; ish < nshell ; ish++ ) {
      iz    = potp[ish].z; 
      wsr   = potp[ish].wsr; 
      a     = potp[ish].a; 
      b     = potp[ish].b; 
      nr    = potp[ish].nr; 
      nspin = potp[ish].nspin ;
      printf("ish=%d z=%d wsr=%lf nr=%d nspin=%d\n", ish, iz, wsr, nr, nspin );

      rofi  = ( double * ) malloc( nr * sizeof( double ) ) ;
      vofi  = ( double * ) malloc( nr * sizeof( double ) ) ;
      fwork = ( double * ) calloc( nr, sizeof( double ) ) ;
      gfac  = ( double * ) calloc( nr, sizeof( double ) ) ;
      pofi  = ( double * ) calloc( 2 * nr, sizeof( double ) ) ;

      fprintf(fp,"# r_i, r_i*v(r_i) (v=H+XC-2z/r), shell no %d\n", ish ) ;

      for (ir = 0; ir < nr; ir++) {
/* make mesh */
	rofi[ir] = b * ( exp( a * ir ) - 1. ) ;	
/* copy and printout potential */
	vofi[ir] = potp[ish].v[ir] ;
	fprintf(fp,"%25.10E ", rofi[ir]) ;
	fprintf(fp,"%25.10E ", rofi[ir] * vofi[ir]        - 2.*iz ) ;
	fprintf(fp,"%25.10E ", rofi[ir] * potp[ish].v[ir] - 2.*iz ) ;
	fprintf(fp,"\n");
      }	

/* energy loop: for (energy0 = emin ; energy0 <= emax ; energy0 += deltae) */
      for (iep = 0 ; iep < nep ; iep++ ) {   
	energy0 = energy[iep] ;
	energyf[iep] = energy[iep] ;
	ke = sqrt( energy0 - vmtzero ) ;
/*    printf("energy=%lf, vmtzero=%lf, ke=%lf\n", energy0, vmtzero, ke ) ; */
	for ( l = 0 ; l <= LMAX_cont ; l++ ) {
	  besneu( ke*wsr, l, &bes, &neu, &dbes, &dneu );
	  dbes *= ke ; dneu *= ke ; 
/* dbes := dj/dr = dx/dr * dj/dx = ke*dbes */
/*	printf("ke*wsr=%lf, bes=%lf, dbes=%lf, neu=%lf, dneu=%lf\n", 
	        ke*wsr, bes, dbes, neu, dneu );  */
	  rsq1(a,b,energy0,pofi,l,nr-1,&nod,nr,rofi,&dpr0,vofi,&pr0,iz);
/* calculate wnwj = -(kK)^{-1} = W[neu,P/r]/W[bes,P/r] */
	  rWjR = bes * dpr0 - dbes * pr0 - bes * pr0 / wsr ;
	  rWnR = neu * dpr0 - dneu * pr0 - neu * pr0 / wsr ;
	  wnwjf[iep][ish][l] = rWnR / rWjR ;
	  if ( ish == 0 ) {  /* calculate radial dipole matrix element */
	    wnwj0 = rWnR / rWjR ;
	    mkgfac(energy0, gfac, l, nr, rofi, vofi, iz) ; 
	    for ( ir = 0 ; ir < nr ; ir++ ) {
	      pofi[ir] *= sqrt(gfac[ir]) ;
	      fwork[ir] = 
		( cofi[ir] * pofi[ir] + cofi[nr+ir] * pofi[nr+ir] ) * rofi[ir];
	    }
	    radipme0 = defintegr(fwork,rofi,nr) ;
/*	printf("pr0=%lf, dpr0=%lf,radipme0=%lf\n", pr0, dpr0, radipme0);*/
/* calculate ramsauer factor W[bes,P/r] * k * r * r */
/*	ramf = rWjR * ke * wsr ;                    */  
/* my way:            zfac/r0 = j_l * t^{-1} -/+ i * k  * h(+/-)_l  
   is correct amplitude of radial wavefunction at r0 
*/ 
	    zfac = wsr * ke * ( - bes * wnwj0 + neu ) ;
	    radipme0 *= zfac / pr0 ;             /* T-matrix normalization */
	    radipmef[iep][l] = radipme0 ;
	  }
	}
      }
      free(fwork); free(gfac); free(pofi); free(vofi); free(rofi);
      free(potp[ish].v) ;
    }
    free(potp) ;

    fclose(fp) ;

/* call continuum */
    avmtz = (float) vmtzero ;
    aimvz = (float) imvmtzero ;
    ylmisx = 0 ;        /* ylmisx==1  means _ylm_ _is_ comple_x_ */
    if (!ylmisx) printf("ATTENTION: real Ylm -- NO trafo to complex Ylm\n") ;  
    continuum_(&avmtz,&aimvz,&nep,energyf,wnwjf,&ntasc,rhof,chif,&ylmisx) ;

    printf("CONTINUUM successfully completed.\n") ;
    printf("Reflectivity expressed in ") ;
    if (ylmisx) printf("complex ") ; else printf("real ") ;
    printf("spherical harm. basis\n") ;

    printf(" #basis states on absorber ntasc = %d",ntasc ) ;
    if ( ntasc != ntas ) {
      printf(" != ntas = %d.\n", ntas ) ;
      exit(1) ;
    }
    printf("\n");
    

 /* copy rho and chi and write to file "rhochi.dat" */
    fp = fopen("rhochi.dat","w") ;
    fprintf(fp,"# comments\n");
    fprintf(fp,"# comments\n");
    fprintf(fp,"# 1.L: neenv, ntasc, nnzele. Rest: energy\\n ");
    fprintf(fp,"lm1 lm2 Rerho Imrho Rechi Imchi\n") ;
    fprintf(fp,"%3d  %3d  %6d\n", nep, ntas, ntas*ntas); /*nnzele 2b revised*/
    for ( iep = 0 ; iep < nep; iep++ ) {
      fprintf(fp,"%17.10lg\n", energy[iep]);
      ke = sqrt(energy[iep] - vmtzero) ;
      for ( lm1 = 0 ; lm1 < ntas; lm1++ ) 
	for ( lm2 = 0 ; lm2 < ntas; lm2++ ) {
	  for ( i = 0 ; i < 2 ; i++ ) {
	    rho[lm1][lm2][i][iep] = - ke * (double) rhof[iep][lm1][lm2][i] ;
	    chi[lm1][lm2][i][iep] = - (double) chif[iep][lm1][lm2][i] ;
	    /* using my definitions with tau^{-1} rather than (tau*k)^{-1} */
	    /* might do transposition lm1<->lm2 needed because of different 
	       array orders in C and Fortran here rather than in continuum.*/
	  }
	  fprintf(fp,"%3d %3d %17.10lg %17.10lg %17.10lg %17.10lg\n",
		  lm1,lm2, 
		  rho[lm1][lm2][0][iep] ,rho[lm1][lm2][1][iep],
		  chi[lm1][lm2][0][iep] ,chi[lm1][lm2][1][iep] );
	}
    }
    fclose(fp);

/* define tau^{-1} = - k *( wnwj + i - rho ) , but call it tau */
    tau = calloc4double( nep, ntas, ntas, 2 ) ;

    printf(" lsym = %d, ntas = %d\n", lsym, ntas) ;
    if ( lsym < 0 && ntas >= 9 ) { 
      for ( iep=0; iep<nep; iep++) {
	ke = sqrt(energy[iep] - vmtzero) ;
	for (lm1=0; lm1<ntas; lm1++) {
	  l = (int)sqrt(lm1+0.1) ; 
	  for (lm2=0; lm2<ntas; lm2++) {
	    tau[iep][lm1][lm2][0] = 
	      ( lm1==lm2 ? -ke * wnwjf[iep][0][l] : 0. )
	      - rho[lm1][lm2][0][iep] ;
	    tau[iep][lm1][lm2][1] =
	      ( lm1==lm2 ? -ke : 0. )
	      - rho[lm1][lm2][1][iep] ;
	  }
	}
      }    
    }
    else if ( lsym >= 0 && ntas == 1 ) {
      for ( iep=0; iep<nep; iep++) {
	ke = sqrt(energy[iep] - vmtzero) ;
	tau[iep][0][0][0] =  -ke * wnwjf[iep][0][lsym] - rho[0][0][0][iep] ;
	tau[iep][0][0][1] =  -ke - rho[0][0][1][iep] ;
      }
    }  
    else {
      printf("Don't know how to define tau in single channel calc.\n");
      return ;
    } 

    invtau(nep, ntas, tau) ;

/*  calculate sigma in spherical approx to tau and neglecting s,p,f */
    fp = fopen("sigma1ch.dat","w") ;
    fprintf(fp,"#    energy   sigma_d(cluster) sigma_d(atom) Im_tau_av_d\n") ; 
/* non-symmetrized basis => assume normal lm-order; make trace over d-states*/
    if ( lsym < 0 && ntas >= 9 ) { 
      for ( iep=0; iep<nep; iep++) {
	imtauavd = 0. ;
	for ( lm1 = 4; lm1 < 9; lm1++ )
	  imtauavd += tau[iep][lm1][lm1][1];
	imtauavd /= 5. ;
	ke = sqrt( energy[iep] - vmtzero ) ;
	radipme0 = radipmef[iep][2] ;
	wnwj0 = wnwjf[iep][0][2] ;
	sigma0 = radipme0 * radipme0 * imtauavd ;
	sigmaat = radipme0 * radipme0 / ( ke * ( wnwj0 * wnwj0 + 1.) ) ;
	fprintf(fp,"%15.6lg %15.6lg %15.6lg %15.6lg\n", 
		energy[iep], sigma0, sigmaat, imtauavd );
	/*
	fprintf(fp,"%15.6lg %15.6lg %15.6lg\n", energy[iep], sigma0, sigmaat );
	*/
      }
    }
/* symmetrized basis and only 1 basis state on absorber with l=lsym */
    else if ( lsym >= 0 && ntas == 1 ) {
      for ( iep=0; iep<nep; iep++) {
	imtauavd = tau[iep][0][0][1];
	ke = sqrt( energy[iep] - vmtzero ) ;
	radipme0 = radipmef[iep][lsym] ;
	wnwj0 = wnwjf[iep][0][lsym] ;
	sigma0 = radipme0 * radipme0 * imtauavd ;
	sigmaat = radipme0 * radipme0 / ( ke * ( wnwj0 * wnwj0 + 1.) ) ;
	fprintf(fp,"%15.6lg %15.6lg %15.6lg %15.6lg\n", 
		energy[iep], sigma0, sigmaat, imtauavd );
      }     
    }
    else {
      fprintf(fp,"Don't know how to calc. 1ch. spec\n"); 
    }
    fclose(fp);

    printf("Single channel spectrum (neglecting p->s transition) written ") ;
    printf("into file %s.\n End of environment.\n", "sigma1ch.dat" ) ;

    free4double(nep,ntas,ntas,tau); 
}


