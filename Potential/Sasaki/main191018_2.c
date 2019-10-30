#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

/*
#include "w3j.h"
#include "ndimarraycalloc.h"
#include "globals.h"
struct Lmtopot {
  int z; double wsr; double a; double b; int nr; int nspin; double *v;
};
struct Lmtopot readlmtopot(char *filename) ;
double overlap( int nr, double *rofi, double *pofi1, double *pofi2 ) ;
*/

double  defintegr(double *f,  double *x, int Nx ) ;
double complex cylm( int l, int m, double x, double phi);
void sp_sym_3dinterpolar(int nin, double *xin, double *yin, int nout, double ***xout, double ***yout);

#define V0 -1.

void rseq(double a, double b, double *e, double eb1, double eb2,
	  double *g, double *gfac, int l, int nod, int nr, int *nre,
	  double *rofi, double *slo, double *v, double *val, int iz) ;
void definepotential( int nr, double *rofi, double *v , double wsr, int gridn , double ***cellpot , double *gridx) {
  int i, j, k ;
  double r , range;

  // uniform cubic potential spherical part (v[i])
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
  }
  
  // before devide potential(cellpot[x][y][z])
  range = wsr / sqrt(3);
  for ( i = 0 ; i < gridn ; i++) {
    for ( j = 0 ; j < gridn ; j++) {
      for ( k = 0 ; k < gridn ; k++) {
	if(fabs(gridx[i]) <= range && fabs(gridx[j]) <= range && fabs(gridx[k]) <= range){
	  cellpot[i][j][k] = V0;
	}
      } 
    } 
  }
  
}

void main()
{
  const double a = 0.03 ;
  int l, nod, nr, nre, iz, i, j, k ;
  int l1, l2, n1, n2, m1, m2, x_u, y_u, z_u;
  int lmax, cln, opn, qmax, nodn, gridn;
  double sgnf, pol_phi, delta_v; 
  double wsr, b, e, emin, emax, *g, *gfac, *rofi, slo, *v, val ;   
  double ****qmat, ****lmat, ****hsmat, ***allg, **allslo, **allval, **alle, *hsg;
  double complex *****ally, ******umat;
  double *gridx, ***gridr, ***gridv, ***pol_r, ***ipol_v, ***upot, ***ipol_g, *****ipol_allg, ***cellpot;
  clock_t time_start, time_end;
  FILE *sf;

  sf = fopen("test2.dat","w");
  if(sf==NULL){
    printf("file open error\n");
    exit(0);
  }
  
  //  time_start = clock();
  
  for(gridn=21;gridn<180;gridn++){
    time_start = clock();
    printf("%d\t",gridn);
    
    nr = 201 ;
    // gridn = 171;
    wsr = 2. ;
  
    // match grid and x_n of nmat
    gridx = ( double * ) calloc( gridn, sizeof( double )) ;
  
    for (x_u = 0 ; x_u < gridn ; x_u++ ){
      gridx[x_u] = ((double)x_u - ((double)gridn - 1.) / 2.) * 2. * wsr / (double)gridn;      
    }
    /*  for(i=0;i<gridn;i++){
	printf("%lf\t",gridx[i]);
	}*/

    //  printf("\n");
    // spherical potential and before devided potential calculate
    rofi = ( double * ) calloc( nr, sizeof( double ) ) ;
    v    = ( double * ) calloc( nr, sizeof( double ) ) ;
    gfac = ( double * ) calloc( nr, sizeof( double ) ) ;
    g    = ( double * ) calloc( 2 * nr, sizeof( double ) ) ;
  
    cellpot = ( double *** ) calloc( gridn, sizeof( double ** ));
    for ( x_u = 0; x_u < gridn ; x_u++ ){
      cellpot[x_u]= ( double ** ) calloc( gridn, sizeof( double* ));
      for ( y_u = 0; y_u < gridn ; y_u++ ){
	cellpot[x_u][y_u]= ( double * ) calloc( gridn, sizeof( double ));
      } 
    }
  
    b = wsr / ( exp(a*(nr-1)) -1. ) ;
    for ( i = 0 ; i < nr ; i++ ) rofi[i] = b *( exp(a*i) -1. ) ;
    definepotential( nr, rofi, v, wsr, gridn, cellpot, gridx ) ;

    /*  for (x_u = 0 ; x_u < gridn ; x_u++ ){
	printf("\n");
	for (y_u = 0 ; y_u < gridn ; y_u++ ){
	printf("\n");
	for (z_u = 0 ; z_u < gridn ; z_u++ ){
	printf("%f", cellpot[x_u][y_u][z_u]);
	}      
	}      
	}*/
    
    // setting
    lmax = 4 ; cln = 2 ; opn = 1 ; emin = V0; emax = 1000. ; e = 0. ;
    nodn = cln + opn;
    qmax = lmax * nodn;
    delta_v = pow( 2. * wsr / ( (double)gridn ) , 3. );
  
    /*  printf("\nl max = %d (0 ~ %d)",lmax-1,lmax-1);
	printf("\nclose channel = %d",cln);
	printf("\nopen channel = %d",opn);
	printf("\nenergy min = %lf",emin);
	printf("\nenergy max = %lf\n",emax);
    */
    qmat = ( double **** ) calloc( lmax, sizeof( double *** ));
    for ( l1 = 0; l1 < lmax ; l1++ ){
      qmat[l1]= ( double *** ) calloc( lmax, sizeof( double **));
      for ( l2 = 0; l2 < lmax ; l2++ ){
	qmat[l1][l2]= ( double ** ) calloc( nodn, sizeof( double *));
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  qmat[l1][l2][n1]= ( double * ) calloc( nodn, sizeof( double ));
	}
      }
    }
  
    lmat = ( double **** ) calloc( lmax, sizeof( double *** ));
    for ( l1 = 0; l1 < lmax ; l1++ ){
      lmat[l1]= ( double *** ) calloc( lmax, sizeof( double **));
      for ( l2 = 0; l2 < lmax ; l2++ ){
	lmat[l1][l2]= ( double ** ) calloc( nodn, sizeof( double *));
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  lmat[l1][l2][n1]= ( double * ) calloc( nodn, sizeof( double ));
	}
      }
    }
  
    hsmat = ( double **** ) calloc( lmax, sizeof( double *** ));
    for ( l1 = 0; l1 < lmax ; l1++ ){
      hsmat[l1]= ( double *** ) calloc( lmax, sizeof( double **));
      for ( l2 = 0; l2 < lmax ; l2++ ){
	hsmat[l1][l2]= ( double ** ) calloc( nodn, sizeof( double *));
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  hsmat[l1][l2][n1]= ( double * ) calloc( nodn, sizeof( double ));
	}
      }
    }
  
    umat = ( double complex ****** ) calloc( lmax, sizeof( double complex ***** ));
    for ( l1 = 0; l1 < lmax ; l1++ ){
      umat[l1]= ( double complex ***** ) calloc( lmax, sizeof( double complex ****));
      for ( l2 = 0; l2 < lmax ; l2++ ){
	umat[l1][l2]= ( double complex **** ) calloc( nodn, sizeof( double complex ***));
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  umat[l1][l2][n1]= ( double complex *** ) calloc( nodn, sizeof( double complex ** ));
	  for ( n2 = 0; n2 < nodn ; n2++ ){
	    umat[l1][l2][n1][n2]= ( double complex ** ) calloc( l1 * 2 + 1, sizeof( double complex * ));
	    for ( m1 = 0; m1 < l1 * 2 + 1 ; m1++ ){
	      umat[l1][l2][n1][n2][m1]= ( double complex * ) calloc( l2 * 2 + 1, sizeof( double complex ));
	    }
	  }
	}
      }
    }
  
    allslo = ( double ** ) calloc( lmax, sizeof( double * ));
    for ( l1 = 0; l1 < lmax; l1++ ){
      allslo[l1]= ( double * ) calloc( nodn, sizeof( double ));
    }
  
    allval = ( double ** ) calloc( lmax, sizeof( double * ));
    for ( l1 = 0; l1 < lmax; l1++ ){
      allval[l1]= ( double * ) calloc( nodn, sizeof( double ));
    }
  
    alle = ( double ** ) calloc( lmax, sizeof( double * ));
    for ( l1 = 0; l1 < lmax; l1++ ){
      alle[l1]= ( double * ) calloc( nodn, sizeof( double ));
    }
  
    allg = ( double *** ) calloc( lmax, sizeof( double ** ));
    for ( l1 = 0; l1 < lmax ; l1++ ){
      allg[l1]= ( double ** ) calloc( nodn, sizeof( double* ));
      for ( n1 = 0; n1 < nodn ; n1++ ){
	allg[l1][n1]= ( double * ) calloc( nr, sizeof( double ));
      } 
    }
  
    hsg = ( double * ) calloc( nr, sizeof( double ) ) ;
  
    ally = ( double complex ***** ) calloc( lmax, sizeof( double complex **** ));
    for ( l1 = 0; l1 < lmax; l1++ ){
      ally[l1]= ( double complex **** ) calloc( l1 * 2 + 1, sizeof( double complex *** ));
      for ( m1 = 0; m1 < l1 * 2 + 1; m1++ ){
	ally[l1][m1]= ( double complex *** ) calloc( gridn, sizeof( double complex ** ));
	for ( x_u = 0; x_u < gridn; x_u++ ){
	  ally[l1][m1][x_u]= ( double complex ** ) calloc( gridn, sizeof( double complex * ));
	  for ( y_u = 0; y_u < gridn; y_u++ ){
	    ally[l1][m1][x_u][y_u]= ( double complex * ) calloc( gridn, sizeof( double complex ));
	  } 
	} 
      }
    }
  
    ipol_allg = ( double ***** ) calloc( lmax, sizeof( double **** ));
    for ( l1 = 0; l1 < lmax; l1++ ){
      ipol_allg[l1]= ( double **** ) calloc( nodn, sizeof( double *** ));
      for ( n1 = 0; n1 < nodn; n1++ ){
	ipol_allg[l1][n1]= ( double *** ) calloc( gridn, sizeof( double ** ));
	for ( x_u = 0; x_u < gridn; x_u++ ){
	  ipol_allg[l1][n1][x_u]= ( double ** ) calloc( gridn, sizeof( double * ));
	  for ( y_u = 0; y_u < gridn; y_u++ ){
	    ipol_allg[l1][n1][x_u][y_u]= ( double * ) calloc( gridn, sizeof( double ));
	  } 
	} 
      }
    }
  
    // interpolate spherical potential
    pol_r = ( double *** ) calloc( gridn, sizeof( double ** ));
    ipol_v = ( double *** ) calloc( gridn, sizeof( double ** ));
    upot = ( double *** ) calloc( gridn, sizeof( double ** ));
    ipol_g = ( double *** ) calloc( gridn, sizeof( double ** ));
    for ( x_u = 0; x_u < gridn ; x_u++ ){
      pol_r[x_u] = ( double ** ) calloc( gridn, sizeof( double* ));
      ipol_v[x_u] = ( double ** ) calloc( gridn, sizeof( double * ));
      upot[x_u] = ( double ** ) calloc( gridn, sizeof( double * ));
      ipol_g[x_u] = ( double ** ) calloc( gridn, sizeof( double * ));
      for ( y_u = 0; y_u < gridn ; y_u++ ){
	pol_r[x_u][y_u]= ( double * ) calloc( gridn, sizeof( double ));
	ipol_v[x_u][y_u] = ( double * ) calloc( gridn, sizeof( double ));
	upot[x_u][y_u] = ( double * ) calloc( gridn, sizeof( double ));
	ipol_g[x_u][y_u] = ( double * ) calloc( gridn, sizeof( double ));
      } 
    }
  
    for (x_u = 0 ; x_u < gridn ; x_u++ ){
      for (y_u = 0 ; y_u < gridn ; y_u++ ){
	for (z_u = 0 ; z_u < gridn ; z_u++ ){
	  pol_r[x_u][y_u][z_u] = sqrt(pow(gridx[x_u] , 2.) + pow(gridx[y_u] , 2.) + pow(gridx[z_u] , 2.));
	}      
      }      
    }
  
    sp_sym_3dinterpolar(nr, rofi, v, gridn, pol_r, ipol_v);
  
    /*  for (x_u = 0 ; x_u < gridn ; x_u++ ){
	printf("\n");
	for (y_u = 0 ; y_u < gridn ; y_u++ ){
	printf("\n");
	for (z_u = 0 ; z_u < gridn ; z_u++ ){
	printf("%lf", ipol_v[x_u][y_u][z_u]);
	}      
	}      
	}*/
    // rseq data
    emin = V0 ;
    emax = 100. ;
    iz = 0;
  
    for( l1 = 0 ; l1 < lmax ; l1++ ){
      l = l1 ;
      for( n1 = 0 ; n1 < nodn ; n1++ ){
	if(n1 % nodn < cln){
	  nod = n1 % nodn ; val = 0. ;  slo = -1. ; e = 0. ;
	}else if(n1 % nodn >= cln){
	  nod = n1 % nodn - cln ; val = 1. ;  slo = 0. ; e = 0. ;
	}else{
	  printf("fault\n");
	  exit(0);
	}
	rseq(a, b, &e, emin, emax, g, gfac, l, nod, nr, &nre, rofi, &slo, v, &val, iz ) ;
	allslo[l1][n1] = slo;
	allval[l1][n1] = val;
	alle[l1][n1] = e;
	for( i = 0 ; i < nr ; i++ ){
	  allg[ l1 ][ n1 ][ i ] = g[ i ] * sqrt( gfac[ i ] );
	}
	sp_sym_3dinterpolar(nr, rofi, allg[l1][n1], gridn, pol_r, ipol_g);
	for( x_u = 0 ; x_u < gridn ; x_u++ ){
	  for( y_u = 0 ; y_u < gridn ; y_u++ ){
	    for( z_u = 0 ; z_u < gridn ; z_u++ ){
	      ipol_allg[l1][n1][x_u][y_u][z_u] = ipol_g[x_u][y_u][z_u];
	    }
	  }
	}
      }
      for( m1 = 0 ; m1 < 2 * l1 + 1 ; m1++ ){
	// printf("\n");
	for( x_u = 0 ; x_u < gridn ; x_u++ ){
	  for( y_u = 0 ; y_u < gridn ; y_u++ ){
	    for( z_u = 0 ; z_u < gridn ; z_u++ ){
	      if(pol_r[x_u][y_u][z_u] != 0 && pow(gridx[x_u],2.) + pow(gridx[y_u],2.) != 0){
		if(gridx[y_u] < 0){
		  sgnf = -1.;
		}else if(gridx[y_u] == 0){
		  sgnf = 0.;
		}else if(gridx[y_u] > 0){
		  sgnf = 1.;
		}
		pol_phi = sgnf * acos(gridx[x_u]/sqrt(pow(gridx[x_u],2.) + pow(gridx[y_u],2.)));
		ally[l1][m1][x_u][y_u][z_u] = cylm(l1, m1-l1, gridx[z_u] / pol_r[x_u][y_u][z_u], pol_phi);
		// printf("%f", creal(ally[l1][m1][x_u][y_u][z_u]));
	      }
	      if(pol_r[x_u][y_u][z_u] <= wsr ){
		upot[x_u][y_u][z_u] = cellpot[x_u][y_u][z_u] - ipol_v[x_u][y_u][z_u];
	      
	      }
	    }
	  }
	}
      }
    }
  
    /*  for (x_u = 0 ; x_u < gridn ; x_u++ ){
	printf("\n");
	for (y_u = 0 ; y_u < gridn ; y_u++ ){
	printf("\n");
	for (z_u = 0 ; z_u < gridn ; z_u++ ){
	printf("%lf", upot[x_u][y_u][z_u]);
	}      
	}      
	}*/
  
    // matrix calculation
  
    for ( l1 = 0 ; l1 < lmax ; l1++){
      for ( n1 = 0 ; n1 < nodn ; n1++){
	for( l2 = 0 ; l2 < lmax ; l2++){
	  for ( n2 = 0 ; n2 < nodn ; n2++){
	    if ( l1 == l2 ){
	      qmat[l1][l2][n1][n2] = allval[ l1 ][ n1 ] * allval[ l2 ][ n2 ];
	      lmat[l1][l2][n1][n2] = allval[ l1 ][ n1 ] * allslo[ l2 ][ n2 ];
	      for(i = 0; i < nr; i++){
		hsg[ i ] = allg[ l1 ][ n1 ][ i ] * allg[ l2 ][ n2 ][ i ];
	      }
	      hsmat[ l1 ][ l2 ][ n1 ][ n2 ] = defintegr( hsg, rofi, nr ) * alle[ l2 ][ n2 ];
	    }
	    for(m1 = 0; m1 < l1 * 2 + 1; m1++){
	      for(m2 = 0; m2 < l2 * 2 + 1; m2++){
		for( x_u = 0 ; x_u < gridn ; x_u++ ){
		  for( y_u = 0 ; y_u < gridn ; y_u++ ){
		    for( z_u = 0 ; z_u < gridn ; z_u++ ){
		      if (pol_r[x_u][y_u][z_u] != 0){
			umat[l1][l2][n1][n2][m1][m2] += ipol_allg[l1][n1][x_u][y_u][z_u] * ipol_allg[l2][n2][x_u][y_u][z_u] * conj(ally[l1][m1][x_u][y_u][z_u]) * ally[l2][m2][x_u][y_u][z_u] * upot[x_u][y_u][z_u] * delta_v / pow(pol_r[x_u][y_u][z_u], 2.);
			// printf("%lf %lf %lf %lf %lf %lf %lf\n", ipol_allg[l1][n1][x_u][y_u][z_u], ipol_allg[l2][n2][x_u][y_u][z_u], creal(conj(ally[l1][m1][x_u][y_u][z_u])), creal(ally[l2][m2][x_u][y_u][z_u]), upot[x_u][y_u][z_u], delta_v, pow(pol_r[x_u][y_u][z_u], 2.))  ;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  
  
  
    /*  printf("\nL=0,close1\tL=0,close2\tL=0,open1...\n");
	printf("\nQ matrix\n");
	for ( l1 = 0 ; l1 < lmax ; l1++){
	for ( n1 = 0 ; n1 < nodn ; n1++){
	for( l2 = 0 ; l2 < lmax ; l2++){
	for ( n2 = 0 ; n2 < nodn ; n2++){
	if (l2 == lmax - 1 && n2 == nodn - 1){
	printf("%9.6lf\n", qmat[l1][l2][n1][n2]);
	}else{
	printf("%9.6lf\t", qmat[l1][l2][n1][n2]);
	}
	}
	}
	}
	}
  
	printf("\nL matrix\n");
  
	for ( l1 = 0 ; l1 < lmax ; l1++){
	for ( n1 = 0 ; n1 < nodn ; n1++){
	for( l2 = 0 ; l2 < lmax ; l2++){
	for ( n2 = 0 ; n2 < nodn ; n2++){
	if (l2 == lmax - 1 && n2 == nodn - 1){
	printf("%9.6lf\n", lmat[l1][l2][n1][n2]);
	}else{
	printf("%9.6lf\t", lmat[l1][l2][n1][n2]);
	}
	}
	}
	}
	}
  
	printf("\nHs matrix\n");
	for ( l1 = 0 ; l1 < lmax ; l1++){
	for ( n1 = 0 ; n1 < nodn ; n1++){
	for( l2 = 0 ; l2 < lmax ; l2++){
	for ( n2 = 0 ; n2 < nodn ; n2++){
	if (l2 == lmax - 1 && n2 == nodn - 1){
	printf("%9.6lf\n", hsmat[l1][l2][n1][n2]);
	}else{
	printf("%9.6lf\t", hsmat[l1][l2][n1][n2]);
	}
	}
	}
	}
	}
  
	printf("\nu matrix\n");
	for ( l1 = 0 ; l1 < lmax ; l1++){
	for ( n1 = 0 ; n1 < nodn ; n1++){
	for ( m1 = 0 ; m1 < l1 * 2 + 1; m1++){
	for( l2 = 0 ; l2 < lmax ; l2++){
	for ( n2 = 0 ; n2 < nodn ; n2++){
	for ( m2 = 0 ; m2 < l2 * 2 + 1; m2++){
	printf("l1 = %d\tn1 = %d\tm1 = %d\tl2 = %d\tn2 = %d\tm2 = %d\n", l1, n1, m1, l2, n2, m2);
	printf("%9.6e\n", creal(umat[l1][l2][n1][n2][m1][m2]));
	}
	}
	}
	}
	}
	}
    */

    //  printf("\nl1 = 3\tn1 = 2\tm1 = 3\tl2 = 3\tn2 = 2\tm2 = 3\n");
    //  printf("%9.6e\n", creal(umat[3][3][2][2][6][6]));

    fprintf(sf,"%d\t",gridn);
    for ( l1 = 0 ; l1 < lmax ; l1++){
      for ( n1 = 0 ; n1 < nodn ; n1++){
	for ( m1 = 0 ; m1 < l1 * 2 + 1; m1++){
	  for( l2 = 0 ; l2 < lmax ; l2++){
	    for ( n2 = 0 ; n2 < nodn ; n2++){
	      for ( m2 = 0 ; m2 < l2 * 2 + 1; m2++){
		fprintf(sf,"%9.6e\t", creal(umat[l1][l2][n1][n2][m1][m2]));
	      }
	    }
	  }
	}
      }
    }
    
    free(gridx);
    free(rofi); 
    free(v);    
    free(gfac);
    free(g);   
  
    for ( x_u = 0; x_u < gridn ; x_u++ ){
      for ( y_u = 0; y_u < gridn ; y_u++ ){
	free(cellpot[x_u][y_u]);
      }
      free(cellpot[x_u]);
    }
    free(cellpot);

  
    for ( l1 = 0; l1 < lmax ; l1++ ){
      for ( l2 = 0; l2 < lmax ; l2++ ){
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  free(qmat[l1][l2][n1]);
	}
	free(qmat[l1][l2]);
      }
      free(qmat[l1]);
    }
    free(qmat);

    for ( l1 = 0; l1 < lmax ; l1++ ){
      for ( l2 = 0; l2 < lmax ; l2++ ){
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  free(lmat[l1][l2][n1]);
	}
	free(lmat[l1][l2]);
      }
      free(lmat[l1]);
    }
    free(lmat);
  
    for ( l1 = 0; l1 < lmax ; l1++ ){
      for ( l2 = 0; l2 < lmax ; l2++ ){
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  free(hsmat[l1][l2][n1]);
	}
	free(hsmat[l1][l2]);
      }
      free(hsmat[l1]);
    }
    free(hsmat);
  
    for ( l1 = 0; l1 < lmax ; l1++ ){
      for ( l2 = 0; l2 < lmax ; l2++ ){
	for ( n1 = 0; n1 < nodn ; n1++ ){
	  for ( n2 = 0; n2 < nodn ; n2++ ){
	    for ( m1 = 0; m1 < l1 * 2 + 1 ; m1++ ){
	      free(umat[l1][l2][n1][n2][m1]);
	    }
	    free(umat[l1][l2][n1][n2]);
	  }
	  free(umat[l1][l2][n1]);
	}
	free(umat[l1][l2]);
      }
      free(umat[l1]);
    }
    free(umat);
  
    for ( l1 = 0; l1 < lmax; l1++ ){
      free(allslo[l1]);
    }
    free(allslo);
  
    for ( l1 = 0; l1 < lmax; l1++ ){
      free(allval[l1]);
    }
    free(allval);
  
    for ( l1 = 0; l1 < lmax; l1++ ){
      free(alle[l1]);
    }
    free(alle);
  
    for ( l1 = 0; l1 < lmax ; l1++ ){
      for ( n1 = 0; n1 < nodn ; n1++ ){
	free(allg[l1][n1]);
      }
      free(allg[l1]);
    }
    free(allg);
  
    free(hsg) ;
  
    for ( l1 = 0; l1 < lmax; l1++ ){
      for ( m1 = 0; m1 < l1 * 2 + 1; m1++ ){
	for ( x_u = 0; x_u < gridn; x_u++ ){
	  for ( y_u = 0; y_u < gridn; y_u++ ){
	    free(ally[l1][m1][x_u][y_u]);
	  }
	  free(ally[l1][m1][x_u]);
	}
	free(ally[l1][m1]);
      }
      free(ally[l1]);
    }
    free(ally);
  
    for ( l1 = 0; l1 < lmax; l1++ ){
      for ( n1 = 0; n1 < nodn; n1++ ){
	for ( x_u = 0; x_u < gridn; x_u++ ){
	  for ( y_u = 0; y_u < gridn; y_u++ ){
	    free(ipol_allg[l1][n1][x_u][y_u]);
	  }
	  free(ipol_allg[l1][n1][x_u]);
	}
	free(ipol_allg[l1][n1]);
      }
      free(ipol_allg[l1]);
    }
    free(ipol_allg);
  
    for ( x_u = 0; x_u < gridn ; x_u++ ){
      for ( y_u = 0; y_u < gridn ; y_u++ ){
	free(pol_r[x_u][y_u]);
	free(ipol_v[x_u][y_u]);
	free(upot[x_u][y_u]);
	free(ipol_g[x_u][y_u]);
      }
      free(pol_r[x_u]);
      free(ipol_v[x_u]);
      free(upot[x_u]);
      free(ipol_g[x_u]);
    }
    free(pol_r);
    free(ipol_v);
    free(upot);
    free(ipol_g);

    
    fprintf(sf,"\n");
    time_end = clock();
    printf("total time = %.6f sec\n", (double)(time_end - time_start) / CLOCKS_PER_SEC );
  }

  fclose(sf);
  
  /*  time_end = clock();
      printf("total time = %.6f sec\n", (double)(time_end - time_start) / CLOCKS_PER_SEC );*/
}
