#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define max(a,b) ( (a) > (b) ? (a) : (b) )

char * addnumber( char * str1, int n ) {
  int i , m = strlen(str1) ;
  char * str2 = malloc( m + 4 ) ;
  strcpy(str2, str1 ) ;
  strcat(str2, "_00" ) ;
  for ( i = 0 ; i < 2 ; i++ ) { str2[m+2-i] = '0' + n%10 ; n /= 10 ; }
  return str2 ;
}

void lmtomeshparam(int z, double wsr, double *a, double *b, int *nr) {
  *b = 1. / ( z + z + 1. ) ;
  *a = 0.03 ;
  *nr = 2. * ( 0.5 + log( 1. + wsr / *b ) / *a ) ;
  *nr = max( 25, ( (*nr - 1) / 2 ) * 2 + 1 ) ;
  *b = wsr / ( exp( *a * ( *nr - 1 ) ) - 1. ) ;
}

void interpolar(int nin, double *xin, double *yin,
                int nout, double *xout, double *yout) ;


int main( int argc, char *argv[] )
{
  char   symbl[6] ;
  int    neq, nz, kmax, kplace, chg[10], k, i, ch, lastline ;
  int    nr, ineqat, nineqat, nat ;
  double x, y, z, mtr, *rhs, *rvhs, vhs, vcon, hh ;
  double a, b, *rlmto, *rvlmto ;
  FILE *fp, *fpop, *fpstr ;
 
  if ( argc < 2 ) { printf("usage: >readhspot filename\n"); exit(1) ; }

  
  fp = fopen(argv[1],"r") ;


  ineqat = 0 ;
  fscanf(fp,"%d%d%lf", &nineqat, &nat, &vcon) ; 
/*    printf("%d %d %lf\n", nineqat, nat, vcon) ; */

  fpstr = fopen("instr.dat","w") ;

  fprintf(fpstr," %4d %4d    2\n", nat, nat ) ;  

  fscanf(fp,"%s", symbl) ;   /*    printf("%s ",symbl) ;   */
  while ( strcmp(symbl,"XXXXX") ) { 
    fscanf(fp,"%d %d %d %d %lf %lf %lf %lf",
	   &neq,&nz,&kmax,&kplace,&x,&y,&z,&mtr ) ;  
/*      printf("%d %d %d %d %lf %lf %lf %lf \n", */
/*  	     neq,nz,kmax,kplace,x,y,z,mtr ) ; */
    if ( neq == 0 ) {
      rhs    = ( double * ) malloc( ( kmax + 1 ) * sizeof( double ) ) ;    
      rvhs   = ( double * ) malloc( ( kmax + 1 ) * sizeof( double ) ) ;

      for ( i = 0 ; i < 10 ; i++ ) {
	fscanf(fp,"%d", &chg[i] ) ;
      }
/*       for ( i = 0 ; i < 10 ; i++ ) printf("%d ", chg[i] ) ; printf("\n");*/

      rhs[0] = 0. ;
      for ( i = 1 ; i <= 5 ; i++ ) {
	fscanf(fp,"%lf", &rhs[i] ) ;
      }
/*        for ( i = 0 ; i <= 5 ; i++ ) printf("%lf ", rhs[i] ); printf("\n");*/
      ch = 0 ;
      hh = rhs[2] - rhs[1] ;
      for ( k = 3 ; k <= kmax ; k++ ) {
	rhs[k] = rhs[k-1] + hh ;
	if ( k == chg[ch] ) {
	  hh += hh ;
	  ch++ ;
	}
      }
      
      rvhs[0] = 0. ; 
      for ( k = 1 ; k <= kmax ; k++ ) {
	fscanf(fp,"%lf", &vhs ) ;
	rvhs[k] = rhs[k] * vhs + 2. * nz ;
        printf("%14.7e %14.7e %14.7e\n", rhs[k], rvhs[k], vhs ) ;
      }

      
      lmtomeshparam(nz, mtr, &a, &b, &nr) ;
      rvlmto = ( double * ) malloc( nr * sizeof( double ) ) ;
      rlmto  = ( double * ) malloc( nr * sizeof( double ) ) ;


      for (i=0; i<nr; i++) {
	rlmto[i] = b*(exp(a*i)-1.) ; 
      }

      rvlmto[0] = 0. ;
      interpolar(kmax, rhs, rvhs, nr-1, rlmto+1, rvlmto+1) ;
      
      ineqat++;
      
      fpop = fopen( addnumber(symbl,ineqat), "w") ;
      fprintf(fpop,"LMTO-type pot-file generated from H.S. pot-file\n") ;
      fprintf(fpop,"Z= %d\nPOT:\n", nz ) ;
      fprintf(fpop," %4d %4d %11.5lf %11.5lf\n", nr, 1, a, mtr ) ;
      fprintf(fpop," %15.9e",0.) ;
      for ( i = 1 ; i < nr ; i++ ) { 
        if (i%5==0) fprintf(fpop,"\n");
        fprintf(fpop," %15.9e", rvlmto[i] / rlmto[i] ) ;
      }
      fprintf(fpop,"\n");
      fclose(fpop);

      free( rvlmto ) ; free( rlmto ) ; free( rvhs ) ; free( rhs ) ;

      fprintf(fpstr,"%6s",addnumber(symbl,ineqat)) ;
    } else {
      fprintf(fpstr,"%6s",addnumber(symbl,neq)) ;
    }
    fprintf(fpstr," %2d %3d %11.6lf %11.6lf %11.6lf\n", 
      ( (ineqat==1 || nz>38) ? 3 : ( nz > 12 ? 2 : 1 ) ), nz, x, y, z ) ; 
    fscanf(fp,"%s", symbl) ; /* printf("%s ",symbl) ;  */
  }
  fclose(fpstr) ;  fclose(fp) ;
  if ( nineqat != ineqat ) {
    printf("ERROR in readrinopot.c: nineqat != ineqat\n") ;
    exit(1) ;
  }

  return 0 ;
}





double cubicinterpol(double x0, double y0, double x1, double y1,
		     double x2, double y2, double x3, double y3, double x ){
/* given four points (x_i,y_i) i=0,3, the exact polynomial (cubic) inter-
   polation is calculated and the y-value corresponding to x is returned */
  double yd0, yd1, yd2, yd3, a0, a1, a2, a3;
  yd0 = y0/(x0-x1)/(x0-x2)/(x0-x3) ;
  yd1 = y1/(x1-x0)/(x1-x2)/(x1-x3) ;
  yd2 = y2/(x2-x0)/(x2-x1)/(x2-x3) ;
  yd3 = y3/(x3-x0)/(x3-x1)/(x3-x2) ;
  a0 = - x1*x2*x3*yd0 - x0*x2*x3*yd1 - x0*x1*x3*yd2 - x0*x1*x2*yd3 ;
  a1 =  (x1*x2 + x1*x3 + x2*x3) * yd0
      + (x0*x2 + x0*x3 + x2*x3) * yd1 
      + (x0*x1 + x0*x3 + x1*x3) * yd2 
      + (x0*x1 + x0*x2 + x1*x2) * yd3 ;
  a2 = - (x1+x2+x3)*yd0 - (x0+x2+x3)*yd1 - (x0+x1+x3)*yd2 - (x0+x1+x2)*yd3 ;
  a3 = yd0 + yd1 + yd2 + yd3 ;
  return a0 + a1*x + a2*x*x + a3*x*x*x ;
}

int findindex(double x0, int n, double *x) {
  /* returns index i0 such that x[i0] <= x0 < x[i0+1] 
     if x0 < x[i] forall i  => return  -1
     only works with pos. ordered *x, i.e. i < j => xin[i] < xin[j].
  */
  int low, high, mid ;

  if ( x0 < x[0] ) return -1 ;
  if ( x[n-1] <= x0 ) return n-1 ;

  low = 0 ; high = n - 1 ;  
  while ( low < high - 1 ) {
    mid = ( low + high ) / 2 ;
    if ( x[mid] <= x0 ) low = mid ;
    else high = mid ;
  }
  return low ;
}

void interpolar(int nin, double *xin, double *yin, 
		int nout, double *xout, double *yout){
/*
   calculates interpolated y-values yout[i] for x-values xout[i] (i=0,nin-1)
   from dataset (xin[j], yin[j]), j=0,nout-1.
   input:   nin, *xin, *yin, nout, *xout
   output:  *yout
   Uses exact cubic interpolation given in function "cubicinterpol" with 
   the 4 points j=i0-1,...,i0+2, where i0 is def. by xin[i0] <= x < xin[i0+1].
   IMPORTANT CONDITIONS:
   nin >= 4 AND *xin must be positively ordered, i.e. i < j => xin[i] < xin[j]
*/
  int i, i0 ;

/* if ( nin < 4 ) exit(1) */           /* check if nin < 4 */
/* check if  xin is positively ordered, i < j => xin[i] < xin[j] */
/* for ( i = 1 ; i < nin ; i++ ) if ( xin[i-1] > xin[i] ) exit(1) ;   */

  for (i = 0 ;  i < nout ; i++ ){
    i0 = findindex(xout[i], nin, xin) ;  /* xin[i0] <= x < xin[i0+1] */
    if ( i0 < 1 ) i0 = 1 ;
    else if ( i0 > nin - 3 ) i0 = nin - 3 ;
    yout[i] = cubicinterpol(xin[i0-1], yin[i0-1], xin[i0], yin[i0], 
	       xin[i0+1], yin[i0+1], xin[i0+2], yin[i0+2], xout[i] ) ;
  }
}

#undef max
