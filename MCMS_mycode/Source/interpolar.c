double cubicinterpol(double x0, double y0, double x1, double y1,
		     double x2, double y2, double x3, double y3, double x ) ;

int findindex(double x0, int n, double *x) ;

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
