#include <stdio.h>
#include <stdlib.h>

double derivquadr( double x0, double y0, double x1, double y1, 
                   double x2, double y2, double x ) 
{
/* given three points (x_i,y_i) i=0,1,2 the quadratic interpolation 
  ( y = a0 + a1 * x + a2 * x * x ) is calculated and the derivative at 
   x is returned 
*/ 
  double yd0, yd1, yd2, a1, a2 ;

  yd0 = y0/(x0-x1)/(x0-x2) ;
  yd1 = y1/(x1-x0)/(x1-x2) ;
  yd2 = y2/(x2-x0)/(x2-x1) ;

/*    a0 is given as follows but it is not needed for derivative  */
/*    a0 =   x1*x2 * yd0 + x0*x2 * yd1 + x0*x1 * yd2 ;            */

  a1 = - (x1+x2)*yd0 - (x0+x2)*yd1 - (x0+x1)*yd2 ;
  a2 = yd0 + yd1 + yd2 ;

  return  2. * a2 * x  + a1 ;
}

void derivative( int n, double *x, double *y, double *dy )
{
/* input:  (x_i,y_i). 
   output: derivative ( dy_i = dy/dx(x_i) ) using quadratic interpotation 
*/
  int i ;
  if ( n < 3 ) {
    printf("ERROR in void derivative: #pts = n = %d < 3\n", n ) ; exit(1) ;
  } else {
    dy[0] = derivquadr( x[0], y[0], x[1], y[1], x[2], y[2], x[0] ) ;
    for ( i = 1 ; i < n-1 ; i++ ) 
      dy[i] = derivquadr( x[i-1], y[i-1], x[i], y[i], x[i+1], y[i+1], x[i] );
    dy[n-1] = derivquadr( x[n-3],y[n-3],x[n-2],y[n-2],x[n-1],y[n-1], x[n-1]);
  }
}
