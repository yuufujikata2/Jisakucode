#include<stdio.h>
#include<math.h>


double ***calloc3double(int n1, int n2, int n3){
  int i1, i2;
  double ***p;
  p = ( double *** ) calloc( n1, sizeof( double ** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( double ** ) calloc( n2, sizeof( double * ) );
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( double * ) calloc( n3, sizeof( double ) );
    }
  }
  return p ;
}


int main()
{
  double*** v ;
  int LX = 20, LY = 20, LZ = 20 ;
  int PX = 10, PY = 10, PZ = 10 ;
  int i, j, k, l ;

  v = calloc3double(2*LX+1,2*LY+1,2*LZ+1) ;
  for (i = ; i < 31; i++){
    for (j = 10; j < 31; j++){
      for (k = 10; k < 31; k++){
	v[i][j][k] = 1.
      }
    }
  }







  return 1;
}
