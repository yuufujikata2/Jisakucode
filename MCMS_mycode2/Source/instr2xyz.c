#include <stdio.h>
#define BOHR 0.529177
int main()
{
  char symbl[31][3] = { "  ", " H", "He", 
   "Li", "Be", " B", " C", " N", " O", " F", "Ne",
   "Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar",
   " K", "Ca", "Sc", "Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
   "Zn" } ;

  int nat, i, iz ;
  double x, y, z ;
  FILE *fpstr, *fpxyz ;

  fpstr = fopen("instr.dat","r") ;
  fpxyz = fopen("str.xyz","w") ;
  fscanf(fpstr,"%d%*d%*d", &nat ) ; 
  fprintf(fpxyz,"%d\n\n", nat ) ;
  for ( i = 0 ; i < nat ; i++ ) {   
   fscanf(fpstr,"%*s%*d%d%lf%lf%lf", &iz, &x, &y, &z) ;
   if ( iz <= 30 ) fprintf(fpxyz," %2s ", symbl[iz] ) ;
   else fprintf(fpxyz,"%3d ", iz ) ;
   fprintf(fpxyz," %10.5lf %10.5lf %10.5lf\n",BOHR*x,BOHR*y,BOHR*z) ;
  }
}
