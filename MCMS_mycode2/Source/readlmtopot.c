#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MAXLINE 100

struct Lmtopot { 
  int z; double wsr; double a; double b; int nr; int nspin; double *v; 
};


struct Lmtopot readlmtopot(char *filename) {
/* 
   reads potential from file potfilename into struct LmtoAtPot, 
   which is returned. Allocates sufficient memory for v. 
*/
    int c, i, z, nr, nspin ;
    char line[MAXLINE] ;
    double wsr, a, vi ;
    struct Lmtopot atom ;
    FILE * fppot ;

    fppot = fopen(filename,"r");
    while((c=fgetc(fppot))!=EOF&&c!='\n');  
    while((c=fgetc(fppot))!='='); 
    fscanf(fppot,"%d",&z);
    printf("in readlmtopot: just read: z= %d\n", z);
    while(fgets(line,MAXLINE,fppot) && strncmp(line,"POT:",4));
    fscanf(fppot,"%d%d%lf%lf",&nr,&nspin,&a,&wsr);
    atom.z     = z ;
    atom.wsr   = wsr ;
    atom.a     = a ;
    atom.b     = wsr/(exp(a*(nr-1))-1.) ;
    atom.nr    = nr ;
    atom.nspin = nspin ;
    atom.v     = ( double * ) malloc( nr * nspin * sizeof( double ) ) ;
    for ( i = 0 ;  i < nr * nspin ; i++ ) {
      fscanf(fppot,"%lf", &vi );
      atom.v[i] = vi ;
    }
    fclose(fppot);

    return atom ;
}

    /*    
	  printf("# z=%d, wsr=%lf, a=%lf, b=%lf, nr=%d, nspin=%d\n", 
	  atom.z, atom.wsr, atom.a, atom.b, atom.nr, atom.nspin ) ;
    */
