#include <stdio.h>
#include <string.h>
#define MAXWORD 132

void readinput(char* absorbergsfile, char* absorberchfile,
               int* nshell, char* corewffile, double* eenvmin, double* eenvmax,
               double* deenv, int* lmax, int* lsym, double* ksicor, 
	       double *ksival, double *tendqval, double *hz, double* vmtzero, 
	       double* imvmtzero, double* efermi, double* omegamin, 
               double* omegamax, double* domega, 
               double* vmixpar, double* rkscale,
               double* e0min, double* e0max, int* nrwf1,
	       int* lvalsh, int* nvalel, int* withrefl, int myid)
/* omega = photon E = total final E. Photo-el E = omega + ecore, since Eg:=0 */
{
  char token[MAXWORD] ;

  FILE *fp ; 

  *lsym = -1 ;
  *ksival = 0. ;
  *ksicor = 0. ;
  *tendqval = 0. ; 
  *hz = 0. ;
  *imvmtzero = 0. ;
  *vmixpar = 0. ;
  *rkscale = 1. ;
  *e0min = -1. ;
  *e0max = 1. ;
  *nrwf1 = 1 ;
  *lvalsh = 2 ;
  *nvalel = 0 ;
  *withrefl = 1 ;
  *efermi = -1.e10 ;
  *nshell = 0 ;
  fp = fopen("instr.dat","r") ;
  fscanf(fp, "%d%*d%*d%s%d", nshell, absorbergsfile,lmax ) ;
  fclose( fp ) ;
  if ( *nshell > 0 ) (*nshell)-- ;
  strcpy(absorberchfile,absorbergsfile) ; 

  fp = fopen("mcms.in","r") ;
  while( fscanf(fp,"%s", token) != EOF ) {
    if ( strcmp(token,"AbsCHPotFile")==0 )     fscanf(fp,"%s",  absorberchfile ) ;
    else if ( strcmp(token,"CoreWFFile")==0 )       fscanf(fp,"%s",  corewffile ) ;
    else if ( strcmp(token,"PhotoelEMin")==0 )      fscanf(fp,"%lf", eenvmin ) ;  
    else if ( strcmp(token,"PhotoelEMax")==0 )      fscanf(fp,"%lf", eenvmax ) ;  
    else if ( strcmp(token,"DeltaPhotoelE")==0 )    fscanf(fp,"%lf", deenv ) ;    
    else if ( strcmp(token,"Lsym")==0 )             fscanf(fp,"%d",  lsym ) ;     
    else if ( strcmp(token,"KsiCor")==0 )           fscanf(fp,"%lf", ksicor ) ;
    else if ( strcmp(token,"KsiVal")==0 )           fscanf(fp,"%lf", ksival ) ;
    else if ( strcmp(token,"TendqVal")==0 )         fscanf(fp,"%lf", tendqval ) ;
    else if ( strcmp(token,"Hz")==0 )               fscanf(fp,"%lf", hz ) ;
    else if ( strcmp(token,"ReVmtz")==0 )           fscanf(fp,"%lf", vmtzero ) ;  
    else if ( strcmp(token,"ImVmtz")==0 )           fscanf(fp,"%lf", imvmtzero ) ;
    else if ( strcmp(token,"EFermi")==0 )           fscanf(fp,"%lf", efermi ) ;
    else if ( strcmp(token,"OmegaMin")==0 )         fscanf(fp,"%lf", omegamin ) ; 
    else if ( strcmp(token,"OmegaMax")==0 )         fscanf(fp,"%lf", omegamax ) ;
    else if ( strcmp(token,"DeltaOmega")==0 )       fscanf(fp,"%lf", domega ) ;
    else if ( strcmp(token,"UnscreenedWeight")==0 ) fscanf(fp,"%lf", vmixpar ) ;
    else if ( strcmp(token,"SlaterRkScale")==0 )    fscanf(fp,"%lf", rkscale ) ;  
    else if ( strcmp(token,"BasisEmin")==0 )        fscanf(fp,"%lf", e0min ) ;    
    else if ( strcmp(token,"BasisEmax")==0 )        fscanf(fp,"%lf", e0max ) ;    
    else if ( strcmp(token,"NOpenChannels")==0 )    fscanf(fp,"%d",  nrwf1 ) ;
    else if ( strcmp(token,"LValSh")==0 )           fscanf(fp,"%d",  lvalsh ) ;
    else if ( strcmp(token,"NValEl")==0 )           fscanf(fp,"%d",  nvalel ) ;
    else if ( strcmp(token,"WithRefl")==0 )         fscanf(fp,"%d",  withrefl ) ;
  }

  fclose(fp) ;

/*print it all out */ 
  if (myid == 0){
  printf("%s\t = %s\n", "AbsGSPotFile",    absorbergsfile ) ; 
  printf("%s\t = %s\n", "AbsCHPotFile",    absorberchfile ) ;
  printf("%s\t = %s\n", "CoreWFFile",      corewffile ) ; 
  printf("%s\t\t = %d\n", "Lmax",          *lmax ) ;     
  printf("%s\t\t = %d\n", "Lsym",          *lsym ) ; 
  printf("%s\t\t = %lf\n","KsiCor",        *ksicor ) ;     
  printf("%s\t\t = %lf\n","KsiVal",        *ksival ) ;     
  printf("%s\t   = %lf\n","TendqVal",      *tendqval ) ;     
  printf("%s\t\t = %lf\n","Hz",            *hz ) ;     
  printf("%s\t\t = %lf\n","ReVmtz",        *vmtzero ) ;  
  printf("%s\t\t = %lf\n","ImVmtz",        *imvmtzero ) ;
  printf("%s\t\t = %lf\n","EFermi",        *efermi ) ;
  printf("%s\t = %lf\n","OmegaMin",        *omegamin ) ; 
  printf("%s\t = %lf\n","OmegaMax",        *omegamax ) ; 
  printf("%s\t = %lf\n","DeltaOmega",      *domega ) ;   
  printf("%s  = %lf\n","UnscreenedWeight", *vmixpar ) ;  
  printf("%s\t = %lf\n","SlaterRkScale",   *rkscale ) ;  
  printf("%s\t = %lf\n","BasisEmin",       *e0min ) ;    
  printf("%s\t = %lf\n","BasisEmax",       *e0max ) ;    
  printf("%s\t = %d\n","NOpenChannels",    *nrwf1 ) ; 
  printf("%s\t\t = %d\n","LValSh",         *lvalsh ) ;
  printf("%s\t\t = %d\n","NValEl",         *nvalel ) ;
  printf("%s\t = %d\n","WithRefl",         *withrefl ) ;
  printf("%s\t\t = %d\n", "Nshell",          *nshell ) ;   
  printf("%s\t = %lf\n","PhotoelEMin",     *eenvmin ) ;  
  printf("%s\t = %lf\n","PhotoelEMax",     *eenvmax ) ;  
  printf("%s\t = %lf\n","DeltaPhotoelE",   *deenv ) ;   
  } 
}
