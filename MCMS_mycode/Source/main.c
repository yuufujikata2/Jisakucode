#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "w3j.h"
#include "ndimarraycalloc.h"
#include "globals.h"
#include "headman.h"
#include <mpi.h>

#define FILENAMELENGTH 270
#define NSPIN 1
#define RYDINEV 13.6058
#define QMETRICZERO 1.e-5
#define BETAZERO 1.e-10
#define TEMPERATURE 0.01
#define LMAX_cont 3

int main(int argc , char **argv){
    char c, corewffile[FILENAMELENGTH] ; 
    char absorbergsfile[FILENAMELENGTH], absorberchfile[FILENAMELENGTH] ;
    int nshell, iz, nr, nspin, ncorsh, lcorsh, lvalsh, nvalel, withrefl ;
    int i,j,k, k1, k2, ijre, ir ;
    int ich0, jch0, ich, jch, ndeg, nrwf, nrwf1, nrwfall, nod0 ;
    int neenv, lsym, lmax, ntas, ntasc, nnzele, lm1, lm2 ;
    int  m1, m2, lmmax, ims1, ims2, ispin, jspin ; 
    double eenvmin, eenvmax, deenv, *eenv ;
    double a, b, wsr, ksicor, ksival, tendqval, ecore, vmtzero, imvmtzero ;
    double e0min, e0max, efermi, deltaEoverT  ;
    double vmixpar, rkscale, omegamin, omegamax, domega, omega ; 
    double deltaij, kgai, kgai0, xij ;
    double *vofi, *vcoree, *cofi, *rofi, ****rho, ****chi, **wnwj ;
    double *e0p, *pr0, *dpr0, **ppofi, **e0po, *pr0o, *dpr0o, **e0pomod ;
    double *kga, *rkbes, *drkbes, *rkneu, *drkneu ; 
    double *alpha, *beta, *redipmejj, *redipme, *mat1 ;
    double *tau, *lambda, *qmetric ;
    double r, rerho, imrho, rechi, imchi, imtauij ;
    double sum, sumi, eout, *sigmaqqgst[3][3] ; 
    double *ppofi0, **pqofi, **smat, **tmat ;
    int nopen1, *lsh, *sorb1sh, nsorbs, nconfs, **occ, nstates, nhamele ;
    int nlevels, ngstdeg, ngstbasis, inishell, ndiplistele, iq, q, iqp ;
    int nfstbasis, fst0deg, igstdeg, nqmetricele, ist, info ;
    int nchabasis, ncha0basis, nchannels, valshdeg, contsh1, contsorb1 ; 
    int counter, sign, nwoplistele, nopele, *jk, ims, i1, valsh ;    
    int nfstresbas, ifstresbas, jfstresbas, nfstboundbas, nfstcontbas ;
    int ncontshs, icontsh, ish, jchabas ; 
    int nvalshholes, jms, jorb, ihole, jhole, n1, n2 ;
    double *ksish, *tendqsh, *ham, *eneigval, gstenergy, **gstvec ;
    double *cha0energy, **cha0vec, *chaenergy, ***chavec, ener ;
    double hz, *fstvec, fst0energy, realk, x, jl, nl, jlp, nlp ;
    double *alphalong, *betalong, *revec, *normfac, *wdagger ;
    double **dipmatfbas[3], **dipmatchabas[3], **dipmatele[3] ;
    double **dipmatfresbas[3] ;
    double dum, wij, wik, wkj, mat2kj, *mat1linei ;
    double ***gstholevecs, ***fstbas2res ;
    double imtautrace, *dumvec, **dummat, sumre, sumim, **dummattimesu ;
    double *valshdensmat, *valshdenseval ;
    struct Lmtopot gspot, chpot, pot0 ;
    struct Fock *gstbasis, *cha0basis, *chabasis, *fstbasis, statek ;
    struct Spamaline *hamsparse, *pdipsmline0[3], *pwopsmline0, *pwopsml ;
    struct Spamaline **pwopressml, **hamressparse, **qressparse ;
    struct Spamaline *qmetricsparse, *psmline, *pvalshdensmatsmline0 ;
    struct O1plistitem *po1plistitem ;
    struct O1p wop, dipop[3], valshdensmatop ;
    FILE *fp, *fpm, *fpq ;

    double energy0, bes, dbes, neu, dneu, singlepr0, singledpr0 ;
    double rWjR, rWnR, wnwj0, radipme0, zfac, singlesigmaat ;
    double *singlepofi, *fwork, *gfac, dipmatsqsum, singleimtau ;
    double ekinetic, kappa, fermifct, t0m1re, t0m1im ;
    double sqrthalf, u[5][5], dtaud[5][5][5][5] ; /* **u, ****dtaud ; */
    double sqq[3][3] ;
    int iorb, *modified, modifiedi, lmax_cont, l ;
// for MPI
    int mpit = 0 ;
    int myid, numproces ;

    int neenv_p ;
    double *eenv_p ;

    MPI_Init(&argc,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD,&numproces) ;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid) ;
  

//    if (myid == 0){


/* read input */
    readinput(absorbergsfile, absorberchfile, 
	     &nshell, corewffile, &eenvmin, &eenvmax,
	     &deenv, &lmax, &lsym, &ksicor, &ksival, &tendqval, 
	     &hz, &vmtzero, &imvmtzero, &efermi, &omegamin, 
	     &omegamax, &domega, 
	     &vmixpar, &rkscale, &e0min, &e0max, &nrwf1,
	     &lvalsh, &nvalel, &withrefl, myid) ;   

    //for MPI
    if (myid == 0){

    printf("Here lvalsh = %d\n", lvalsh) ;

/* construct unitarian matrix for basis trafo among d-orbitals */
/*
    u = calloc2double( 2*lvalsh+1, 2*valsh+1 ) ;
*/

    for ( m1 = 0 ; m1 < 2*lvalsh+1 ; m1++ )
      for ( m2 = 0 ; m2 < 2*lvalsh+1 ; m2++ )
        u[m1][m2] = 0. ;

    u[lvalsh][lvalsh] = 1.;
    sqrthalf = sqrt(0.5) ;
    for ( m1 = 1 ; m1 <= lvalsh ; m1++ ) {
      u[ lvalsh + m1 ][ lvalsh + m1 ] =  sqrthalf ;
      u[ lvalsh + m1 ][ lvalsh - m1 ] =  sqrthalf ;
      u[ lvalsh - m1 ][ lvalsh + m1 ] =  sqrthalf ;
      u[ lvalsh - m1 ][ lvalsh - m1 ] = -sqrthalf ;
    }

/*
    dtaud = calloc4double( 2*lvalsh+1, 2*valsh+1, 2*lvalsh+1 , 2*valsh+1 ) ;
*/


/* read absorber potentials */
    gspot = readlmtopot(absorbergsfile) ;
    /*
    for (i=0; i<gspot.nr; i++) printf(" %lf \n", gspot.v[i] ) ;
    */
    chpot = readlmtopot(absorberchfile) ;
    if (!potcompatible(gspot,chpot)) {
      printf("!potcompatible\n"); exit(1); 
    }    
    pot0 = gspot ;
    iz = pot0.z; nr = pot0.nr; nspin = pot0.nspin; 
    a = pot0.a; b = pot0.b;  wsr = pot0.wsr ;
    vofi = ( double * ) malloc( nr * nspin * sizeof( double ) ) ;
    pot0.v = vofi ; 
    /*
    iz = pot0p->z; nr = pot0p->nr; nspin = pot0p->nspin; 
    a = pot0p->a; b = pot0p->b;  wsr = pot0p->wsr ;
    */
    printf("iz=%d nr=%d nspin=%d a=%lf b=%lf wsr=%lf\n",
	   iz, nr, nspin, a, b, wsr ) ;
    rofi = ( double * ) malloc( nr * sizeof( double ) ) ;
    for (ir = 0; ir < nr; ir++) {
      rofi[ir] = b * ( exp( a * ir ) - 1. ) ;
    }
    singlepofi = ( double * ) calloc( 2 * nr , sizeof( double ) ) ;

/* read Ca 2p wavefunction */
    cofi = readcofi( corewffile, iz, nr, a, wsr, &ncorsh, &lcorsh,
	      &ecore ) ;
    printf(" ecore = %9.4lf, lcorsh = %d, ksicor = %9.4lf \n", 
	   ecore, lcorsh, ksicor ) ;

/* calculate bare spherical Hartree potential from 1 core electron */    
    vcoree = ( double * ) malloc( nr * sizeof( double ) ) ;
    calcvhart(nr, rofi, cofi, vcoree) ;

/* mix potentials */
    for ( i = 0; i < nr * nspin; i++ ) {
      vofi[i] = vmixpar * ( gspot.v[i] - vcoree[ i % nr ] ) 
	       + ( 1. - vmixpar ) * chpot.v[i] ;
    }

/* print out different potentials for check */
    fp = fopen("checkpots","w") ;
    fprintf(fp,"# r, r*gspot.v, r*chpot.v, r*vcoree, r*vofi, r*pot0.v\n" );
    for (ir = 0; ir < nr; ir++) {
      r = rofi[ir] ;
      fprintf(fp," %le %le %le %le %le %le\n", 
	      r, r*gspot.v[ir], r*chpot.v[ir], r*vcoree[ir], 
	      r*vofi[ir], r*pot0.v[ir] ) ;
    }
    fclose(fp) ;

    free( chpot.v ) ;



//    mpit=1;
/*  for MPI */
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Bcast(&mpit, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    printf ("myid = %d mpit= %d\n", myid, mpit);


 
/* read or calculate rho, chi */
    fp = fopen("rhochi.dat","r") ;
    if (fp != NULL) {             /* read rho and chi from file */
      // for MPI
      if (myid == 0){
      for (i=0;i<3;i++) while((c=fgetc(fp))!=EOF&&c!='\n'); /*3 comment lines*/
      fscanf(fp,"%d%d%d%d", &neenv, &ntasc, &nnzele, &lmax_cont);
      printf("Read from rhochi.dat: neenv=%d, ntasc=%d, nnzele=%d, lmax_cont=%d\n", 
	     neenv, ntasc, nnzele, lmax_cont );
      ntas = ( lmax + 1 ) * ( lmax + 1 ) ;
      if (ntasc != ntas || lmax != lmax_cont ) {
	printf("Problem: ntasc == %d != %d == ntas || lmax ==%d != %d== lmax_cont\n",
	ntasc, ntas, lmax, lmax_cont );
	exit(1);
      }
      eenv = calloc1double( neenv ) ;
      rho = calloc4double( ntas, ntas, 2, neenv ) ;
      chi = calloc4double( ntas, ntas, 2, neenv ) ;
      wnwj = calloc2double( LMAX_cont + 1, neenv ) ;
      for ( i = 0 ; i < neenv ; i++ ) {
	fscanf(fp,"%lf",&eenv[i]);
        for ( l = 0 ; l <= LMAX_cont ; l++ ) 
          fscanf(fp,"%lf",&wnwj[l][i] ) ;
	for ( j = 0 ; j < nnzele ; j++ ) {
	  fscanf(fp,"%d%d%lf%lf%lf%lf",&lm1,&lm2,&rerho,&imrho,&rechi,&imchi);
	  rho[lm1][lm2][0][i] = rerho ;
	  rho[lm1][lm2][1][i] = imrho ;
	  chi[lm1][lm2][0][i] = rechi ;
	  chi[lm1][lm2][1][i] = imchi ;
	}
      }
      }//MPI end
      fclose(fp);
    }
    else { /* call environment.c that performs standard M.S. calculation */
      if (myid == 0){
      neenv = ( eenvmax - eenvmin ) / deenv + 1.0001 ;
      eenv = calloc1double( neenv ) ;
      if ( lsym >= 0 ) ntas = 1 ; else ntas = ( lmax + 1 ) * ( lmax + 1 ) ;
      rho = calloc4double( ntas, ntas, 2, neenv ) ;
      chi = calloc4double( ntas, ntas, 2, neenv ) ;
      wnwj = calloc2double( LMAX_cont + 1, neenv ) ;
      

      for ( i = 0 ; i < neenv ; i++) eenv[i] = eenvmin + i*deenv ;
      printf("E[i] = ");
      if (neenv < 4 ) {
	for ( i = 0 ; i < neenv ; i++ ) printf("%10lf ",eenv[i]);
	printf("\n");
      }
      else {
	for ( i = 0 ; i < 3 ; i++ ) printf("%10lf ",eenv[i]);
	printf("... E[%d] = %10lf\n", neenv-1, eenv[neenv-1]);
      }
      printf("Here\n");
      // for MPI
      }
      MPI_Barrier(MPI_COMM_WORLD);
      environment( nr, cofi, pot0, nshell, vmtzero, imvmtzero,
		   neenv, eenv, lsym, ntas, rho, chi, wnwj,
                   numproces, myid, MPI_COMM_WORLD );
      if (myid == 0){
      if ( lsym >= 0 ) exit(1) ;
      } 
    }
    
    /* for MPI */
    if (myid == 0){
   

/* call calcppofi */
    nrwfall = calcppofi(nr, a, b, rofi, vofi, iz, lvalsh, e0min, e0max, 
    			1., 0., &nod0, &e0p, &pr0, &dpr0, &ppofi ); 
    nrwf = (nrwfall+1)/2 + nrwf1 ;
    if ( nrwf > nrwfall ) nrwf = nrwfall ;

    w3jtabmake() ; 


    nopen1  = 0 ;  /* nopen1 = 0 <=> first CLOSED shell, so no reordering */

    reorderrwfs( cofi, ecore, &nrwf, nopen1, &ppofi, &e0p, &pr0, &dpr0 ) ;

    gramschmidtplus(nr, rofi, nrwf, &ppofi, &smat ) ;

    invuppertri( nrwf, smat, &tmat ) ;

    transfe0pr0dpr0(nrwf, smat, tmat, e0p, pr0, dpr0, &e0po, &pr0o, &dpr0o ) ;
    
    e0pomod = calloc2double( nrwf, nrwf ) ;
/*    for ( i = 0 ; i < nrwf ; i++ )
      for ( j = 0 ; j < nrwf ; j++ )
	e0pomod[i][j] = ( i < 2 && j < 2 ) ? e0po[i][j] : 0. ;
*/
/* symmetrize ! */
    for ( i = 0 ; i < nrwf ; i++ )
      for ( j = 0 ; j <= i ; j++ ) {
	dum = ( i < 2 && j < 2 ) ? ( (e0po[i][j] + e0po[i][j])/2. ) : 0. ;
	e0pomod[i][j] = dum ;
        e0pomod[j][i] = dum ;
      }

    makeshells( nrwf, lcorsh, lvalsh, ksicor, ksival, tendqval, 
                &lsh, &sorb1sh, &ksish, &tendqsh ) ;

    nsorbs = sorb1sh[nrwf] ;
    if ( sorb1sh[nrwf] > NBITS ) { 
      printf("Too many radial WF. Reducs BasisEmax\n");
      exit(1) ;
    }
    
    for ( i = 0 ; i < nrwf ; i++ ) 
      printf("i = %2d, lsh[i] = %2d, sorb1sh[i] = %2d\n", i,lsh[i],sorb1sh[i]);

    makeconfs( nrwf, lsh, nvalel, 0, &nconfs, &occ, &ngstbasis ) ;

    makestates( nrwf, lsh, sorb1sh, nconfs, 4*lcorsh+2 + nvalel, occ, ngstbasis,
		&gstbasis) ;


    nhamele = calcham( &hamsparse, ngstbasis, gstbasis, nconfs, nrwf, lsh, 
		       sorb1sh, occ, e0pomod, ksish, tendqsh, hz, nr, rofi, 
                       ppofi, rkscale ) ;

    for ( i = 0 ; i < nconfs ; i++ )
      free( occ[i] ) ;
    free( occ ) ;

    ham = ( double * ) calloc( ngstbasis * ngstbasis, sizeof( double ) ) ;
    eneigval = ( double * ) malloc( ngstbasis * sizeof( double ) ) ;
 
    for ( ist = 0 ; ist < ngstbasis ; ist++ ) 
      for ( k = 0 ; k < hamsparse[ist].n ; k++ ) 
	ham[ ist + ngstbasis * hamsparse[ist].j[k] ] = hamsparse[ist].v[k] ;
 

    info =  solve_dsyev( ngstbasis, ham, eneigval ) ; 
    if ( info ) { printf(" DSYEV : INFO = %d\n", info ) ; exit(1) ; } 


    spamadelete( hamsparse, ngstbasis ) ; 
  
    nlevels = printspec( ngstbasis, eneigval, &ngstdeg, &gstenergy ) ; 

    printf("nlevels = %d\n", nlevels ) ;

    /* store degenerate ground states in gstvec[ngstdeg][ngstbasis] */

    gstvec = ( double ** ) malloc( ngstdeg * sizeof( double * ) ) ;
    for ( k = 0 ; k < ngstdeg ; k++ ) {
      gstvec[k] = ( double * ) malloc( ngstbasis * sizeof( double ) ) ; 
      for ( ist = 0 ; ist < ngstbasis ; ist++ )
	gstvec[k][ist] = ham[ ist + ngstbasis * k ] ;
    }
    printf("gstenergy = %lf.  %d degen. gstvec's:\n", gstenergy, ngstdeg ) ;
    for ( k = 0 ; k < ngstdeg ; k++ ) 
      printstatevector( ngstbasis, gstvec[k], gstbasis ) ;
    
    free( eneigval ) ;
    free( ham ) ;

    valsh = 1 ;
    valshdeg = 4 * lvalsh + 2 ;		

    nvalshholes = valshdeg - nvalel ;

    gstholevecs = calloc3double( ngstdeg, nvalshholes, valshdeg ) ;
    valshdensmat = ( double * ) malloc( valshdeg * valshdeg * sizeof(double) );
    valshdenseval = ( double * ) malloc( valshdeg * sizeof(double) ) ;
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) { 
      for ( ims = 0 ; ims < valshdeg ; ims++ ) { 
        for ( jms = 0 ; jms < valshdeg ; jms++ ) { 
          iorb = sorb1sh[ valsh ] + ims ; 
          jorb = sorb1sh[ valsh ] + jms ; 
          po1plistitem = 0 ;
          o1plistitemadd( &po1plistitem, iorb, jorb, 1. ) ;
          valshdensmatop = o1pmake( po1plistitem ) ; 
          o1plistdelete( po1plistitem ) ;
          nopele = o1ptospama( ngstbasis, gstbasis, ngstbasis, gstbasis, 
                   valshdensmatop, &pvalshdensmatsmline0 ) ;
          o1pdelete( valshdensmatop ) ;
/*   printout spama for check 
          printf("ims = %d, jms=%d, numtotele = %d\n", ims, jms, nopele ) ;
          for ( i = 0 ;  i < ngstbasis ; i++ )
            if ( pvalshdensmatsmline0[i].n > 0 ) {
	      printf("<%3d|valshdensmat| . > . :", i ) ;   
	      spamalinewrite( pvalshdensmatsmline0[i] ) ;
	    }
*/
/* calculate < gst | valshdensmat | gst >  */
          sum = 0. ; 
	  for ( i = 0 ; i < ngstbasis ; i++ ) {	
	    psmline = &pvalshdensmatsmline0[i] ;
	    sumi = 0. ;
	    for ( k = 0 ; k < psmline -> n ; k++ )
	      sumi +=  psmline->v[k]  *  gstvec[igstdeg][ psmline -> j[k] ] ;
	    sum += gstvec[igstdeg][i] * sumi ;
	  }
          spamadelete( pvalshdensmatsmline0, ngstbasis ) ;
          valshdensmat[ ims + valshdeg * jms ] = sum ;
        }
      }
      printf("igstdeg = %d. valshdensmat[][] :\n",igstdeg);
      for ( ims = 0 ; ims < valshdeg ; ims++ ) { 
        for ( jms = 0 ; jms < valshdeg ; jms++ )
          printf("%7.3lf ", valshdensmat[ ims + valshdeg * jms] ) ; 
        printf("\n") ;
      }
      solve_dsyev( valshdeg, valshdensmat, valshdenseval ) ;
      printf("  Eval | Evec") ;  
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) {
	printf("\n%6.3lf |", valshdenseval[ ihole ] ) ;
	for ( ims = 0 ; ims < valshdeg ; ims++ ) {
          gstholevecs[igstdeg][ihole][ims]
	     = valshdensmat[ ims + valshdeg * ihole ] ;
          printf(" %6.3lf", gstholevecs[igstdeg][ihole][ims] ) ;
        }
      }
      printf("\n\n");
    }
    free( valshdenseval ) ;
    free( valshdensmat ) ; 

/* calculate photo-emission eigenstates (N-1) electron system */
    makeconfs( nrwf, lsh, nvalel, -1, &nconfs, &occ, &ncha0basis ) ;
    makestates( nrwf, lsh, sorb1sh, nconfs, 4*lcorsh+1 + nvalel, occ, ncha0basis,
		&cha0basis) ;
      
    nhamele = calcham( &hamsparse, ncha0basis, cha0basis, nconfs, nrwf, lsh, 
		       sorb1sh, occ, e0pomod, ksish, tendqsh, hz, nr, rofi, 
                       ppofi, rkscale ) ;

    for ( i = 0 ; i < nconfs ; i++ )
      free( occ[i] ) ;
    free( occ ) ;

    ham = ( double * ) calloc( ncha0basis * ncha0basis, sizeof( double ) ) ;
    eneigval = ( double * ) malloc( ncha0basis * sizeof( double ) ) ;
    
    for ( ist = 0 ; ist < ncha0basis ; ist++ ) 
      for ( k = 0 ; k < hamsparse[ist].n ; k++ ) 
	ham[ ist + ncha0basis * hamsparse[ist].j[k] ] = hamsparse[ist].v[k] ;

    
    info =  solve_dsyev( ncha0basis, ham, eneigval ) ; 
    if ( info ) { printf(" DSYEV : INFO = %d\n", info ) ; exit(1) ; } 

    spamadelete( hamsparse, ncha0basis ) ; 
  
    nlevels = printspec( ncha0basis, eneigval, &fst0deg, &fst0energy ) ; 

    printf("nlevels = %d\n", nlevels ) ;

    /* store eigen energies in cha0energy[ncha0basis] 
       and eigen states in chavec0[ncha0basis][ncha0basis] */
    cha0energy = ( double * ) malloc( ncha0basis * sizeof( double ) ) ;
    cha0vec = ( double ** ) malloc( ncha0basis * sizeof( double * ) ) ;
    for ( k = 0 ; k < ncha0basis ; k++ ) {
      cha0energy[k] = eneigval[k] ;
      cha0vec[k] = ( double * ) malloc( ncha0basis * sizeof( double ) ) ; 
      for ( ist = 0 ; ist < ncha0basis ; ist++ )
	cha0vec[k][ist] = ham[ ist + ncha0basis * k ] ;
    }
    /*
    for ( k = 0 ; k < ncha0basis ; k++ ) {
      printf("cha0energy = %10.6lf\n", cha0energy[k] ) ;
      printstatevector( ncha0basis, cha0vec[k], cha0basis ) ;
    }
    */
    free( eneigval ) ;
    free( ham ) ;

    contsh1 = 2 ;
    if ( nrwf < contsh1 + 1 ) { 
      printf("main PANIC: nrwf < contsh1 + 1\n"); exit(1); 
    }
    contsorb1 = sorb1sh[contsh1] ;

/* construct chabasis */
    nchabasis = ncha0basis * valshdeg ;
    chabasis = ( struct Fock * ) malloc( nchabasis * sizeof( struct Fock ) ) ;
    for ( ims = 0 ; ims < valshdeg ; ims++ ) {
      i1 = contsorb1 + ims ;
      for ( k = 0 ; k < ncha0basis ; k++ ) {
	statek = cha0basis[k] ;
	sign = fockcre( &statek, i1 );
	chabasis[ ncha0basis * ims + k ] = statek ;
	if ( sign != 1 ) {
	  printf("main PANIC: fockcre( sign != 1 \n"); exit(1); 
	}
      }	
    }
    /*
    for ( k = 0 ; k < nchabasis ; k++ ) {
      printf("\n%2d)  ",k); fockdisplay( &chabasis[k] ) ;
    }
    */

/* From (N-1) electron eigenstates to channel functions */
/* in past :    nchannels = ncha0basis ; */ 
    nchannels = ncha0basis * nvalshholes ;
    chaenergy = ( double * ) malloc( nchannels * sizeof( double ) ) ;
    for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) 
      for ( k1 = 0 ; k1 < ncha0basis ; k1++ )
        chaenergy[ ncha0basis * ihole + k1 ] = cha0energy[k1] ;
  
    chavec = calloc3double( ngstdeg, nchannels, nchabasis ) ;
    for ( k1 = 0 ; k1 < ncha0basis ; k1++ )
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ )
        for ( ims = 0 ; ims < valshdeg ; ims++ ) 
          for ( k2 = 0 ; k2 < ncha0basis ; k2++ )  
	    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )
              chavec[igstdeg][ ncha0basis*ihole + k1 ][ ncha0basis*ims + k2 ] 
                = cha0vec[k1][k2] * gstholevecs[igstdeg][ihole][ims] ;
   
/* construct N electron final state basis "state" */
    makeconfs( nrwf, lsh, nvalel, 1, &nconfs, &occ, &nfstbasis ) ;
    makestates( nrwf, lsh, sorb1sh, nconfs, 4*lcorsh+2 + nvalel, occ, 
		nfstbasis, &fstbasis ) ;


/* construct transformation matrix urestrfst between original and 
   restricted final state basis:  < fstbasis(i) | fstresbas(j) > .  
   For each i, there will be at most one non-zero j and we have
   < fstbasis(i) | fstresbas(j) > = gstholevec(ims)  where
   ims is the occupied orbital in the continuum shell of | fstbasis(i) >.
*/       

/* first a little check */
    for ( ish = contsh1 ; ish < nrwf ; ish++ ) 
      if ( sorb1sh[ish+1] - sorb1sh[ish] != valshdeg ) {
        printf("main PANIC: sorb1sh[] vs valshdeg\n" ) ; exit(1) ; } 
    contsorb1 = sorb1sh[contsh1] ;

/* count bound basis states and continuum basis states and 
   check that all bound states come before continuum states 
*/
    nfstboundbas = 0 ;
    nfstcontbas  = 0 ;
    for ( ist = 0 ; ist < nfstbasis ; ist++ ) { 
      iorb = contsorb1 ; 
      while ( iorb < nsorbs && fockocc( &fstbasis[ist], iorb ) == 0 ) 
        iorb++ ; 
      if ( iorb < nsorbs )
        nfstcontbas++ ; 
      else {               /* then, it is a "bound" state (3d^{n+1}) */
        nfstboundbas++ ;  
/*        if ( nfstcontbas != 0 ) { 
          printf("some cont.states come before bound states\n"); exit(1) ; }
*/
      }
    }
    if ( nfstbasis != nfstboundbas + nfstcontbas ) {
      printf("nfstbasis != nfstboundbas + nfstcontbas\n"); exit(1) ; } 
    printf("nfstboundbas=%d, nfstcontbas=%d\n",nfstboundbas,nfstcontbas ) ;

/*    printf("\nWARNING: so far only works for 3d0 and 3d9 conf\n\n") ; */

    ncontshs = nrwf - contsh1 ;
    printf(" ncontshs = %d\n", ncontshs ) ;
    nfstresbas = nfstboundbas + nchannels * ncontshs ;
    fstbas2res = calloc3double( ngstdeg, nfstresbas, nfstbasis ) ; 

    for ( ist = 0 ; ist < nfstboundbas ; ist++ ) 
      for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )
        fstbas2res[igstdeg][ist][ist] = 1. ;

    for ( ist = nfstboundbas ; ist < nfstbasis ; ist++ ) {
      statek = fstbasis[ist] ;	
      iorb = contsorb1 ;
      while ( iorb < nsorbs && fockocc( &fstbasis[ist], iorb ) == 0 )
        iorb++ ;
      if ( iorb == nsorbs ) { printf("PANIC: iorb == nsorbs\n");exit(1); }
      sign = fockdes( &statek, iorb ) ;
      if ( sign != 1 ){ printf("mainPANIC:fockdes,sign!=1\n");exit(1);} 
      if ( findstate( ncha0basis, cha0basis, &statek, &k ) == 0 ) { 
        printf("main PANIC: findstate==0\n") ; exit(1) ; }
      icontsh = ( iorb - contsorb1 ) / valshdeg ;
      ims = ( iorb - contsorb1 ) % valshdeg ;
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) {
        ifstresbas = nfstboundbas + ncha0basis*(nvalshholes*icontsh+ihole) + k;
        if ( ifstresbas >= nfstresbas ){
	  printf("PANIC: ifstresbas >= nfstresbas\n"); exit(1); }
        for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
	  if ( fstbas2res[igstdeg][ifstresbas][ist] != 0. ) {
	    printf("There is a big big mess 2...\n") ; exit(1) ; }
          fstbas2res[igstdeg][ifstresbas][ist] 
	    = gstholevecs[igstdeg][ihole][ims] ;
        }
      }
    }

/* printout trafo matrix for check 
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
      printf("deg. gstate No %d\n", igstdeg ) ;
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) {
	printf(" | ifstresbas > = | %d > \n", ifstresbas ) ;
	for ( ist = 0 ; ist < nfstbasis ; ist++ ) 
	  if ( fabs(fstbas2res[igstdeg][ifstresbas][ist]) > 1.e-6 ) { 
            printf("%.3lf ",fstbas2res[igstdeg][ifstresbas][ist] ) ;
            fockdisplay( &fstbasis[ist] ) ; printf("\n") ;
          }
        }
    }
*/


/* calculate overlap matrix between channel basis x |lms> and 
   N-electron final state basis */
    
    po1plistitem = 0 ;
    nwoplistele = wop2list( nrwf, lsh, sorb1sh, contsh1, pr0o, 
			    &po1plistitem) ;
    printf("nwoplistele = %d.\n", nwoplistele ) ;
    
    wop = o1pmake( po1plistitem ) ;
    o1pprint( wop ) ;
    o1plistdelete( po1plistitem ) ;

    nopele = o1ptospama( nchabasis, chabasis, nfstbasis, fstbasis,
			 wop, &pwopsmline0 ) ;
    printf("numtotele = %d\n", nopele ) ;

/*
    for ( i = 0 ;  i < nfstbasis ; i++ ) {
      if ( pwopsmline0[i].n > 0 ) {
	printf("<%3d|Wop| . > . :", i ) ;   
	spamalinewrite( pwopsmline0[i] ) ;
      }
      else {
	printf(" %3d", i) ;
      }
    }
*/
/* calculate  Wop  for the restricted final state basis (fstresbas)
   <fres_I || chabasis_k > = \sum_k < fres_I | f_J > < f_J || chabasis_k > 
   where < f_J | fres_I > = fstbas2res[igstbas][I][J]
   ATTENTION: code only correct for real trafo matrix fstbas2res
*/
/* VERY EXPERIMENTAL */
    if ( ngstdeg == 1 && nfstresbas == nfstbasis ) {
      pwopressml = ( struct Spamaline ** )
        malloc( ngstdeg * sizeof( struct Spamaline * ) ) ;  
      for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )
        pwopressml[igstdeg] = ( struct Spamaline * ) 
          malloc( nfstresbas * sizeof( struct Spamaline ) ) ; 
      pwopressml[0] = pwopsmline0 ;
    }
    else {
      printf("NOT ngstdeg == 1 && nfstresbas == nfstbasis\n") ;
      exit(1) ;
    }
/*
    pwopressml = ( struct Spamaline ** )
      malloc( ngstdeg * sizeof( struct Spamaline * ) ) ;  
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )
      pwopressml[igstdeg] = ( struct Spamaline * ) 
        malloc( nfstresbas * sizeof( struct Spamaline ) ) ; 
    dumvec = ( double * ) malloc( nchabasis * sizeof( double ) ) ;
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) {
        printf("%4d",ifstresbas);if(ifstresbas%20==0)printf("\n");
        for ( jchabas = 0 ; jchabas < nchabasis ; jchabas++ ) { 
          sum = 0. ;	
          for ( ist = 0 ; ist < nfstbasis ; ist++ ) {
            pwopsml = pwopsmline0 + ist ;
            for ( k = 0 ; k < pwopsml->n && pwopsml->j[k] != jchabas ; k++ ) ;
            if ( pwopsml->j[k] == jchabas )
              sum += fstbas2res[igstdeg][ifstresbas][ist] * pwopsml->v[k] ;
          }
          dumvec[ jchabas ] = sum ;
        }
        darray2spamaline( &pwopressml[igstdeg][ifstresbas], nchabasis, dumvec );
      }
    }
    free( dumvec ) ;
    spamadelete( pwopsmline0, nfstbasis ) ; 
*/
/*  EXPERIMENTAL another way. Could work with a different chavec ...
    pwopressml = ( struct Spamaline ** )
      malloc( ngstdeg * sizeof( struct Spamaline * ) ) ;  
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )
      pwopressml[igstdeg] = ( struct Spamaline * ) 
        malloc( nfstresbas * sizeof( struct Spamaline ) ) ; 
    dumvec = ( double * ) malloc( nchabasis * sizeof( double ) ) ;

    dummat  = calloc2double( nfstbasis, nchannels ) ;
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) {
        for ( ist = 0 ; ist < nfstbasis ; ist++ ) {
          pwopsml = pwopsmline0 + ist ; 
          sum = 0. ;     
          for ( k = 0 ; k < pwopsml->n ; k++ ) {
*/
 /*  we have:  pwopsml->j[k] == ichabas == ncha0basis * ims + ich0  */
/*            ims = pwopsml->j[k] / ncha0basis ; 
            sum += pwopsml->v[k] * gstholevecs[igstdeg][ihole][ims] ;
          }
          for ( ich0 = 0 ; ich0 < ncha0basis ; ich0++ ) {
            ich = ncha0basis *  ihole + ich0 ;
            dummat[ist][ich] = sum ; 
          }
        }
      }
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
        for ( ich = 0 ; ich < nchannels ; ich++ ) {
          sum = 0. ;
          for ( ist = 0 ; ist < nfstbasis ; ist++ ) 
            sum += fstbas2res[igstdeg][ifstresbas][ist] * dummat[ist][ich] ;
          dumvec[ ich ] = sum ;
        }
        darray2spamaline( &pwopressml[igstdeg][ifstresbas], nchannels, dumvec);
      }
    }
    free2double( nfstbasis, dummat ) ;
    free( dumvec ) ;
    spamadelete( pwopsmline0, nfstbasis ) ; 
*/      	

/*
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) { 
      printf("gst No %d\n", igstdeg ) ;
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
	if ( pwopressml[igstdeg][ifstresbas].n > 0 ) {
	  printf("<%3d|| . > . :", ifstresbas ) ;   
	  spamalinewrite( pwopressml[igstdeg][ifstresbas] ) ;
	}
	else {
	  printf(" %3d", ifstresbas ) ;
	}
      }
    }
*/
	

/* construct dipole operator  P^(1)_q = r C^(1)_q , q=-1,0,1 [Cowan (14.24)] */

    inishell = 0 ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      q = iq - 1 ;
      po1plistitem = 0 ;
      ndiplistele = dipole2list( nrwf, lsh, sorb1sh, inishell, q, 
				 nr, rofi, ppofi, &po1plistitem) ;
      printf("q = %d. ndiplistele = %d.\n", q, ndiplistele ) ;

      dipop[iq] = o1pmake( po1plistitem ) ;
      o1plistdelete( po1plistitem ) ;
      nopele = o1ptospama( ngstbasis, gstbasis, nfstbasis, fstbasis, 
			   dipop[iq], &pdipsmline0[iq] ) ;
      printf("numtotele = %d\n", nopele ) ;
    }
/* printout spama for check 
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      q = iq - 1 ;
      for ( i = 0 ;  i < nfstbasis ; i++ ) {
	if ( pdipsmline0[iq][i].n > 0 ) {
	  printf("<%3d|r_{%2d}| . > . :", i, q ) ;   
	  spamalinewrite( pdipsmline0[iq][i] ) ;
	}
	else {
	  printf(" %3d", i) ;
	}
      }
    }
*/


/* calculate dipole matrix elements between final state basis and gstates 
   < f_I | D[iq] | gst[igstdeg] >  for all final basis states f_I = state[I] 
   | gst[igstdeg] > = \sum_J | gstbasis[J] > * gstvec[igstdeg][J],   
   igstdeg = 0,ngstdeg-1,  J = 0,ngstbasis-1, I = 0,nfstbasis-1, iq = 0,1,2. 
*/
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      dipmatfbas[iq] = calloc2double( ngstdeg, nfstbasis ) ; /* mem. alloc */
      for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
	for ( i = 0 ; i < nfstbasis ; i++ ) {	
	  psmline = &pdipsmline0[iq][i] ;
	  sum = 0. ;
	  for ( k = 0 ; k < psmline -> n ; k++ ) {
	    sum  +=  psmline -> v[k]  *  gstvec[igstdeg][ psmline -> j[k] ] ;
	  }
	  dipmatfbas[iq][igstdeg][i] = sum ;
	}
      }
    }

/* caculate dipole matrix elements between fstresbas and gststates  
   < fres_I | D | gst >  = \sum_J < fres_I | f_J > < f_J | D | gst >
   where < f_J | fres_I > = fstbas2res[igstbas][I][J]
   ATTENTION: code only correct for real trafo matrix fstbas2res
*/
    for ( iq = 0 ; iq < 3 ; iq++ ) { 
      dipmatfresbas[iq] = calloc2double( ngstdeg, nfstresbas ) ; 
      for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) { 
        for ( i = 0 ; i < nfstresbas ; i++ ) { 
	  sum = 0. ;
	  for ( j = 0 ; j < nfstbasis ; j++ ) 
	    sum += fstbas2res[igstdeg][i][j] * dipmatfbas[iq][igstdeg][j] ; 
	  dipmatfresbas[iq][igstdeg][i] = sum ;
        }
      } 
    } 

/*
    printf("dipmatfresbas. cols j = 3 * igstdeg + iq, lines i = Index fst basis\n");
    for ( i = 0 ; i < nfstresbas ; i++ ) {	
      for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
	for ( iq = 0 ; iq < 3 ; iq++ ) 
	  printf(" %7.3lf", dipmatfresbas[iq][igstdeg][i] ) ;
	printf("\t") ;
      }
      printf("\n");
    }
*/
    
    nstates = nfstresbas ;	

    modified = ( int * ) calloc( nchannels, sizeof( int ) ) ;
    kga    = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    rkbes  = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    drkbes = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    rkneu  = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    drkneu = ( double * ) calloc( nchannels, sizeof( double ) ) ;


    alpha = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    beta  = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    wdagger = ( double * ) calloc( nchannels * nchannels, sizeof( double ) ) ;
    normfac = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    mat1 = ( double * ) calloc( nchannels * nchannels, sizeof( double ) ) ; 
    mat1linei = ( double * ) calloc( nchannels, sizeof( double ) ) ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      dipmatchabas[iq] = calloc2double( ngstdeg, nchannels ) ;
      dipmatele[iq]  = calloc2double( ngstdeg, nchannels ) ;
    }
    tau = ( double * ) calloc( nchannels * nchannels * 2, sizeof( double ) ) ;

    jk = ( int * ) calloc( nstates, sizeof( int ) ) ;
 
    lambda  = ( double * ) calloc( nstates * nstates, sizeof( double ) ) ;
    if ( lambda == 0 ) { printf("lambda alloc failed\n"); exit(1);}
    qmetric = ( double * ) calloc( nstates * nstates, sizeof( double ) ) ;
    if ( qmetric == 0 ) { printf("qmetric alloc failed\n"); exit(1);}


    alphalong = ( double * ) calloc( nstates, sizeof( double ) ) ;
    betalong  = ( double * ) calloc( nstates, sizeof( double ) ) ;
    revec = ( double * ) calloc( nstates * nstates, sizeof( double ) ) ;
    normfac = ( double * ) calloc ( nchannels, sizeof( double ) ) ;


    fp = fopen("spectrum.dat","w");
    fprintf(fp,"# omega, Ix+Iy+Iz, Ix, Iy, Iz, I-, I+, omega+ecore\n");
    fpq = fopen("specqqprime.dat","w");
    fprintf(fpq,"# omega, I_qq' (qq'= -1-1 -10 -11 0-1 00 01 1-1 10 11)\n");
    fpm = fopen("specmmprime.dat","w");
    fprintf(fpm,"# omega, I_mm' (mm'= -2-2 -2-1 -20 -21 -22 -1-2 .. .. 22) [-2=ixy -1=iyz 0=z2 1=zx 2=x2]\n");

    /* this is actually H+L */
/*
    printf("\n h+L\n") ;
    for ( i = 0 ; i < nrwf ; i++ ) {
      for ( j = 0 ; j < nrwf ; j++ ) {
	e0pomod[i][j] = e0po[i][j] + pr0o[i]*dpr0o[j] ;
	printf(" %10.6lf", e0pomod[i][j] ) ;
      }
      printf("\n") ;
    }
*/
    printf("\n h+L symmetrized !! \n") ;
    for ( i = 0 ; i < nrwf ; i++ ) {
      for ( j = 0 ; j < nrwf ; j++ ) {
	e0pomod[i][j] = 
          (e0po[i][j] + pr0o[i]*dpr0o[j] + e0po[j][i] + pr0o[j]*dpr0o[i]) / 2.;
	printf(" %10.6lf", e0pomod[i][j] ) ;
      }
      printf("\n") ;
    }

    nhamele = calcham( &hamsparse, nfstbasis, fstbasis, nconfs, nrwf, lsh, 
		       sorb1sh, occ, e0pomod, ksish, tendqsh, hz, nr, rofi, 
                       ppofi, rkscale ) ;

    /* this is actually Q  symmetrized !! */
    printf("\n Q\n") ;
    for ( i = 0 ; i < nrwf ; i++ ) {
      for ( j = 0 ; j < nrwf ; j++ ) {
	dum = ( pr0o[i]*pr0o[j] + pr0o[j]*pr0o[i] ) / 2. ;
	e0pomod[i][j]  = ( fabs( dum ) > QMETRICZERO ? dum : 0. ) ;
	printf(" %10.6lf", e0pomod[i][j] ) ;
      }
      printf("\n") ;
    }    

    dummat  = calloc2double( nfstbasis, nfstbasis ) ;
    dummattimesu = calloc2double( nfstbasis, nfstresbas ) ;

    for ( ist = 0 ; ist < nfstbasis ; ist++ )
      for ( k = 0 ; k < hamsparse[ist].n ; k++ ) 
	dummat[ ist ][ hamsparse[ist].j[k] ] 
	  = hamsparse[ist].v[k] ;
    spamadelete( hamsparse, nfstbasis ) ;

    hamressparse = ( struct Spamaline ** ) 
      malloc( ngstdeg * sizeof( struct Spamaline * ) ) ;
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )  
      hamressparse[igstdeg] = ( struct Spamaline * )
        malloc( nfstresbas * sizeof( struct Spamaline ) ) ; 
    dumvec = ( double * ) malloc( nfstresbas * sizeof( double ) ) ;


   printf("here\n") ;

/* caculate hamiltonian in fstresbas 
   < fres_K | H | fres_L >  = \sum_{IJ} 
      < fres_K | f_I > < f_I | H | f_J > < f_J | fres_L>
   where < f_J | fres_I > = fstbas2res[igstdeg][I][J]
   ATTENTION: code only correct for real trafo matrix fstbas2res
*/
    printf("nfstbasis=%4d, nfstresbas=%4d\n",nfstbasis,nfstresbas); 
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
      for ( i = 0 ; i < nfstbasis ; i++ ) 
        for ( jfstresbas = 0 ; jfstresbas < nfstresbas ; jfstresbas++ ) {
	  sum = 0. ; 
          for ( j = 0 ; j < nfstbasis ; j++ )
            sum += dummat[i][j] * fstbas2res[igstdeg][jfstresbas][j] ;
          dummattimesu[i][jfstresbas] = sum ; 
        }
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
        for ( jfstresbas = 0 ; jfstresbas < nfstresbas ; jfstresbas++ ) { 
          sum = 0. ;  
          for ( i = 0 ; i < nfstbasis ; i++ )
            sum += fstbas2res[igstdeg][ifstresbas][i] 
                    * dummattimesu[i][jfstresbas] ;  
          dumvec[jfstresbas] = sum ;
        }
        darray2spamaline( &hamressparse[igstdeg][ifstresbas],  
                          nfstresbas, dumvec ); 
      }
    } 



/* printout spama for check */
/*
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) { 
      printf("\n hamilton matrix in fstresbas. gstate No %d\n\n", igstdeg ) ;
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
	if ( hamressparse[igstdeg][ifstresbas].n > 0 ) {
	  printf("<%3d|H| . > . : ", ifstresbas ) ;   
	  spamalinewrite( hamressparse[igstdeg][ifstresbas] ) ;
	}
	else {
	  printf(" %3d", ifstresbas ) ;
	}
      }
    }
*/

    nqmetricele = calcham1p( &qmetricsparse, nfstbasis, fstbasis, nconfs, nrwf, 
			     lsh, sorb1sh, occ, e0pomod, nr, rofi, ppofi ) ;

    for ( i = 0 ; i < nconfs ; i++ )
      free( occ[i] ) ;
    free( occ ) ;

    for ( i = 0 ; i < nfstbasis ; i++ ) 
      for ( j = 0 ; j < nfstbasis ; j++ ) 
        dummat[i][j] = 0. ;

    for ( ist = 0 ; ist < nfstbasis ; ist++ ) 
      for ( k = 0 ; k < qmetricsparse[ist].n ; k++ ) 
	dummat[ ist ][ qmetricsparse[ist].j[k] ] 
	  = qmetricsparse[ist].v[k] ;
    spamadelete( qmetricsparse, nfstbasis ) ;

    qressparse = ( struct Spamaline ** ) 
      malloc( ngstdeg * sizeof( struct Spamaline * ) ) ;
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ )  
      qressparse[igstdeg] = ( struct Spamaline * )
        malloc( nfstresbas * sizeof( struct Spamaline ) ) ; 


/* caculate qmetric in fstresbas  */

    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) {
      for ( i = 0 ; i < nfstbasis ; i++ ) 
        for ( jfstresbas = 0 ; jfstresbas < nfstresbas ; jfstresbas++ ) {
	  sum = 0. ; 
          for ( j = 0 ; j < nfstbasis ; j++ )
            sum += dummat[i][j] * fstbas2res[igstdeg][jfstresbas][j] ;
          dummattimesu[i][jfstresbas] = sum ; 
        }
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
        for ( jfstresbas = 0 ; jfstresbas < nfstresbas ; jfstresbas++ ) { 
          sum = 0. ;  
          for ( i = 0 ; i < nfstbasis ; i++ )
            sum += fstbas2res[igstdeg][ifstresbas][i] 
                    * dummattimesu[i][jfstresbas] ;  
          dumvec[jfstresbas] = sum ;
        }
        darray2spamaline( &qressparse[igstdeg][ifstresbas],  
                          nfstresbas, dumvec ); 
      }
    } 

    free( dumvec ) ;
    free2double( nfstbasis, dummattimesu ) ; 
    free2double( nfstbasis, dummat ) ;

/* printout spama for check */
/*
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) { 
      printf("\n Q-matrix in fstresbas. gstate No %d\n\n", igstdeg ) ;
      for ( ifstresbas = 0 ; ifstresbas < nfstresbas ; ifstresbas++ ) { 
	if ( qressparse[igstdeg][ifstresbas].n > 0 ) {
	  printf("<%3d|Q| . > . : ", ifstresbas ) ;   
	  spamalinewrite( qressparse[igstdeg][ifstresbas] ) ;
	}
	else {
	  printf(" %3d", ifstresbas ) ;
	}
      }
    }
*/

    for ( iq = 0 ; iq < 3 ; iq++ )
      for ( iqp = 0 ; iqp < 3 ; iqp++ )
        sigmaqqgst[iq][iqp] = ( double * ) malloc( ngstdeg * sizeof( double ) );

/*
    for ( i=0 ; i < neenv ; i++ ) {
      printf("\n E=%lf ", eenv[i] ) ;
      for ( l = 0 ; l <= LMAX_cont ; l++ )
        printf(" %lf", wnwj[l][i] ) ;
    }
    printf("\n") ;
*/
/* (total final state) energy loop */
    for (omega = omegamin ; omega <= omegamax ; omega += domega) 
    {
    for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) 
    {
      for ( i = 0 ; i < nchannels ; i++ ) {
	ekinetic = gstenergy + omega - chaenergy[i] - vmtzero ;
	if ( ekinetic >= 0. ) {
	  realk = sqrt( ekinetic ) ;
	  modified[i] = 0 ;
	} 
	else {
	  realk = sqrt( -ekinetic ) ;
	  modified[i] = 1 ;
	}
	kga[i] = realk ;
	x = realk * wsr ;
	besneu( modified[i], x, lvalsh, &jl, &nl, &jlp, &nlp ) ;
	rkbes[i]  = x * jl ;
	drkbes[i] = realk * jl + realk * x * jlp ;	
	rkneu[i]  = x * nl ;
	drkneu[i] = realk * nl + realk * x * nlp ;	  
      }

      if ( omega < omegamin + 1.e-5 || omega > omegamax - 1.e-5) {
        printf("\nomega=%lf, omin=%lf, omax=%lf, do=%lf, kga[0]=%lf ",
	     omega, omegamin, omegamax, domega, kga[0]) ;
        printf("%s\n",( modified[0] ? "(imaginary)" : "" ) ) ;
      }	else printf(".") ;
   
/* define qmetric and lambda matrix */

    for ( i = 0 ; i < nstates * nstates ; i++ ) {
      qmetric[i] = 0. ;
      lambda[i] = 0. ;
    }

    for ( ist = 0 ; ist < nstates ; ist++ ) 
      for ( k = 0 ; k < qressparse[igstdeg][ist].n ; k++ )
	qmetric[ ist + nstates * qressparse[igstdeg][ist].j[k] ] 
	  = qressparse[igstdeg][ist].v[k] ;

    for ( ist = 0 ; ist < nstates ; ist++ ) {
      for ( k = 0 ; k < hamressparse[igstdeg][ist].n ; k++ ) 
	lambda[ ist + nstates * hamressparse[igstdeg][ist].j[k] ] 
	  = - hamressparse[igstdeg][ist].v[k] ;
      lambda[ ist + nstates * ist ] += gstenergy + omega ;  
    }

/*
    if ( nstates >= 10 ) {
      printf("Lambda matrix ( upper 10x10 corner)\n") ;
      for ( i = 0 ; i < 10 ; i++ ) {
	for ( j = 0 ; j < 10 ; j++ ) 
	  printf("%lf ", lambda[  i + nstates * j ] ) ;
	printf("\n") ;
      }
      printf("Q matrix ( lower 10x10 corner)\n") ;
      for ( i = nstates -10 ; i < nstates ; i++ ) {
	for ( j = nstates - 10 ; j < nstates ; j++ ) 
	  printf("%lf ", qmetric[  i + nstates * j ] ) ;
	printf("\n") ;
      }
    }
*/
      solve_genev(nstates, lambda, qmetric, alphalong, betalong, revec) ;
/*
   c(p)_{ik} = revec[ i + matdim * k ] 
   here, k labels the "physical" generalized eigenvalues, which are those
   for which betalong[j] != 0., i.e. alphalong[j]/betalong[j] is finite
*/
      k = 0 ;
      for ( j = 0 ;  j < nstates ; j++ ) 
	if ( fabs( betalong[j] ) > BETAZERO )  
	  jk[k++] = j ;
      if ( k != nchannels ) {
	printf(" main PANIC: solution of genev has rank k != nchannels\n") ;
	printf(" nchannels=%d, k=%d -- may change constant BETAZERO\n", 
	       nchannels, k ) ;
	exit(1) ; 
      }

      for ( k = 0 ; k < nchannels ; k++ ) {
	j = jk[k] ;
	alpha[k] = alphalong[j] ;
	beta[k]  = betalong[j]  ;
	sum = 0.;
	for ( ich = 0 ; ich < nchannels ; ich++ ) {
	  dum = spamarealmatele( nchabasis, chavec[igstdeg][ich], nstates, 
				 revec + nstates * j, pwopressml[igstdeg] ) ;
	  wdagger[ k + nchannels * ich ] = dum ;
	  sum += dum * dum ;
	}
	normfac[ k ] = 1. / sqrt( sum ) ;
	for ( ich = 0 ; ich < nchannels ;  ich++ ) 
	  wdagger[ k + nchannels * ich ] *= normfac[ k ] ; /* normalise */
      }


/* normalize   c(p)_{ik} = revec[ i + nstates * jk[k] ] */
      for ( k = 0 ; k < nchannels ; k++ ) 
	for ( i = 0 ; i < nstates ;  i++ ) 
	  revec[ i + nstates * jk[k] ] *= normfac[ k ] ; /* normalise */
/* check 
      for ( k = 0 ; k < nchannels ; k++ ) {
	sum = 0. ;
	for ( i = 0 ; i < nstates ;  i++ ) 
	  sum += revec[ i + nstates * jk[k] ] * revec[ i + nstates * jk[k] ] ;
	printf("%3d %lf\n", k, sum) ;
      }
*/
 
/*
printf("WDAGGER\n") ;
      for ( ich = 0 ; ich < nchannels ; ich++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ )  
	  printf("%c", fabs(wdagger[ ich + nchannels * jch ]) 
                        > 1.e-5 != 0. ? '*' : ' ' ) ;
        printf("\n") ;
      }
      for ( ich = 0 ; ich < nchannels ; ich++ ) 
        for ( jch = 0 ; jch < nchannels ; jch++ )   
          if ( fabs(wdagger[ ich + nchannels * jch ]) > 1.e-5 )
            printf("%2d %2d %lf\n", ich,jch,wdagger[ ich + nchannels * jch ]) ;
*/
/* check that Wdagger * W = 1        
      for ( ich = 0 ; ich < nchannels ; ich++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ ) {
          sum = 0. ;  
	  for ( k = 0 ; k < nchannels ; k++ )  
            sum += wdagger[ich+nchannels*k] * wdagger[jch+nchannels*k] ;
          printf("%.0lf ",sum) ;
	}
        printf("\n");
      }
*/
            
/* calculate the  t^{-1}  matrix, or more precisely the matrix X, which is 
   related to t^{-1} by :   t^{-1} = +/-iK - K X.   (outgoing/incoming BC)
   This is seen as follows (outgoing BC):
   t^{-1} = iK ( beta Wdagger Jdot + alpha Wdagger J )^{-1} 
             * ( beta Wdagger Hdot + alpha Wdagger H ) = iK - K X, 
    where  X = beta Wdagger Jdot + alpha Wdagger J )^{-1}
	   * ( beta Wdagger Ndot + alpha Wdagger N )
   we define
   mat1 := [ alpha Wdagger rkbes + beta  Wdagger drkbes ] 
   mat2 := [ alpha Wdagger rkneu + beta  Wdagger drkneu ]
   and so,  mat1  <-   X =  mat1^{-1} * mat2
   The code is written such that only for mat1 storage is needed.
*/
    for ( i = 0 ; i < nchannels ; i++ ) 
	for ( j = 0 ; j < nchannels ; j++ ) {
 	    wij = wdagger[ i + nchannels * j ] ;
	    mat1[ i + nchannels * j ] = 
		alpha[i] * wij * rkbes[ j ] + beta[i] * wij * drkbes[ j ] ;
	}
/*
printf("MAT1\n") ;
      for ( ich = 0 ; ich < nchannels ; ich++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ )  
	  printf("%c", fabs(mat1[ ich + nchannels * jch ]) 
                        > 1.e-5 != 0. ? '*' : ' ' ) ;
        printf("\n") ;
      }
*/

/* invert mat1 => mat1 = [ alpha Wdagger rkbes + beta  Wdagger drkbes ]^{-1} */

    inv_remat( nchannels, mat1 ) ;
/*
printf("MAT1\n") ;
      for ( ich = 0 ; ich < nchannels ; ich++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ )  
	  printf("%c", fabs(mat1[ ich + nchannels * jch ]) 
                        > 1.e-5 != 0. ? '*' : ' ' ) ;
        printf("\n") ;
      }
*/

/* mat1 = mat1 * mat2  ( on output mat1 contains X, and (tK)^{-1} = -X + i ) */
    for ( i = 0 ; i < nchannels ; i++ ) {
      /* store line i of mat1 in mat1linei */
      for ( j = 0 ; j < nchannels ; j++ ) 
	mat1linei[j] = mat1[ i + nchannels * j ] ; 
      for ( j = 0 ; j < nchannels ; j++ ) {
	sum = 0. ;
	for ( k = 0 ; k < nchannels ; k++ ) {
	  wkj = wdagger[ k + nchannels * j ] ;
	  mat2kj = 
	    alpha[k] * wkj * rkneu[ j ] + beta[k] * wkj * drkneu[ j ] ;
	  sum += mat1linei[k] * mat2kj ;
	}
	mat1[ i + nchannels * j ] = sum ;
      }
    }


/*            calculation of  tau matrix               */
/*                                                     */  
/* define inv. atomic t-matrix  = diag(k) * ( - X +/- i ) and write to tau */
/*                                                 ^ outgoing/incoming BC  */
    for ( ich = 0 ; ich < nchannels ; ich++ ) {
      kgai = kga[ich] ;
      modifiedi = modified[ich] ;
      for ( jch = 0 ; jch < nchannels ; jch++ ) {
	xij  = mat1[ ich + nchannels * jch ] ;
	deltaij = ( ich == jch ? 1. : 0. ) ;
	ijre = 2 * ( ich  +  nchannels * jch ) ;
	tau[ijre  ] = - kgai * xij ;
	if ( modifiedi ) {
	  tau[ijre  ] += - kgai * deltaij ;   /* incoming BC */
	  tau[ijre+1]  = 0. ;
	}
	else {
	  tau[ijre+1]  = - kgai * deltaij ;   /* incoming BC */
	}  /* for outgoing BC change  deltaij  to  -deltaij  in above expressions*/
      }
    }

/*
    for ( ich = 0 ; ich < nchannels ; ich++ ) 
      for ( jch = 0 ; jch < nchannels ; jch++ ) {
	ijre = 2 * ( ich  +  nchannels * jch ) ;
	if ( fabs(tau[ijre])>1.e-6 || fabs(tau[ijre+1])>1.e-6 )
	  printf("tau[%3d,%3d] =  %10lf + %10lf i\n", 
		 ich, jch, tau[ijre], tau[ijre+1] ) ;
      } 
      for ( ich = 0 ; ich < nchannels ; ich++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ ) {
	  ijre = 2 * ( ich  +  nchannels * jch ) ; 
	  printf("%c", fabs(tau[ijre])>1.e-6 || 
	               fabs(tau[ijre+1])>1.e-6  ? '*' : ' ' ) ;
        } 
        printf("\n") ;
      }
*/

/*
    for ( ich0 = 0 ; ich0 < ncha0basis ; ich0++ ) {
      eout = gstenergy + omega - cha0energy[ich0] ;
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) {
        ich = ncha0basis * ihole + ich0 ;
        ijre = 2 * ( ich + nchannels * ich ) ;
        printf("%3d %3d %3d %lf\t%lf\t%lf\n", 
            ich0, ihole, ich, eout, tau[ijre], tau[ijre+1] ) ;
      }
    }  
*/

/* perform tau -= rho,  (with  tau = t^{-1} on input) */
    lmmax = lmorder(lvalsh,lvalsh) ;
    if ( lmmax >= ntas ) { printf("lmmax >= ntas\n") ; exit(1) ;}
    for ( ich0 = 0 ; ich0 < ncha0basis ; ich0++ ) {
      kgai0 = kga[ich0] ;
      eout = gstenergy + omega - cha0energy[ich0] ;
      for ( ihole = 0 ; ihole < nvalshholes ; ihole++ ) {
        ich = ncha0basis * ihole + ich0 ;
        for ( jhole = 0 ; jhole < nvalshholes ; jhole++ ) {
	  jch = ncha0basis * jhole + ich0 ;
          sumre = 0. ;
          sumim = 0. ;
          if ( withrefl == 1 ) {
            for ( m1 = -lvalsh ; m1 <= lvalsh ; m1++ ) {
              lm1 = lmorder( lvalsh, m1 ) ;
              for ( m2 = -lvalsh ; m2 <= lvalsh ; m2++ ) {
                lm2 = lmorder( lvalsh, m2 ) ;
                interpolar( neenv, eenv, rho[lm1][lm2][0], 
          	      1, &eout, &rerho ) ;
                interpolar( neenv, eenv, rho[lm1][lm2][1],
          	      1, &eout, &imrho ) ;
                dum = 0. ;
                for ( ispin = 0 ; ispin < 2 ; ispin++ ) {
                  ims1 = 2 * ( lvalsh + m1 ) + ispin ;
                  ims2 = 2 * ( lvalsh + m2 ) + ispin ;
                  dum += gstholevecs[igstdeg][ihole][ims1] 
	               * gstholevecs[igstdeg][jhole][ims2] ;
                }
                if ( eout < vmtzero - 0.01 || eout  > vmtzero + 0.01 ) {
                  sumre += rerho * dum ;
                  sumim += imrho * dum ;
                }
              }
            }
          }
          ijre = 2 * ( ich + nchannels * jch ) ;        


/*  without cutting */
/*
          tau[ijre  ] -= sumre ;
          tau[ijre+1] -= sumim ;
*/
          deltaEoverT = (eout - efermi)/TEMPERATURE ;
	  if ( deltaEoverT > 10. ) fermifct = 1.e-7  ;
          else if ( deltaEoverT < -10. ) fermifct = 1.-1.e-7  ;
          else fermifct = 1./(1.+exp(deltaEoverT)) ;

          if ( ich == jch ) {
            kappa = sqrt( fabs( eout - vmtzero ) ) ;
/*
CHANGE
            interpolar( neenv, eenv, wnwj[lvalsh], 1, &eout, &wnwj0 ) ;
*/
	    besneu( eout>vmtzero?0:1, kappa*gspot.wsr, lvalsh, &bes, &neu, 
                    &dbes, &dneu );
	    dbes *= kappa ; dneu *= kappa ; 

	    /*
	    rsq1(gspot.a,gspot.b,eout,singlepofi,lvalsh,gspot.nr-1,&nod0,
                 gspot.nr,rofi,&singledpr0,gspot.v,&singlepr0,gspot.z);
	    */
	    rsq1(pot0.a,pot0.b,eout,singlepofi,lvalsh,pot0.nr-1,&nod0,
                 pot0.nr,rofi,&singledpr0,pot0.v,&singlepr0,pot0.z);

/* calculate wnwj = -(kK)^{-1} = W[neu,P/r]/W[bes,P/r] */
	    rWjR = bes * singledpr0 - dbes * singlepr0 - bes * singlepr0 / pot0.wsr ;
	    rWnR = neu * singledpr0 - dneu * singlepr0 - neu * singlepr0 / pot0.wsr ;
	    wnwj0 = rWnR / rWjR ;

            if ( eout >= vmtzero ) {
              t0m1re = -kappa * wnwj0 ;
              t0m1im = -kappa ;
            }
            else {
              t0m1re = -kappa * ( wnwj0 + 1. ) ;
              t0m1im = 0. ;
            }
          }
          else {
            t0m1re = 0.;
            t0m1im = 0. ;
          }
/*
          printf("%lf %lf %lf %d %lf %lf %lf %lf %lf\n",eout,kappa,fermifct,
                          lvalsh,wnwj0,t0m1re,t0m1im,tau[ijre],tau[ijre+1] ) ;
*/
          tau[ijre  ] += ( t0m1re * fermifct - sumre ) / ( 1. - fermifct ) ;
          tau[ijre+1] += ( t0m1im * fermifct - sumim ) / ( 1. - fermifct ) ;
/* */
        }
      }	    
    }

    /* invert tau */
    invdcmplxmat( nchannels, tau ) ;

/*
    for ( ich = 0 ; ich < nchannels ; ich++ ) 
      for ( jch = 0 ; jch < nchannels ; jch++ ) {
	ijre = 2 * ( ich  +  nchannels * jch ) ;
	if ( fabs(tau[ijre]) > 1.e-6 || fabs(tau[ijre+1]) > 1.e-6 ) 
	  printf("tau[%3d,%3d] =  %10lf + %10lf i\n", 
		 ich, jch, tau[ijre], tau[ijre+1] ) ;
      }
*/ 
    imtautrace = 0. ;
    for ( ich = 0 ; ich < nchannels ; ich++ ) {
      ijre = 2 * ( ich  +  nchannels * ich ) ;
      imtautrace +=  tau[ijre+1] ;
    }    
 
/*      calculation of dipole matrix elements    */
/* 
  first,calculate dipole matrix elements for channel functions 
 < g_ig | r_q | phi_ich > = \sum_k  < g_ig | r_q | psi_k > < psi_k || phi_ich >
                     =  \sum_k \sum_I <g_ig|r_q|f_I> c_{Ik} Wdagger_{k,ich} 
        Now we have,       c_{Ik} ==  revec[ I + nstates * jk[k] ]
	and  < g_ig | r_q | f_I > == dipmatfresbas[ iq ] [ igstdeg ] [ I ]^*
   ATTENTION ! the following code is only valid if dipmatfresbas is real !!
*/
      dipmatsqsum = 0. ;
      for ( iq = 0 ; iq < 3 ; iq++ ) {
	for ( ich = 0 ; ich < nchannels ; ich++ ) {
	  sum = 0. ;
	  for ( k = 0 ; k < nchannels ; k++ ) {
	    j = jk[k] ;
	    sumi = 0. ;
	    for ( i = 0 ; i < nstates ; i++ ) 
              sumi +=  dipmatfresbas[iq][igstdeg][i] * revec[i + nstates * j] ;
	    sum += sumi * wdagger[k + nchannels * ich] ;
	  }
	  dipmatchabas[iq][igstdeg][ich] = sum ;
	  dipmatsqsum += sum*sum ;
	}
      }
/* 
   second, we multiply   < g_ig | r_q | phi_ich >  with  Z = -J X + N 
   (note that  Z= -J X + N  both for outgoing/incoming BC) 
   and thus get the correctly normalized dipole matrix elements 
*/
      for ( iq = 0 ; iq < 3 ; iq++ ) {
	for ( jch = 0 ; jch < nchannels ; jch++ ) {
	  sum = dipmatchabas[iq][igstdeg][jch] * rkneu[jch] ;/* N,diagonal*/
	  for ( ich = 0 ; ich < nchannels ; ich++ ) 
	    sum -= dipmatchabas[iq][igstdeg][ich] * rkbes[ich] 
		 * mat1[ ich + nchannels * jch ] ;
	  dipmatele[iq][igstdeg][jch] = sum ;
/*          printf("iq=%d igstdeg=%d jch=%d dipmatele=%lf\n",iq,igstdeg,jch,
                    dipmatele[iq][igstdeg][jch] ) ;
*/
	}
      }
      
      for ( iq = 0 ; iq < 3 ; iq++ ) {
        for (iqp =  0 ; iqp < 3 ; iqp++ ) {
	  sum = 0. ;
	  for ( ich = 0 ; ich < nchannels ; ich++ ) {
	    for ( jch = 0 ; jch < nchannels ; jch++ ) {
	      ijre = 2 * ( ich + nchannels * jch ) ;
	      imtauij = tau[ijre+1] ;
	      sum  += dipmatele[iq][igstdeg][ich] 
		    * imtauij * dipmatele[iqp][igstdeg][jch] ;
	    }
	  }
	  sigmaqqgst[iq][iqp][igstdeg] = sum ;
	}
      }

    } /* end of igstdeg loop */ 


    fprintf(fpq,"%11.5lf", RYDINEV * omega ) ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      for (iqp =  0 ; iqp < 3 ; iqp++ ) {
        sum = 0. ;
        for ( igstdeg = 0 ; igstdeg < ngstdeg ; igstdeg++ ) 
          sum += sigmaqqgst[iq][iqp][igstdeg] ;
	fprintf(fpq," %11.7lf", sum ) ;
	sqq[iq][iqp] = sum ;
      }
    }
    fprintf(fpq,"\n");

    fprintf(fp,"%11.5lf", RYDINEV * omega ) ;
    fprintf(fp," %11.7lf", sqq[0][0] + sqq[1][1]+ sqq[2][2] );
    fprintf(fp," %11.7lf", 
	    0.5*( sqq[0][0] + sqq[2][2] - sqq[0][2] - sqq[2][0] ) ) ;
    fprintf(fp," %11.7lf", 
	    0.5*( sqq[0][0] + sqq[2][2] + sqq[0][2] + sqq[2][0] ) ) ;
    fprintf(fp," %11.7lf", sqq[1][1] ) ;
    fprintf(fp," %11.7lf", sqq[0][0] ) ;
    fprintf(fp," %11.7lf", sqq[2][2] ) ;
    fprintf(fp," %11.5lf\n", RYDINEV * ( omega + ecore ) ) ;
 

    if ( ngstdeg != 1 )
      printf("ngstdeg != 1. No d0 system !! No simple lm decomposition.\n") ;
    else if ( nvalshholes != valshdeg ) 
      printf("nvalshholes != valshdeg. No d0 system !! No simple lm decomposition.\n") ;
    else {
/*
      fprintf(fpm,"%8.4lf", RYDINEV * ( omega + ecore ) ) ;
      for ( m1 = -lvalsh ; m1 <= lvalsh ; m1++ ) {
        for ( m2 = -lvalsh ; m2 <= lvalsh ; m2++ ) { 
          sum = 0. ;
          for ( ich0 = 0 ; ich0 < ncha0basis ; ich0++ ) {
            for ( jch0 = 0 ; jch0 < ncha0basis ; jch0++ ) {
              for ( ispin = 0 ; ispin < 2 ; ispin++ ) {
                ihole = 2 * ( lvalsh + m1 ) + ispin ;
                ich = ncha0basis * ihole + ich0 ;
                for ( jspin = 0 ; jspin < 2 ; jspin++ ) {
                  jhole = 2 * ( lvalsh + m2 ) + jspin ;
                  jch = ncha0basis * jhole + jch0 ;
                  dipmatsqsum = 0. ;
                  for ( iq = 0 ; iq < 3 ; iq++ ) 
                    dipmatsqsum += dipmatele[iq][0][ich] * dipmatele[iq][0][jch] ;
                  ijre = 2 * ( ich + nchannels * jch ) ;
	          sum += dipmatsqsum * tau[ijre+1] ;
		}
	      }
	    }
	  }
          fprintf(fpm," %7.5lf ", sum ) ;
	}
      }
      fprintf(fpm,"\n") ;
*/
      for ( m1 = 0 ; m1 <= 2*lvalsh ; m1++ ) {
        for ( m2 = 0 ; m2 <= 2*lvalsh ; m2++ ) { 
          for ( n1 = 0 ; n1 <= 2*lvalsh ; n1++ ) {
            for ( n2 = 0 ; n2 <= 2*lvalsh ; n2++ ) { 
              sum = 0 ;
              for ( ich0 = 0 ; ich0 < ncha0basis ; ich0++ ) {
                for ( jch0 = 0 ; jch0 < ncha0basis ; jch0++ ) {
                  for ( ispin = 0 ; ispin < 2 ; ispin++ ) {
                    for ( jspin = 0 ; jspin < 2 ; jspin++ ) {
                      ich = ncha0basis * ( 2*m1 + ispin ) + ich0 ;
                      jch = ncha0basis * ( 2*m2 + jspin ) + jch0 ;
                      dipmatsqsum = 0. ;
                      for ( iq = 0 ; iq < 3 ; iq++ ) 
                        dipmatsqsum += dipmatele[iq][0][ich] * dipmatele[iq][0][jch] ;
                      ich = ncha0basis * ( 2*n1 + ispin ) + ich0 ;
                      jch = ncha0basis * ( 2*n2 + jspin ) + jch0 ;
                      ijre = 2 * ( ich + nchannels * jch ) ;
	              sum += dipmatsqsum * tau[ijre+1] ;
		    }
		  }
		}
	      }
                dtaud[m1][n1][n2][m2] = sum ; 
	    }
	  }
	}
      }
      fprintf(fpm,"%11.5lf", RYDINEV * omega ) ;
      for ( k1 = 0 ; k1 <= 2*lvalsh ; k1++ ) {
        for ( k2 = 0 ; k2 <= 2*lvalsh ; k2++ ) { 
	  sum = 0. ;
	  for ( m1 = 0 ; m1 <= 2*lvalsh ; m1++ ) 
	    for ( m2 = 0 ; m2 <= 2*lvalsh ; m2++ ) 
	      for ( n1 = 0 ; n1 <= 2*lvalsh ; n1++ )
		for ( n2 = 0 ; n2 <= 2*lvalsh ; n2++ ) 
                  sum += u[m1][k1] * u[k1][n1] * u[n2][k2] * u[k2][m2]
		    * dtaud[m1][n1][n2][m2] ;
          fprintf(fpm," %8.5lf ", sum ) ;
	}
      }
      fprintf(fpm,"\n") ;
    }

    } /* end of omega loop */
      

    free( lambda ) ;
    free( qmetric ) ;

    fclose(fpm) ;
    fclose(fpq) ;
    fclose(fp) ;




 /* clean up at very end */
    free(tau) ; free(eenv) ;
    free(kga); free(rkbes); free(drkbes); free(rkneu); free(drkneu);
    for ( i=0; i<nrwf; i++) {
      free(ppofi[i]);
    }
    free(ppofi); 
    free(e0p); free(pr0); free(dpr0); 
    free(rofi); free(vofi); free( gspot.v ) ; free(singlepofi) ;
    printf("THE END\n") ;


// for MPI
    }
    else{
      rho = calloc4double( ntas, ntas, 2, neenv ) ;
      chi = calloc4double( ntas, ntas, 2, neenv ) ;
      wnwj = calloc2double( LMAX_cont + 1, neenv ) ;
    }
    MPI_Finalize();


}
