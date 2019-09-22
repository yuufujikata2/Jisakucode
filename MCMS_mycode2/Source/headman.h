#define iabs(i) ( (i) > 0 ? (i) : -(i) )
#define lmorder( l, m ) ( (l)*(l) + ( (m)>0 ? 2*(m)-1 : -2*(m) ) )
#include <mpi.h>

struct Lmtopot {
  int z; double wsr; double a; double b; int nr; int nspin; double *v;
};

struct Lmtopot readlmtopot(char *filename) ;

void readinput(char* absorbergsfile, char* absorberchfile,
               int* nshell, char* corewffile, double* eenvmin, double* eenvmax,
               double* deenv, int* lmax, int* lsym, double* ksicor, 
	       double *ksival, double *tendqval, double *hz, double* vmtzero, 
	       double* imvmtzero, double* efermi, double* omegamin, 
               double* omegamax, double* domega, 
               double* vmixpar, double* rkscale,
               double* e0min, double* e0max, int* nrwf1,
	       int* lvalsh, int* nvalel, int* withrefl ,int myid) ;

double* readcofi(char* cowffilename, int iz, int nr, 
		 double a, double wsr, int* ncorsh, int* lcorsh, 
		 double* ecore) ;
int potcompatible(struct Lmtopot p1, struct Lmtopot p2) ;

void calcvhart(int nr, double* rofi, double* pofi, double* vhart) ;

void makestates( int nshells, int *lsh, int *sorb1sh, int nconfs, 
		 int nelectrons, int **occ, int nstates0, 
		 struct Fock **pstate) ;

void environment(int nr, double *cofi, struct Lmtopot pot0, 
		 int nshell, double vmtzero, double imvmtzero, 
		 int nep, double *energy, int lsym, 
		 int ntas, double ****rho, double ****chi, double **wnwj,
                 int numproces, int myid, MPI_Comm comm ) ;
     

/* double* makerofi(int nr, double wsr, double* b, double a); */
int calcppofi(int nr, double a, double b, double *rofi, double *vofi, 
	      int iz, int l, double emin, double emax, double valin, 
	      double sloin, int *nod0, double **ep, double **val,  
	      double **slo, double ***ppofi) ; 
int gramschmidtplus(int nr, double *rofi, int nrwf, double ***pppofi,
		    double ***psmat ) ;
double overlap( int nr, double *rofi, double *pofi1, double *pofi2 ) ;

int invuppertri(int dim, double **mat, double ***pinvmat ) ;

int transfe0pr0dpr0( int nrwf, double **smat, double **tmat, 
		     double *e0p, double *pr0, double *dpr0, 
		     double ***pe0po, double **ppr0o, double **pdpr0o ) ;

void reorderrwfs( double *cofi, double ecore, int *pnrwf, int nopen1, 
		  double ***pppofi, double **pe0p, double **ppr0, 
		  double **pdpr0 ) ;
void makeshells( int nshells, int lcorsh, int lvalsh, double ksicor, 
		 double ksival, double tendqval, int **plsh, int **psorb1sh, 
                 double **pksish, double **ptendqsh) ;
void makeconfs( int nshells, int *lsh, int nvalel, int which,
		int *pnconfs, int ***pocc, int *pnstates ) ;
void besneu(int modified, double x, int l, double* pjl, double* pnl, double* pjlp, 
	    double* pnlp) ;
int solve_genev(int n, double* a, double* b,
                double* alpha, double* beta, double* revec) ;
double  defintegr(double *f,  double *x, int Nx ) ;
void interpolar(int nin, double *xin, double *yin, 
		int nout, double *xout, double *yout) ;
int inv_remat( int matdim, double* mat ) ;
void rsq1(double a, double b, double e, double* g, int l, int nc,
	  int* nod, int nr, double* rofi, double* slonc, double* v, 
	  double* valnc, int iz) ;
int invdcmplxmat( int matdim , double* mat ) ;



