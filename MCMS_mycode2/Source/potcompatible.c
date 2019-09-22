#define RELDIFF 1.E-10
#define abs(a) ( a > 0 ? a : -a ) 
#define verynear(a,b) ( abs((a)-(b)) <= RELDIFF * abs((a)+(b)) )

struct Lmtopot {
  int z; double wsr; double a; double b; int nr; int nspin; double *v;
};

int potcompatible(struct Lmtopot p1, struct Lmtopot p2) {
  return (p1.z == p2.z) && verynear(p1.wsr, p2.wsr)
    && verynear(p1.a, p2.a) && verynear(p1.b, p2.b)
    && (p1.nr == p2.nr) && (p1.nspin == p2.nspin) ;
}

#undef verynear
#undef abs
#undef RELDIFF
