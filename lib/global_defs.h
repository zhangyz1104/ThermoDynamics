#include <stdio.h>

#define MAX_NCOMP 13
#define max(A, B) (((A) > (B)) ? (A) : (B))
#define min(A, B) (((A) < (B)) ? (A) : (B))

struct GRID {
  int nr;
  int nT;
  double r0;
  double r1;

  double *r;
  double *P;
};

struct COMP {
  int ncomp;
  double T0[MAX_NCOMP];
  double A[MAX_NCOMP];
  double B[MAX_NCOMP];
  double L[MAX_NCOMP];
  double r[MAX_NCOMP];

  double c0[MAX_NCOMP];

  double *Tm[MAX_NCOMP];
  double *K[MAX_NCOMP];
  double *cl[MAX_NCOMP];
  double *cs[MAX_NCOMP];
  double *Tsol;
  double *Tliq;
  double *phi;

  FILE *solliq_file;
  FILE *equilibrium_file;
};

struct All_variables {
  struct GRID grid;
  struct COMP comp;
};

#include "prototypes.h"