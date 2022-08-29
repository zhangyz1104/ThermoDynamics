#include "global_defs.h"
#include <stdlib.h>
#include <stdio.h>

static void set_parameters(struct All_variables *E) {
  const double dS = 320;
  // set composition parameters
  E->comp.ncomp = 3;
  E->comp.c0[1] = 0.5;
  E->comp.c0[2] = 0.35;
  E->comp.c0[3] = 0.15;
  E->comp.T0[1] = 2050;
  E->comp.T0[2] = 1500;
  E->comp.T0[3] = 1350;
  E->comp.A[1] = 60;
  E->comp.A[2] = 100;
  E->comp.A[3] = 120;
  E->comp.B[1] = 0;
  E->comp.B[2] = 0;
  E->comp.B[3] = 0;
  E->comp.r[1] = 50;
  E->comp.r[2] = 20;
  E->comp.r[3] = 10;
  E->comp.L[1] = dS * E->comp.T0[1];
  E->comp.L[2] = dS * E->comp.T0[2];
  E->comp.L[3] = dS * E->comp.T0[3];

  // equilibrium state
  E->grid.nr = 141;
  E->grid.r0 = 340;
  E->grid.r1 = 1740;
  E->grid.nT = 100;

  // crystallization model
//  E->grid.nr = 2;
//  E->grid.r0 = 0;
//  E->grid.r1 = 0;
}

static void init_field(struct All_variables *E) {
  int i;

  const int nno = E->grid.nr;
  const int ncomp = E->comp.ncomp;
  const double k = 2 * 3.14 * 6.67e-11 / 3 * 3400 * 3400 * 1e6 * 1e-9;

  for (i = 1; i <= ncomp; i++) {
    E->comp.K[i] = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.cl[i] = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.cs[i] = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.Tm[i] = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.phi = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.Tsol = (double *) malloc((nno + 1) * sizeof(double));
    E->comp.Tliq = (double *) malloc((nno + 1) * sizeof(double));

    E->grid.r = (double *) malloc((nno + 1) * sizeof(double));
    E->grid.P = (double *) malloc((nno + 1) * sizeof(double));
  }

  // initialize grid configuration
  for (i = 1; i <= nno; i++) {
    E->grid.r[i] = E->grid.r0 + (E->grid.r1 - E->grid.r0) / (E->grid.nr - 1) * (i - 1);
    E->grid.P[i] = k * (E->grid.r1 * E->grid.r1 - E->grid.r[i] * E->grid.r[i]);
  }

  // initialize Tm(P)
  Tm(E);
}

void init(struct All_variables *E) {
  set_parameters(E);
  init_field(E);
}