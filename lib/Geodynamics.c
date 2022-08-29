#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "global_defs.h"

/* equilibrium distribution coefficients K(T, P) */
double K(struct All_variables *E, double T, int i, int n) {
  double Tm, K;

  Tm = E->comp.Tm[n][i];
  K = exp(E->comp.L[n] / E->comp.r[n] * (1 / T - 1 / Tm));

  return K;
}

/* Newton residual */
static double residual(struct All_variables *E, double T, double phi, int i) {
  int n;
  double r;

  r = 0;
  for (n = 1; n <= E->comp.ncomp; n++)
    r += E->comp.c0[n] / (phi + (1 - phi) * K(E, T, i, n));
  for (n = 1; n <= E->comp.ncomp; n++)
    r -= E->comp.c0[n] / (phi / K(E, T, i, n) + (1 - phi));

  return r;
}

/* Newton iteration for T */
static double Newton_iteration_T(struct All_variables *E, double phi, int i) {
  int n, its;
  double T, r, rp, rm, drdT;
  double Tmin, Tmax;

  const double r_tol = 1e-10; // tolerance for Newton residual
  const double its_tol = 100; // maximum number of iterations
  const double eps_T = 5; // temperature perturbation for finite differencing, degrees

  // set starting guess for Tsol
  T = 0;
  for (n = 1; n <= E->comp.ncomp; n++)
    T += E->comp.c0[n] * E->comp.Tm[n][i];
  Tmax = Tmin = E->comp.Tm[1][i];
  for (n = 1; n <= E->comp.ncomp; n++) {
    if (E->comp.Tm[n][i] > Tmax)
      Tmax = E->comp.Tm[n][i];
    if (E->comp.Tm[n][i] < Tmin)
      Tmin = E->comp.Tm[n][i];
  }

  //  get residual for sum(c0_a/K_a) = 1
  r = residual(E, T, phi, i);
  its = 0;

  // Newton iteration
  while ((fabs(r) > r_tol) && (its < its_tol)) {
    rp = residual(E, T + eps_T, phi, i);
    rm = residual(E, T - eps_T, phi, i);
    drdT = (rp - rm) / eps_T / 2;
    T -= r / drdT;
    r = residual(E, T, phi, i);

    its += 1;
    if (its == its_tol)
      fprintf(stderr, "!!! Newton solver for phi = %.1f has not converged"
                      " after %d iterations when P = %.2f GPa!!!\n", phi, its, E->grid.P[i]);
  }

  return T;
}

/* Newton iteration for phi */
double Newton_iteration_phi(struct All_variables *E, double T, int i) {
  int n, its;
  double phi, r, drdphi;

  const double r_tol = 1e-10; // tolerance for Newton residual
  const double its_tol = 100; // maximum number of iterations

  // set starting guess for Tsol
  phi = (T - E->comp.Tsol[i]) / (E->comp.Tliq[i] - E->comp.Tsol[i]);

  //  get residual for sum(c0_a/K_a) = 1
  r = residual(E, T, phi, i);
  its = 0;

  // Newton iteration
  while ((fabs(r) > r_tol) && (its < its_tol)) {
    drdphi = 0;
    for (n = 1; n <= E->comp.ncomp; n++)
      drdphi -= E->comp.c0[n] * (1 - K(E, T, i, n)) / pow(phi + (1 - phi) * K(E, T, i, n), 2);
    for (n = 1; n <= E->comp.ncomp; n++)
      drdphi += E->comp.c0[n] * (1 / K(E, T, i, n) - 1) / pow(phi / K(E, T, i, n) + (1 - phi), 2);

    phi -= r / drdphi;
    r = residual(E, T, phi, i);

    its += 1;
    if (its == its_tol)
      fprintf(stderr, "!!! Newton solver for equilibrium has not converged"
                      " after %d iterations !!!\n", its);
  }

  return phi;
}

/* ************************************************ */

/* compute melt temperature Tm(P) */
double Tm(struct All_variables *E) {
  int i, n;

  for (i = 1; i <= E->grid.nr; i++)
    for (n = 1; n <= E->comp.ncomp; n++)
      E->comp.Tm[n][i] = E->comp.T0[n] +
          E->comp.A[n] * E->grid.P[i] +
          E->comp.B[n] * E->grid.P[i] * E->grid.P[i];
}

/* compute solidus Tsol(P, C) and liquidus Tliq(P, C) */
void solidus_liquidus(struct All_variables *E) {
  int i;

  for (i = 1; i <= E->grid.nr; i++) {
    E->comp.Tsol[i] = Newton_iteration_T(E, 0, i);
    E->comp.Tliq[i] = Newton_iteration_T(E, 1, i);
    fprintf(E->comp.solliq_file, "%.0f %.3f %.3f\n",
            E->grid.r[i], E->comp.Tsol[i], E->comp.Tliq[i]);
  }
}

/* compute melt fraction phi(T, P, C) and composition c^{s, l}(T, P, C) */
void Equilibrium(struct All_variables *E) {
  int i, j, n;
  double T, phi, Kn;
  double cl, cs;

  for (i = 1; i <= E->grid.nr; i++)
    for (j = 1; j <= E->grid.nT; j++) {
      T = E->comp.Tsol[i] + (E->comp.Tliq[i] - E->comp.Tsol[i]) / (E->grid.nT - 1) * (j - 1);
      phi = max(0, min(1, Newton_iteration_phi(E, T, i)));

      fprintf(E->comp.equilibrium_file, "%.3f ", T);

      for (n = 1; n <= E->comp.ncomp; n++) {
        Kn = K(E, T, i, n);
        cl = E->comp.c0[n] / (phi + (1 - phi) * Kn);
        cl = max(0, min(1, cl));
        cs = E->comp.c0[n] / (phi / Kn + (1 - phi));
        cs = max(0, min(1, cs));

        fprintf(E->comp.equilibrium_file, "%.3f %.3f ", cs, cl);
      }
      fprintf(E->comp.equilibrium_file, "%.3f\n", phi);
    }
}