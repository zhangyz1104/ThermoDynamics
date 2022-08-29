#include <stdlib.h>
#include "global_defs.h"

/* fractional crystallization model */
void Frac_crystal(struct All_variables *E) {
  int continue_cooling, its;
  double T, phi, Kn, cl, cs;
  int nT, n, i;
  FILE *fp;

  const double dT = 1;
  const int nc = E->comp.ncomp;
  E->grid.nr = 1;
  E->grid.P[1] = 0;
  Tm(E);

  // initialization
  solidus_liquidus(E);
  T = E->comp.Tliq[1];
  nT = (int) (T / dT);

  double rc[nT], cumulate_comp[nc + 1][nT], melt_comp[nc + 1][nT];
  double Tliq[nT], Tsol[nT];

  rc[0] = 0;
  its = 1;
  continue_cooling = 1;

  // cooling iteration
  while (continue_cooling) {

    // compute equilibrium state for a given T
    phi = max(0, min(1, Newton_iteration_phi(E, T, 1)));

    for (n = 1; n <= nc; n++) {
      Kn = K(E, T, 1, n);
      cl = E->comp.c0[n] / (phi + (1 - phi) * Kn);
      cl = max(0, min(1, cl));
      cs = E->comp.c0[n] / (phi / Kn + (1 - phi));
      cs = max(0, min(1, cs));

      // update segregated melt and solid
      E->comp.c0[n] = melt_comp[n][its] = cl;
      cumulate_comp[n][its] = cs;
    }
    rc[its] = rc[its - 1] + (1 - phi) * (1 - rc[its - 1]);
    Tliq[its] = E->comp.Tliq[1];
    Tsol[its] = E->comp.Tsol[1];

    T -= dT;
    solidus_liquidus(E);
    if (T < E->comp.Tsol[1])
      continue_cooling = 0;
    else
      its += 1;
  }

  // norm cumulate radius
  fprintf(stderr, "  (error = %.3e)", rc[its] - 1);
  for (i = 1; i <= its; i++) {
    rc[i] /= rc[its];
  }

  // output
  fp = fopen("../data/fractional_crystal", "w");
  for (i = 1; i <= its; i++) {
    fprintf(fp, "%.3f %.3f %.3f ", T + (its - i + 1) * dT, Tliq[i], Tsol[i]);
    for (n = 1; n <= nc; n++) {
      fprintf(fp, "%.3f %.3f ", cumulate_comp[n][i], melt_comp[n][i]);
    }
    fprintf(fp, "%.3f\n", rc[i]);
  }
}