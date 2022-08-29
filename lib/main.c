#include <stdlib.h>
#include "global_defs.h"

int main() {
  struct All_variables *E = (struct All_variables *) malloc(sizeof(struct All_variables));

  init(E);
  // fractional crystallization model
  // Frac_crystal(E);

  // equilibrium state {phi, c^{s, l}}(P, T, c0)
  // solidus_liquidus(E);
  // Equilibrium(E);

  // solidus and liquidus for a binary system
  binary_solidus_liquidus(E);

  return 1;
}
