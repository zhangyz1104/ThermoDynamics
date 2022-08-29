#include <stdlib.h>
#include "global_defs.h"

int main() {
  struct All_variables *E = (struct All_variables *) malloc(sizeof(struct All_variables));

  // set parameters and allocate memory
  init(E);

  switch (E->mode) {
    // equilibrium state {phi, c^{s, l}}(P, T, c0)
  case 1:fprintf(stderr, "Equilibrium state\n");
    Equilibrium(E);
    break;

    // fractional crystallization model
  case 2: fprintf(stderr, "Fractional crystallization\n");
    Frac_crystal(E);
    break;

    // solidus and liquidus for a binary system
  case 3: fprintf(stderr, "Binary system\n");
    binary_solidus_liquidus(E);
    break;

  default:fprintf(stderr, "Invalid value of this method\n");
    return 0;
  }

  return 1;
}
