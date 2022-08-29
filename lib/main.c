#include <stdlib.h>
#include "global_defs.h"

int main() {
  struct All_variables *E = (struct All_variables *) malloc(sizeof(struct All_variables));

  init(E);
  // Frac_crystal(E);
  solidus_liquidus(E);
  Equilibrium(E);

  return 1;
}
