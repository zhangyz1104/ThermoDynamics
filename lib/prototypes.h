/* Instructions.c */
void init(struct All_variables *E);
/* Geodynamics.c */
double Tm(struct All_variables *E);
void solidus_liquidus(struct All_variables *E);
void Equilibrium(struct All_variables *E);
void binary_solidus_liquidus(struct All_variables *E);
void single_solidus_liquidus(struct All_variables *E);
double Newton_iteration_phi(struct All_variables *E, double T, int i);
double K(struct All_variables *E, double T, int i, int n);
/* Crystallization.c */
void Frac_crystal(struct All_variables *E);