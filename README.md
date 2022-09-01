# ThermoDynamics
thermodynamic calculations based on R_DMC

## Parameter setting and Outputs
set mode, composition and grid parameters in `Instructions.c`  

`mode==1`: equilibrium states for different pressure and temperature -> `equilibrium_state` `solidus_liquidus`  

`mode==2`: fractional crystallization model -> `fractional_crystal`  

`mode==3`: solidus and liquidus for a binary system -> `binary_sol_liq`   

`mode==4`: debug solidi and liquidi for a single system -> `olv_sol_liq` `pxn_sol_liq` `fsp_sol_liq`

Default: olv-pxn-fsp whole-mantle LMO (lunar magma ocean) conditions
