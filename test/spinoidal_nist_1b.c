// or something like that
#include "../include/equation.h"
#include "../include/cahn_hillard.h"

f64 chem_pot(f64 *params) {
  return param[0]*(param[1]-param[2])*(param[1]-param[2])
    *(param[3]-param[1])*(param[3]-param[1]);
}

f64 custom_init(f64 *coords) {
  f64 c0 = 0.5;
  f64 eps = 0.01;
  f64 x = coords[0];
  f64 y = coords[1];
  return c0 + eps*(cos(0.105*x)*cos(0.11*y) +
		   (cos(0.13*x)*cos(0.087*y) *
		    cos(0.13*x)*cos(0.087*y)) +
		   cos(0.025*x-0.15*y) *
		   cos(0.07*x-0.02*y));
}

int main() {

  uint64_t N[2] = {128, 128}; // gridpoints
  f64 L[2] = {200.0, 200.0};

  f64 rho = 5.0;
  f64 kappa = 2.0;
  f64 M = 5.0;

  equation_t *eq = create_cahn_hilliard(2, N, L, M, kappa,
					rho, BC_NEUMANN);

  // set time step and duration
  f64 dt = 0.1;
  uint64_t max_iter = 10; // run to 1
  set_iter(eq, dt, max_iter, 0);

  // set custom_init

  // loop over time intervals for output
  equation_run(eq);
  // write output

  equation_destroy(eq);

  return 0;
}

  
