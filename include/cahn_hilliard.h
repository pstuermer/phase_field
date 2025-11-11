#ifndef CAHN_HILLIARD_H
#define CAHN_HILLIARD_H

#include "types.h"
#include "memory_management.h"
#include <fftw3.h>


typedef struct cahn_hilliard_data_t {
  f64 M;
  f64 kappa;
  f64 A;

  f64 *restrict c;
  f64 (*nonlinear)(f64 *params); // df/dc
  f64 (*mu)(f64 *params);

  f64 *restrict wrk1;
  f64 *restrict wrk2;

  fftw_plan fwd_plan;
  fftw_plan bwd_plan;

  f64 *restrict linear_op;

} cahn_hilliard_data_t;



equation_t *create_cahn_hilliard(uint64_t dim, uint64_t *N,
				 f64 *L, f64 M, f64 kappa, f64 A,
				 bc_type_t bc);

// set function pointers for nonlinear term and chemical potential
void set_nonlinear_cahn_hilliard(equation_t *eq,
				 f64 (*nonlinear)(f64 *params));
void set_mu_cahn_hilliard(equation_t *eq, f64 (*mu)(f64 *params));

// accesor functions for interface
f64 *cahn_hilliard_get_field(equation_t *eq);
void cahn_hilliard_set_field_value(equation_t *eq, uint64_t i, f64 value);
f64 cahn_hilliard_get_field_value(equation_t *eq, uint64_t i);

#endif // CAHN_HILLIARD_H
