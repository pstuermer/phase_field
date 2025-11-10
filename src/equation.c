#include "equation.h"
#include "cartesian.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void create_fftw_plans(equation_t *eq) {
  grid_t *grid = eq->grid;
  uint64_t dim = grid->dim;
  uint64_t *N = grid->N;

  uint64_t fftw_N[dim];
  for (uint64_t i = 0; i < dim; i++) fftw_N[i] = N[dim-1-i];

  boundary_type_t bc_type = grid->bc.type;

  switch (bc_type) {
  case BC_NEUMANN:
    data->fwd_plan = fftw_plan_many_r2r(dim, fftw_N, 1, eq->data->wrk1,
					NULL, 1, 1, eq->data->wrk2,
					NULL, 1, 1, FFTW_REDFT10,
					FFTW_MEASURE);
    data->bwd_plan = fftw_plan_many_r2r(dim, fftw_N, 1, eq->data->wrk1,
					NULL, 1, 1, eq->data->wrk2,
					NULL, 1, 1, FFTW_REDFT01,
					FFTW_MEASURE);
    break;
  }

  default:
    fprintf(stderr, ANSI_COLOR_RED "ERROR" ANSI_COLOR_RESET
	    ": Unknown boundary condition type\n");
    exit(1);
}


equation_t *create_cahn_hilliard(uint64_t dim, uint64_t *N,
				 f64 *L, f64 M f64 kappa, f64 A,
				 bondary_type_t bc) {

  // Allocate equation structure
  equation_t *eq = alloc(1, sizeof(equation_t));

  // set type and function table
  eq->type = EQUATION_CAHN_HILLIARD;
  eq->ftable = &cahn_hilliard_ftable;

  // Create grid
  eq->grid = create_grid(dim, N, L);
  if (BC_PERIDIC == bc)
    set_periodic_boundary_condition(eq->grid);
  if (BC_NEUMANN == bc)
    set_neumann_boundary_condition(eq->grid);
  if (BC_DIRICHLET == bc)
    set_dirichlet_boundary_condition(eq->grid, 0.0);
    

  // Initialize time stepping
  *(f64*)eq->dt = 0.01;
  *(uint64_t*)&eq->max_iter = 1000;
  eq->iter = 0;
  eq->time_method = TIME_SEMI_IMPLICIT;

  // Allocate Cahn-Hilliard specific data
  cahn_hilliard_data_t *data = alloc(1, sizeof(cahn_hilliard_data_t));

  data->M = M;
  data->kappa = kappa;
  data->A = A;

  // Allocate fields
  data->c = alloc(eq->grid->size, sizeof(f64));
  data->wrk1 = alloc(eq->grid->size, sizeof(f64));
  data->wrk2 = alloc(eq->grid->size, sizeof(f64));

  data->nonlinear = default_nonlinear;
  data->mu = default_mu;

  data->linear_op = alloc(eq->grid->size, sizeof(f64));

  eq->data = data;

  create_fftw_plans(eq);

  if (eq->ftable->setup->spectral) {
    eq->ftable->setup_spectral(eq);
  }

  return eq;
}

/* ------------------------------------------------------------------ */
void equation_destroy_internal(equation_t **eq) {
  if (NULL == eq || NULL == *eq) return;

  if ((*eq)->ftable && (*eq)->ftable->cleanup) {
    (*eq)->ftable->cleanup(*eq);
  }

  grid_destroy((*eq)->grid);
  safe_free( *eq );
}

static void cahn_hilliard_cleanup(equation_t *eq) {
  cahn_hilliard_data_t *data = (cahn_hilliard_data_t*)eq->data;

  if (data) {
    safe_free( data->c );
    safe_free( data->wrk1 );
    safe_free( data->wrk2 );

    safe_free( data->linear_op );

    if (data->fwd_plan) fftw_destroy_plan( fwd_plan );
    if (data->bwd_plan) fftw_destroy_plan( bwd_plan );

    safe_free( data );
  }
}

/* ------------------------------------------------------------------ */
void set_nonlinear_cahn_hilliard(equation_t *eq,
				 f64 (*nonlinear)(f64 *params)) {
  if (eq->type != EQUATION_CAHN_HILLIARD) {
    fprintf(stderr, ANSI_COLOR_RED "Error" ANSI_COLOR_RESET 
            ": Equation is not Cahn-Hilliard type\n");
    return;
  }
  
  cahn_hilliard_data_t *data = (cahn_hilliard_data_t*)eq->data;
  data->nonlinear = nonlinear;
}

void set_mu_cahn_hilliard(equation_t *eq, f64 (*mu)(f64 *params)) {
  if (eq->type != EQUATION_CAHN_HILLIARD) {
    fprintf(stderr, ANSI_COLOR_RED "Error" ANSI_COLOR_RESET 
            ": Equation is not Cahn-Hilliard type\n");
    return;
  }
  
  cahn_hilliard_data_t *data = (cahn_hilliard_data_t*)eq->data;
  data->mu = mu;
}


/* ------------------------------------------------------------------ */
void set_iter(equation_t *eq, f64 dt,
	      uint64_t max_iter, uint64_t iter) {
  *(f64*)&eq->dt = dt;
  *(uint64_t*)&eq->max_iter = max_iter;
  eq->iter = iter;
}

void set_semi_implicit_prop(equation_t *eq) {
  eq->time_method = TIME_SEMI_IMPLICIT;
  
  if (eq->ftable->setup_spectral) {
    eq->ftable->setup_spectral(eq);
  }
}

/* ------------------------------------------------------------------ */
    


