#include "equation.h"
#include "cartesian.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct equation_ftable_t {
  void (*step)(equation_t *eq);
  void (*setup_spectral)(equation_t *eq);
  void (*cleanup)(equation_t *eq);
} equation_ftable_t;

static void cahn_hilliard_step_semi_implicit(equation_t *eq);
static void cahn_hilliard_setup_spectral(equation_t *eq);
static void cahn_hilliard_cleanup(equation_t *eq);

static const equation_ftable_t cahn_hilliard_ftable =
  {
   .step = cahn_hilliard_step_semi_implicit,
   .setup_spectral = cahn_hilliard_setup_spectral,
   .cleanup = cahn_hilliard_cleanup
  };

static f64 default_nonlinear(f64 *params) {
  f64 c = params[0];
  f64 A = params[1];
  return 2.0*A*c * (1.0-c) * (1.0-2.0*c);
}

static f64 default_mu(f64 *params) {
  return default_nonlinear(params);
}


static void create_fftw_plans(equation_t *eq) {
  grid_t *grid = eq->grid;
  // need to cast this to that data type for fftw-plans
  // need to think about how to do this differently
  cahn_hilliard_data_t *data = eq->data;
  uint64_t dim = grid->dim;
  uint64_t *N = grid->N;

  int32_t fftw_N[dim];
  for (uint64_t i = 0; i < dim; i++) fftw_N[i] = N[dim-1-i];

  bc_type_t bc_type = grid->bc.type;

  fftw_r2r_kind *kind_fwd = alloc(dim, sizeof(fftw_r2r_kind));
  fftw_r2r_kind *kind_bwd = alloc(dim, sizeof(fftw_r2r_kind));

  for (uint64_t d = 0; d < dim; d++) {
    kind_fwd[d] = FFTW_REDFT10;
    kind_bwd[d] = FFTW_REDFT01;
  }

  switch (bc_type) {
  case BC_NEUMANN:
    data->fwd_plan = fftw_plan_many_r2r(dim, fftw_N, 1,
					data->wrk1,
					NULL, 1, 1,
					data->wrk2,
					NULL, 1, 1, kind_fwd,
					FFTW_MEASURE);
    data->bwd_plan = fftw_plan_many_r2r(dim, fftw_N, 1,
					data->wrk1,
					NULL, 1, 1,
					data->wrk2,
					NULL, 1, 1, kind_bwd,
					FFTW_MEASURE);
    break;

  default:
    fprintf(stderr, "\e[1:31m Error\e[0m"
	    ": Unknown boundary condition type\n");
    exit(1);
    break;
  }
}

equation_t *create_cahn_hilliard(uint64_t dim, uint64_t *N,
				 f64 *L, f64 M, f64 kappa, f64 A,
				 bc_type_t bc) {

  // Allocate equation structure
  equation_t *eq = alloc(1, sizeof(equation_t));

  // set type and function table
  eq->type = EQUATION_CAHN_HILLIARD;
  eq->ftable = &cahn_hilliard_ftable;

  // Create grid
  eq->grid = create_grid(dim, N, L);
  if (BC_PERIODIC == bc)
    set_periodic_boundary_condition(eq->grid);
  if (BC_NEUMANN == bc)
    set_neumann_boundary_condition(eq->grid);
  if (BC_DIRICHLET == bc)
    set_dirichlet_boundary_condition(eq->grid, 0.0);
    

  // Initialize time stepping
  eq->dt = 0.01;
  eq->max_iter = 1000;
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

  if (eq->ftable->setup_spectral) {
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
  cahn_hilliard_data_t *data = (cahn_hilliard_data_t *)eq->data;

  if (data) {
    safe_free( data->c );
    safe_free( data->wrk1 );
    safe_free( data->wrk2 );

    safe_free( data->linear_op );

    if (data->fwd_plan) fftw_destroy_plan( data->fwd_plan );
    if (data->bwd_plan) fftw_destroy_plan( data->bwd_plan );

    safe_free( *data );
  }

  safe_free( *eq );
}

/* ------------------------------------------------------------------ */
void set_nonlinear_cahn_hilliard(equation_t *eq,
				 f64 (*nonlinear)(f64 *params)) {
  if (eq->type != EQUATION_CAHN_HILLIARD) {
    fprintf(stderr, "\e[1:31m Error\e[0m"
            ": Equation is not Cahn-Hilliard type\n");
    return;
  }
  
  cahn_hilliard_data_t *data = (cahn_hilliard_data_t*)eq->data;
  data->nonlinear = nonlinear;
}

void set_mu_cahn_hilliard(equation_t *eq, f64 (*mu)(f64 *params)) {
  if (eq->type != EQUATION_CAHN_HILLIARD) {
    fprintf(stderr, "\e[1:31m Error\e[0m"
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
    


