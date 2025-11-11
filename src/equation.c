#include "../include/equation.h"
#include "../include/cartesian.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
void equation_destroy_internal(equation_t **eq) {
  if (NULL == eq || NULL == *eq) return;

  if ((*eq)->ftable && (*eq)->ftable->cleanup) {
    (*eq)->ftable->cleanup(*eq);
  }

  grid_destroy((*eq)->grid);
  safe_free( *eq );
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

void equation_init_custom(equation_t *eq, f64(*init_func)(f64 *coords)) {
  if (!eq->ftable || !eq->ftable->set_field) {
    fprintf(stderr, "\e[1:31m ERROR\e[0m: No set field function.\n");
    exit(1);
  }
  
  grid_t *grid = eq->grid;
  f64 coords[3] = {0.0, 0.0, 0.0};

  for (uint64_t i = 0; i < grid->size; i++) {
    // Convert linear index to coordinates
    coords[0] = i%(grid->N[0]);
    if (grid->dim > 1) coords[1] = (i/(grid->N[0]))%(grid->N[1]);
    if (grid->dim > 2) coords[2] = (i/((grid->N[0])*(grid->N[1])));
    
    eq->ftable->set_field(eq, i, init_func(coords));
  }
}



/* ------------------------------------------------------------------ */
void time_prop(equation_t *eq) {
  if (!eq->ftable || !eq->ftable->step) {
    fprintf(stderr, "\e[1:31m ERROR\e[0m: No time stepping function.\n");
    exit(1);
  }

  eq->ftable->step(eq);
}

void equation_run(equation_t *eq)
{
  printf("Starting simulation: ");
  if (EQUATION_CAHN_HILLIARD == eq->type) printf("Cahn-Hilliard\n");
  printf("dt = %.2e, max_itere = %lu\n", eq->dt, eq->max_iter);

  for (eq->iter = 0; eq->iter < eq->max_iter; eq->iter++) {
    time_prop(eq);

    if (!(eq->iter % 100)) {
      printf("Iteration %lu/%lu\n", eq->iter, eq->max_iter);
    }
  }

  printf("Simulation complete.\n");
  // could put some output values here
}

/* ------------------------------------------------------------------ */
// mass conservation check

f64 equation_compute_mass(equation_t *eq) {
  if (!eq->ftable || !eq->ftable->get_field) {
    fprintf(stderr, "\e[1:31m ERROR\e[0m: No get field function.\n");
    exit(1);
  }
  
  register f64 val = 0.0;
  f64 *field = eq->ftable->get_field(eq);
  for (uint64_t i = 0; i < eq->grid->size; i++)
    val += field[i];

  return val;
  
}
