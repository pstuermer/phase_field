#ifndef EQUATION_H
#define EQUATION_H

#include "types.h"
#include "memory_management.h"
#include <fftw3.h>

// forward declare
struct grid_t;
struct equation_ftable_t;
//enum bc_type_t;

typedef enum {
	      EQUATION_NONE,
	      EQUATION_CAHN_HILLIARD
} equation_type_t;

typedef enum {
	      TIME_NONE,
	      TIME_SEMI_IMPLICIT
} time_integration_t;

typedef struct equation_t {

  equation_type_t type; // what type of equation
  const struct equation_ftable_t *ftable; // function-pointer table
  void *data; // equation-specific data

  // assigned grid
  struct grid_t *grid;

  // time-stepping
  f64 dt;
  uint64_t max_iter;
  uint64_t iter;
  time_integration_t time_method;

  
} equation_t;

typedef struct equation_ftable_t {
  void (*step)(equation_t *eq);
  void (*setup_spectral)(equation_t *eq);
  f64* (*get_field)(equation_t *eq);
  void (*set_field)(equation_t *eq, uint64_t i, f64 val);
  f64 (*get_free_energy)(equation_t *eq);
  void (*cleanup)(equation_t *eq);
} equation_ftable_t;

void equation_destroy_interal(equation_t **eq);
#define equation_destroy(eq) equation_destroy_internal((equation_t **) &(eq))


// boilerplate code can be extended to other types via X-macros
void set_iter(equation_t *eq, f64 dt,
	      uint64_t max_iter, uint64_t iter);


// Initialize concentration fields
void equation_init_uniform(equation_t *eq, f64 value);
void equation_init_random(equation_t *eq, f64 amplitude);
void equation_init_custom(equation_t *eq, f64 (*init_funct)(f64 *coords));


// time propagation
void set_semi_implicit_prop(equation_t *eq);
void time_prop(equation_t *eq);
void equation_run(equation_t *eq);

// conservation checks
f64 equation_compute_mass(equation_t *eq);
f64 equation_compute_free_energy(equation_t *eq);

// I/O
void equation_output_csv(equation_t *eq, const char *filename, uint8_t append);
void equation_output_field(equation_t *eq, const char *filename);

#endif // EQUATION_H
