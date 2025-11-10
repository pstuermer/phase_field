#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "types.h"
#include "memory_management.h"

// forward declare
struct equation_t

typedef enum {
	      BC_PERIODIC,
	      BC_NEUMANN,
	      BC_DIRICHLET
} boundary_type_t

typedef struct bc_t {
  boundary_type_t type;
  f64 value; // for Dirichlet BC
} bc_t;


typedef struct grid_t {

  uint64_t dim;
  uint64_t *N;
  uint64_t size;
  f64 *L;
  double *dx;

  // so far global, but can easily be extended to individual
  // along each axis and on each end
  boundary_condition_t bc;

  // an easily extensible, user_defined struct for different
  // kind of systems
  struct equation_t;
} grid_t;

grid_t create_grid(uint64_t dim, uint64_t *N, f64 *L);
void grid_destroy_interal(grid_t **grid);
#define grid_destroy(grid) grid_destroy_internal((grid_t **) &(grid))

void set_periodic_boundary_condition(grid_t *grid);
void set_neumann_boundary_condition(grid_t *grid);
void set_dirichlet_boundary_condition(grid_t *grid, f64 val);

#endif CARTESIAN_H
