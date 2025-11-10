#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "types.h"
#include "memory_management.h"

// forward declare
typedef struct grid_t grid_t;

typedef enum {
	      BC_PERIODIC,
	      BC_NEUMANN,
	      BC_DIRICHLET
} boundary_type_t;    

typedef struct bc_t {
  boundary_type_t type;
  f64 value; // for Dirichlet BC
} bc_t;


typedef struct grid_t {

  uint64_t dim;
  uint64_t *N;
  uint64_t size;
  f64 *L;
  f64 *restrict k2;
  f64 *dx;

  // boundary conditions - can be extended to per-axis
  boundary_condition_t bc;
} grid_t;

grid_t create_grid(uint64_t dim, uint64_t *N, f64 *L);
void grid_destroy_interal(grid_t **grid);
#define grid_destroy(grid) grid_destroy_internal((grid_t **) &(grid))

void set_periodic_boundary_condition(grid_t *grid);
void set_neumann_boundary_condition(grid_t *grid);
void set_dirichlet_boundary_condition(grid_t *grid, f64 val);

#endif // CARTESIAN_H
