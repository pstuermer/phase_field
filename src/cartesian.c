#include "cartesian.h"
#include <math.h>
#include <stdio.h>

grid_t *create_grid(uint64_t dim, uint64_t *N, f64 *L) {

  if (dim < 1 || dim >3) {
    fprintf(stderr, "\e[1;31m ERROR \e[0m"
	    ": Grid dimension mst be 1, 2 or 3. Got %lu instead\n",
	    dim);
    exit(1);
  }

  for (uint64_t d = 0; d < dim; d++) {
    if (0 == N[d]) {
      fprintf(stderr, "\e[1;31m ERROR \e[0m"
	      ": Grid size cannot be 0 in dimension %lu\n", d);
      exit(1);
    }

    if (0.0 >= L[d]) {
      fprintf(stderr, "\e[1;31m ERROR \e[0m"
	      ": Grid length must be positive in dimension %lu\n", d);
      exit(1);
    }
  }

  // Allocate grid structure
  grid_t *grid = alloc(1, sizeof(grid_t));

  grid->dim = dim;
  grid->size = 1;
  for (uint64_t d = 0; d < dim; d++) grid->size *= N[d];

  grid->N = alloc(dim, sizeof(uint64_t));
  grid->L = alloc(dim, sizeof(f64));
  grid->dx = alloc(dim, sizeof(f64));

  grid->size = 1;
  for (uint64_t d = 0; d < dim; d++) {
    grid->size *= N[d];
    grid->N[d] = N[d];
    grid->L[d] = L[d];
    grid->dx[d] = L[d]/N[d];
  }

  grid->k2 = alloc(grid->size, sizeof(f64));
  compute_k2(grid);

  // set to periodic boundary conditions as default
  grid->bc.type = BC_PERIODIC;
  grid->bc.value = 0.0;

  return grid;
}

void grid_destroy_internal(grid_t **grid) {

  if (NULL == grid || NULL == *grid) return;

  safe_free( (*grid)->N );
  safe_free( (*grid)->L );
  safe_free( (*grid)->dx );
  safe_free( (*grid)->k2 );

  safe_free( *grid );

}

/* ------------------------------------------------------------------ */

void set_periodic_boundary_conditions(grid_t *grid) {
  if (NULL == grid) {
    fprintf(stderr, "\e[1;31m ERROR\e[0m" 
	    ": NULL grid pointer when setting BC.\n");
    return;
  }

  grid->bc.type = BC_PERIODIC;
  grid->bc.value = 0.0;
}

void set_neumann_boundary_conditions(grid_t *grid) {
  if (NULL == grid) {
    fprintf(stderr, "\e[1;31m ERROR\e[0m"
	    ": NULL grid pointer when setting BC.\n");
    return;
  }

  grid->bc.type = BC_NEUMANN;
  grid->bc.value = 0.0;
}


void set_dirichlet_boundary_conditions(grid_t *grid, f64 val) {
  if (NULL == grid) {
    fprintf(stderr, "\e[1;31m ERROR\e[0m"
	    ": NULL grid pointer when setting BC.\n");
    return;
  }

  grid->bc.type = BC_DIRICHLET;
  grid->bc.value = val;
}

/* ------------------------------------------------------------------ */

void compute_k2(grid_t *grid) {

  register f64 dkx, dky, dkz = 0.0;
  register f64 kx, ky, kz = 0.0;
  register uint64_t Nx, Ny, Nz = 0;
  
  if (1 == grid->dim) {
    dkx = 2.0*M_PI/grid->L[0];
    Nx = grid->N[0];
    
    for (uint64_t i = 0; i < grid->N[0]; i++) {
      kx = (i <= 0.5*Nx) ? (dkx*i) :
	(dkx*((int32_t)(i-Nx))); grid->k2[i] = kx*kx; // be aware of 0 value later on
    }
  }

  if (2 == grid->dim) {
    dkx = 2.0*M_PI/grid->L[0];
    dky = 2.0*M_PI/grid->L[1];
    Nx = grid->N[0];
    Ny = grid->N[1];
    
    for (uint64_t j = 0; j < Ny; j++) {
      ky = (j <= 0.5*Ny) ? (dky*j) : (dky*((int32_t)(j-Ny)));
      for (uint64_t i = 0; j < Nx; i++) {
	f64 kx = (i <= 0.5*Nx) ? (dkx*i) : (dkx*((int32_t)(i-Nx)));
	grid->k2[i] = kx*kx + ky*ky;
      }
    }
  }

  if (3 == grid->dim) {
    dkx = 2.0*M_PI/grid->L[0];
    dky = 2.0*M_PI/grid->L[1];
    dkx = 2.0*M_PI/grid->L[2];
    Nx = grid->N[0];
    Ny = grid->N[1];
    Nz = grid->N[2];

    for (uint64_t k = 0; k < Nz; k++) {
      kz = (k <= 0.5*Nz) ? (dkz*k) : (dkz*((int32_t)(k-Nz)));
      for (uint64_t j = 0; j < Ny; j++) {
	ky = (j <= 0.5*Ny) ? (dky*j) : (dky*((int32_t)(j-Ny)));
	for (uint64_t i = 0; j < Nx; i++) {
	  kx = (i <= 0.5*Nx) ? (dkx*i) : (dkx*((int32_t)(i-Nx)));
	  grid->k2[i] = kx*kx + ky*ky + kz*kz;
	}
      }
    }
  }
}

