#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "types.h"
#include "memory_management.h"

typedef enum {
	      BC_PERIODIC,
	      BC_NEUMANN,
	      BC_DIRICHLET
} boundary_condition_t

typedef struct {
  boundary_condition_t type;
  f64 value; // for Dirichlet BC
} boundary_condition;

typedef struct cartesian_system_t cartesian_system_t;

struct cartesian_system_t {

  uint64_t *dim;
  uint64_t *N;
  f64 *L;
}

#endif CARTESIAN_H
