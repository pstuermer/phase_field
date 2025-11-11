#ifndef EQUATION_H
#define EQUATION_H

#include "types.h"
#include "memory_management.h"
#include <fftw3.h>

// Forward declarations
struct grid_t;
struct equation_ftable_t;

/**
 * @enum equation_type_t
 * @brief Identifies the type of equation being solved
 * 
 * Used for type checking and dispatching to equation-specific implementations.
 */
typedef enum {
  EQUATION_NONE,           ///< No equation assigned
  EQUATION_CAHN_HILLIARD   ///< Cahn-Hilliard equation for phase separation
} equation_type_t;

/**
 * @enum time_integration_t
 * @brief Time integration schemes available
 */
typedef enum {
  TIME_NONE,              ///< No time integration
  TIME_SEMI_IMPLICIT      ///< Semi-implicit spectral method (linear implicit, nonlinear explicit)
} time_integration_t;

/**
 * @struct equation_t
 * @brief Generic equation container - holds equation-agnostic data and interfaces to specific implementations
 * 
 * This structure uses a vtable pattern (function pointers) to enable polymorphic behavior.
 * Different equation types (Cahn-Hilliard, Allen-Cahn, etc.) implement the same interface.
 */
typedef struct equation_t {
  equation_type_t type;                  ///< Identifies equation type (for runtime checking)
  const struct equation_ftable_t *ftable; ///< Virtual function table for equation-specific operations
  void *data;                            ///< Opaque pointer to equation-specific data structure
  
  struct grid_t *grid;                   ///< Spatial grid (owned by equation)
  
  // Time integration parameters
  f64 dt;                                ///< Time step size
  uint64_t max_iter;                     ///< Maximum number of iterations
  uint64_t iter;                         ///< Current iteration number
  time_integration_t time_method;        ///< Time integration method
} equation_t;

/**
 * @struct equation_ftable_t
 * @brief Virtual function table - defines interface that all equations must implement
 * 
 * This enables polymorphism in C: different equation types provide their own
 * implementations of these functions, called through the same interface.
 */
typedef struct equation_ftable_t {
  void (*step)(equation_t *eq);            ///< Advance one time step
  void (*setup_spectral)(equation_t *eq);  ///< Precompute spectral operators (k-space)
  void (*cleanup)(equation_t *eq);         ///< Free equation-specific memory
  f64* (*get_field)(equation_t *eq);       ///< Get pointer to primary field
  f64 (*get_free_energy)(equation_t *eq);   ///< Get value of current free energy
  void (*set_field)(equation_t *eq, uint64_t i, f64 val); ///< Set field value at index
} equation_ftable_t;

// ============================================================================
// Lifecycle Management
// ============================================================================

/**
 * @brief Internal destructor - use equation_destroy() macro instead
 * @param eq Pointer to pointer to equation (set to NULL after freeing)
 */
void equation_destroy_internal(equation_t **eq);

/**
 * @def equation_destroy
 * @brief Safely destroy equation and free all memory
 * 
 * Calls equation-specific cleanup, destroys grid, and frees equation structure.
 * Sets pointer to NULL after freeing.
 */
#define equation_destroy(eq) equation_destroy_internal((equation_t **) &(eq))

// ============================================================================
// Configuration
// ============================================================================

/**
 * @brief Set time integration parameters
 * @param eq Equation to configure
 * @param dt Time step size
 * @param max_iter Maximum number of iterations
 * @param iter Starting iteration number (usually 0)
 */
void set_iter(equation_t *eq, f64 dt, uint64_t max_iter, uint64_t iter);

/**
 * @brief Configure semi-implicit time integration and rebuild spectral operators
 * @param eq Equation to configure
 * 
 * Recomputes spectral operators (e.g., exp(-k^4*dt)) based on current dt.
 */
void set_semi_implicit_prop(equation_t *eq);

// ============================================================================
// Field Initialization
// ============================================================================

/**
 * @brief Initialize field with uniform value
 * @param eq Equation whose field to initialize
 * @param value Constant value for entire field
 */
void equation_init_uniform(equation_t *eq, f64 value);

/**
 * @brief Initialize field with random noise
 * @param eq Equation whose field to initialize
 * @param amplitude Random values in range [-amplitude, +amplitude]
 * 
 * Uses rand(), so call srand() first to set seed.
 */
void equation_init_random(equation_t *eq, f64 amplitude);

/**
 * @brief Initialize field with custom function
 * @param eq Equation whose field to initialize
 * @param init_func Function f(x,y,z) that returns field value at coordinates
 * 
 * Example: f64 my_init(f64 *coords) { return sin(coords[0]) * cos(coords[1]); }
 */
void equation_init_custom(equation_t *eq, f64 (*init_func)(f64 *coords));

// ============================================================================
// Time Integration
// ============================================================================

/**
 * @brief Advance equation one time step (calls vtable function)
 * @param eq Equation to advance
 */
void time_prop(equation_t *eq);

/**
 * @brief Run simulation for max_iter steps with periodic console output
 * @param eq Equation to solve
 */
void equation_run(equation_t *eq);

// ============================================================================
// Analysis
// ============================================================================

/**
 * @brief Compute total mass (integral of concentration field)
 * @param eq Equation to analyze
 * @return Total mass = sum(c) over all grid points
 * 
 * For Cahn-Hilliard, mass should be conserved (constant in time).
 */
f64 equation_compute_mass(equation_t *eq);

// ============================================================================
// Output
// ============================================================================

/**
 * @brief Write time and free energy to CSV file
 * @param eq Equation to output
 * @param filename CSV file path
 * @param append If true, append to existing file; if false, overwrite
 */
void equation_output_csv(equation_t *eq, const char *filename, uint8_t append);

/**
 * @brief Write field data in binary format
 * @param eq Equation to output
 * @param filename Output file path
 */
void equation_output_field(equation_t *eq, const char *filename);

#endif // EQUATION_H
