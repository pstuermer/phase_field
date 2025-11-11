# Phase-Field Simulation Playground
A modular playground for simulating microstructure evolution using phase-fields

## Overview
Implements (currently) the Cahn-Hilliard equation for spinodal decomposotion using:
- Spectral methods (FFTW) for spatial derivatives
- semi-implicit time stepping
- extensible framework via function pointer tables
- free energy tracking for validation (spatial derivative needs fixing)
- mass conservation verification
- NIST benchmark 1b

## Installation

### Prerequisites

- C11 compiler (gcc/clang)
- FFTW3 library
- Make

### Build
```bash
# Clone repository
git clone 
cd phase-field-framework

# Compile
make
(Don't forget to link fftw3 appropiately if not installed globally)

include/
  types.h              - Basic type definitions
  memory_management.h  - Memory allocation with alignment
  cartesian.h          - Grid and boundary conditions
  equation.h           - Generic equation interface
  cahn_hilliard.h      - Cahn-Hilliard specific code
src/
  cartesian.c          - Grid implementation, k² computation
  equation.c           - Generic functions (init, run, output)
  cahn_hilliard.c      - Semi-implicit spectral solver
tests/
  benchmark1b.c        - NIST benchmark driver
```


### Design Decisions

**Modular architecture:** Equations interface through vtables (function pointers), allowing new equation types without modifying core code.

**Spectral methods:** FFTW for spatial derivatives. Uses DCT-II/III for Neumann (no-flux) boundary conditions.

**Semi-implicit time integration:** 
```
ĉ^(n+1) = [ĉ^n - Mdt·k²·F{f'(c^n)}] / [1 + Mdt·κ·k⁴]
```
Linear ∇⁴ term treated implicitly (stable), nonlinear term explicit. Precompute `1/(1+Mκk⁴dt)` for efficiency.

**Memory management:** Custom allocator with 64-byte alignment for SIMD optimization.



### Assumptions and shortcut
- assume cartesian grid for all systems
- currently only no-flux (von-Neumann) boundary conditions, but easily extendable
- no adaptive timestepping
- currently no unittests
- only Cahn-Hilliard equation
