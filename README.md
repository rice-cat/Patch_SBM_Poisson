# Patch SBM Poisson

A deal.II-based application for solving Poisson equations using the Shifted Boundary Method (SBM) with patch-based smoothers.

## Structure

This project follows the standard deal.II application structure:

- `include/` - Header files
  - `step-85.h` - Main LaplaceSolver class declaration and standalone assembly function
  - `shy_patches.h` - Shy vertex patches utility functions
- `source/` - Implementation files
  - `step-85.cc` - LaplaceSolver class implementation and assembly function
  - `shy_patches.cc` - Shy vertex patches implementation
- `execs/` - Executable files with main() functions
  - `smoother_tester.cc` - Tests patch-based smoothers
  - `sbm_mg_solver.cc` - Full multigrid solver with parameter handling
- `build/` - Build directory (not tracked in git)

## Building

To build the project:

```bash
mkdir -p build
cd build
cmake ..
make
```

## Requirements

- CMake >= 3.13.4
- deal.II >= 9.5

## Running

After building, run the executables:

### Smoother Tester
```bash
./smoother_tester
```

### MG Solver
```bash
./sbm_mg_solver [parameter_file.prm]
```

If no parameter file is provided, default parameters will be used. See `mg_solver_params.prm` for an example parameter file.

### Parameter File

The MG solver accepts a parameter file with the following options:
- `omega`: Relaxation parameter for smoother (default: 1.0)
- `shyness`: Shyness parameter for patch construction (default: 3)
- `multiplicative`: Use multiplicative (true) or additive (false) smoother (default: true)
- `smoother_type`: Type of patch smoother - "shy_patches" or "full_residual" (default: shy_patches)
- `fe_degree`: Finite element polynomial degree (default: 2)
- `n_refinements`: Number of global refinements (default: 2)
- `max_iterations`: Maximum solver iterations (default: 100)
- `solver_tolerance`: Solver tolerance (default: 1e-10)

