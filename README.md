# Patch SBM Poisson

A deal.II-based application for solving Poisson equations using the Shifted Boundary Method (SBM) with patch-based smoothers.

## Structure

This project follows the standard deal.II application structure:

- `include/` - Header files
  - `step-85.h` - Main LaplaceSolver class declaration
  - `shy_patches.h` - Shy vertex patches utility functions
- `source/` - Implementation files
  - `step-85.cc` - LaplaceSolver class implementation
  - `shy_patches.cc` - Shy vertex patches implementation
- `execs/` - Executable files with main() functions
  - `step-85.cc` - Main executable
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

After building, run the executable:

```bash
./step-85
```

or from the build directory:

```bash
cd build
./step-85
```
