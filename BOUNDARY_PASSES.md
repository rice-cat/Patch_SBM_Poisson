# Boundary Passes Feature

## Overview

Both `make_shy_vertex_patches` and `make_full_residual_vertex_patches` now support an additional parameter `n_boundary_passes` to control how many times boundary patches are applied in the smoother.

## Rationale

Patches that touch the boundary of the computational domain often have different properties than interior patches. Specifically, they typically have fewer degrees of freedom than a standard interior patch, which has `(2p-1)^dim` DoFs where `p` is the polynomial degree and `dim` is the spatial dimension.

For better convergence in multigrid methods, it can be beneficial to apply boundary patches multiple times compared to interior patches.

## Implementation

### Boundary Patch Identification

A patch is identified as a boundary patch if it does not have the standard number of DoFs:
- Standard DoF count: `(2*degree - 1)^dim`
- Boundary patches: Any patch with DoF count ≠ standard count

### Patch Duplication

When `n_boundary_passes > 1`, boundary patches are duplicated in the block list. For example, with `n_boundary_passes = 3`, each boundary patch appears 3 times in the patch list, meaning the smoother will apply it 3 times per smoothing step.

## Usage

### Function Signatures

```cpp
template <int dim, int spacedim>
void make_shy_vertex_patches(
  SparsityPattern                 &block_list,
  const DoFHandler<dim, spacedim> &dof_handler,
  const unsigned int               level,
  const std::function<bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
    &cell_is_in_domain,
  unsigned int                     shyness = 3,
  unsigned int                     n_boundary_passes = 1);  // NEW PARAMETER

template <int dim, int spacedim>
void make_full_residual_vertex_patches(
  SparsityPattern                 &block_list,
  const DoFHandler<dim, spacedim> &dof_handler,
  const unsigned int               level,
  const std::function<bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
    &cell_is_in_domain,
  unsigned int shyness = numbers::invalid_unsigned_int,
  unsigned int n_boundary_passes = 1);  // NEW PARAMETER
```

### Example

```cpp
// Create shy patches with boundary patches applied 3 times
SparsityPattern block_list;
make_shy_vertex_patches(
  block_list,
  dof_handler,
  level,
  cell_is_in_domain,
  3,    // shyness
  3);   // n_boundary_passes - boundary patches applied 3x

// Create full residual patches with default behavior (boundary patches applied once)
SparsityPattern block_list_default;
make_full_residual_vertex_patches(
  block_list_default,
  dof_handler,
  level,
  cell_is_in_domain);  // Uses defaults: shyness=invalid, n_boundary_passes=1
```

## Testing

A comprehensive test (`tests/boundary_passes.cc`) verifies that:
1. Boundary patches are correctly identified
2. The number of patches increases with `n_boundary_passes > 1`
3. Both shy and full residual patches work correctly with the new parameter

Run tests with:
```bash
cd build
ctest --output-on-failure
```

## Docker Build

The project includes Docker support for easy building and testing:

```bash
# Build the Docker image
docker build -t patch_sbm_poisson:latest .

# Run tests in Docker
docker run --rm patch_sbm_poisson:latest bash -c "cd build && ctest --output-on-failure"
```

## Default Behavior

The default value of `n_boundary_passes = 1` maintains backward compatibility - existing code will work without changes.
