// Test for boundary patch identification with n_boundary_passes parameter
// Creates a triangulation with boundary and interior patches,
// verifies that boundary patches are duplicated n_boundary_passes times

#include "shy_patches.h"

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <iostream>
#include <sstream>

#include "tests.h"

using namespace dealii;

int
main()
{
  initlog();

  const int dim = 2;

  // Create a simple triangulation with proper smoothing
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.refine_global(2);

  // Create DoFHandler with FE_Q(2)
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  const unsigned int level = triangulation.n_levels() - 1;

  // Mark cells on the target level with user flag based on distance to origin
  // Create a subdomain that excludes some boundary cells
  for (const auto &cell : dof_handler.cell_iterators_on_level(level))
    {
      if (cell->center().norm() < 1.2)
        cell->set_user_flag();
      else
        cell->clear_user_flag();
    }

  // Create a predicate function that uses user flags
  auto cell_is_in_domain =
    [](const typename DoFHandler<dim>::cell_iterator &cell) -> bool {
    return cell->user_flag_set();
  };

  // Test 1: Create shy patches with n_boundary_passes = 1 (default)
  SparsityPattern    block_list_1;
  const unsigned int shyness = 3;

  deallog << "Test 1: Creating shy patches with n_boundary_passes = 1"
          << std::endl;
  Step85::make_shy_vertex_patches(
    block_list_1, dof_handler, level, cell_is_in_domain, shyness, 1);

  deallog << "Number of patches with n_boundary_passes = 1: "
          << block_list_1.n_rows() << std::endl;

  // Test 2: Create shy patches with n_boundary_passes = 3
  SparsityPattern block_list_3;

  deallog << "Test 2: Creating shy patches with n_boundary_passes = 3"
          << std::endl;
  Step85::make_shy_vertex_patches(
    block_list_3, dof_handler, level, cell_is_in_domain, shyness, 3);

  deallog << "Number of patches with n_boundary_passes = 3: "
          << block_list_3.n_rows() << std::endl;

  // Verify that with n_boundary_passes = 3, we have more patches due to duplication
  deallog << "Patches increased due to boundary duplication: "
          << (block_list_3.n_rows() > block_list_1.n_rows() ? "YES" : "NO")
          << std::endl;

  // Test 3: Create full residual patches with n_boundary_passes = 1
  SparsityPattern block_list_fr_1;

  deallog << "Test 3: Creating full residual patches with n_boundary_passes = 1"
          << std::endl;
  Step85::make_full_residual_vertex_patches(
    block_list_fr_1, dof_handler, level, cell_is_in_domain, numbers::invalid_unsigned_int, 1);

  deallog << "Number of full residual patches with n_boundary_passes = 1: "
          << block_list_fr_1.n_rows() << std::endl;

  // Test 4: Create full residual patches with n_boundary_passes = 2
  SparsityPattern block_list_fr_2;

  deallog << "Test 4: Creating full residual patches with n_boundary_passes = 2"
          << std::endl;
  Step85::make_full_residual_vertex_patches(
    block_list_fr_2, dof_handler, level, cell_is_in_domain, numbers::invalid_unsigned_int, 2);

  deallog << "Number of full residual patches with n_boundary_passes = 2: "
          << block_list_fr_2.n_rows() << std::endl;

  // Verify that with n_boundary_passes = 2, we have more patches due to duplication
  deallog << "Full residual patches increased due to boundary duplication: "
          << (block_list_fr_2.n_rows() > block_list_fr_1.n_rows() ? "YES" : "NO")
          << std::endl;

  deallog << "All tests completed successfully" << std::endl;

  return 0;
}
