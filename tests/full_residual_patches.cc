// Test for full_residual_patches functionality
// Creates a simple triangulation, distributes dofs, and creates full residual patches
// This is analogous to shy_patches.cc but uses user flags to determine cell activity

#include "shy_patches.h"

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

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
  triangulation.refine_global(3);

  // Create DoFHandler with FE_Q(2)
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // Mark cells with user flag based on distance to origin
  // Only cells within a certain radius are considered "in the domain"
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->center().norm() < 0.75)
        cell->set_user_flag();
      else
        cell->clear_user_flag();
    }

  // Propagate user flags to all levels
  for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
    {
      for (const auto &cell : dof_handler.cell_iterators_on_level(level))
        {
          // A cell on a coarser level has user flag if any of its children has it
          bool has_active_descendant = false;
          if (cell->has_children())
            {
              for (unsigned int child = 0; child < cell->n_children(); ++child)
                {
                  if (cell->child(child)->user_flag_set())
                    {
                      has_active_descendant = true;
                      break;
                    }
                }
            }
          else
            {
              // Leaf cell - check if it's marked
              has_active_descendant = cell->user_flag_set();
            }
          
          if (has_active_descendant)
            cell->set_user_flag();
          else
            cell->clear_user_flag();
        }
    }

  // Create a predicate function that uses user flags
  auto cell_is_in_domain = 
    [](const typename DoFHandler<dim>::cell_iterator &cell) -> bool {
      return cell->user_flag_set();
    };

  // Create full residual patches
  SparsityPattern    block_list;
  const unsigned int level = triangulation.n_levels() - 1;

  deallog << "Creating full residual patches for level " << level << std::endl;
  deallog << "Number of cells: " << triangulation.n_cells() << std::endl;
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
  deallog << "Number of DoFs at level " << level << ": "
          << dof_handler.n_dofs(level) << std::endl;

  Step85::make_full_residual_vertex_patches(block_list, dof_handler, level, 
                                             cell_is_in_domain);

  // Output patches to VTK file for visualization
  Step85::output_patches_vtk(block_list, dof_handler, level, 
                             "full_residual_patches.vtk");

  // Output mesh with user flags for visualization
  Vector<float> user_flags(triangulation.n_active_cells());
  unsigned int  cell_index = 0;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      user_flags(cell_index) = cell->user_flag_set() ? 1.0 : 0.0;
      ++cell_index;
    }

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  data_out.add_data_vector(user_flags, "user_flags");
  data_out.build_patches();

  std::ofstream output("user_flags_full_residual.vtk");
  deallog << "Full residual patches created successfully" << std::endl;

  {
    std::ostringstream oss;
    block_list.print(oss);
    deallog << oss.str();
    deallog << std::endl;
  }

  return 0;
}
