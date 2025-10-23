// Test for shy_patches functionality
// Creates a simple triangulation, distributes dofs, and creates shy patches

#include "tests.h"
#include "shy_patches.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparsity_pattern.h>

using namespace dealii;

int
main()
{
  initlog();
  
  const int dim = 2;
  
  // Create a simple triangulation with proper smoothing
  Triangulation<dim> triangulation(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.refine_global(2);
  
  // Create DoFHandler with FE_Q(2)
  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();
  
  // Mark all cells with user flag (required by make_shy_vertex_patches)
  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_user_flag();
  
  // Create shy patches
  SparsityPattern block_list;
  const unsigned int level = 2;
  const unsigned int shyness = 3;
  
  deallog << "Creating shy patches for level " << level 
          << " with shyness " << shyness << std::endl;
  deallog << "Number of cells: " << triangulation.n_cells() << std::endl;
  deallog << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
  deallog << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
  deallog << "Number of DoFs at level " << level << ": " 
          << dof_handler.n_dofs(level) << std::endl;
  
  Step85::make_shy_vertex_patches(block_list, dof_handler, level, shyness);
  
  deallog << "Shy patches created successfully" << std::endl;
  
  return 0;
}
