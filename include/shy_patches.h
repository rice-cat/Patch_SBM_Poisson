#ifndef SHY_PATCHES_H
#define SHY_PATCHES_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <string>

namespace Step85
{
  using namespace dealii;

  template <int dim, int spacedim>
  void
  make_shy_vertex_patches(SparsityPattern                 &block_list,
                          const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int               level,
                          unsigned int                     shyness = 3);

  template <int dim, int spacedim>
  void
  output_patches_vtk(const SparsityPattern               &block_list,
                     const DoFHandler<dim, spacedim>     &dof_handler,
                     const unsigned int                   level,
                     const std::string                   &filename);

} // namespace Step85

#endif // SHY_PATCHES_H
