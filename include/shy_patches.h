#ifndef SHY_PATCHES_H
#define SHY_PATCHES_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <functional>
#include <string>

namespace Step85
{
  using namespace dealii;

  template <int dim, int spacedim>
  void
  make_shy_vertex_patches(
    SparsityPattern                 &block_list,
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    const std::function<bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
      &cell_is_in_domain,
    unsigned int                     shyness = 3);

  template <int dim, int spacedim>
  void
  output_patches_vtk(const SparsityPattern               &block_list,
                     const DoFHandler<dim, spacedim>     &dof_handler,
                     const unsigned int                   level,
                     const std::string                   &filename);

  /**
   * Create vertex patches where each patch includes all DoFs for which the
   * residual can be fully evaluated.
   *
   * For each vertex in the domain, this function creates a patch that includes
   * all DoFs from cells touching that vertex, but only if all cells using each
   * DoF are part of the vertex patch. This ensures that the residual
   * contribution of each DoF in the patch can be fully computed using only
   * cells in the patch.
   *
   * @param block_list The output sparsity pattern where each row represents a
   *        patch and contains the DoF indices in that patch.
   * @param dof_handler The DoF handler for the triangulation.
   * @param level The multigrid level to build patches for.
   * @param cell_is_in_domain A predicate function that returns true if a cell
   *        should be considered part of the domain. Only cells for which this
   *        returns true are included in the patches.
   * @param shyness Minimum number of cells a vertex must have to form a patch.
   *        If set to numbers::invalid_unsigned_int (default), all vertices
   *        with at least one cell can form patches.
   */
  template <int dim, int spacedim>
  void
  make_full_residual_vertex_patches(
    SparsityPattern                 &block_list,
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    const std::function<bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
      &cell_is_in_domain,
    unsigned int shyness = numbers::invalid_unsigned_int);

} // namespace Step85

#endif // SHY_PATCHES_H
