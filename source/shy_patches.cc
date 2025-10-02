#include "shy_patches.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace Step85
{
  using namespace dealii;

  template <int dim, int spacedim>
  void
  make_shy_vertex_patches(
    SparsityPattern                 &block_list,
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    unsigned int shyness
  )
  {
    AssertDimension(dim, 2);

    bool is_active         = level == numbers::invalid_unsigned_int;
    using patch_index_type = unsigned int;

    patch_index_type i = 0;
    std::map<types::global_vertex_index, patch_index_type>
      index_from_global_vertex_to_patch;


    std::vector<bool> vertex_on_patch(
      dof_handler.get_triangulation().n_vertices(), false);
    std::vector<unsigned int> vertex_cell_count(
      dof_handler.get_triangulation().n_vertices());
    std::vector<bool> vertex_is_sheltered;
    vertex_is_sheltered.resize(dof_handler.get_triangulation().n_vertices());

    std::vector<std::set<dealii::CellId>> cells_of_vertex;
    cells_of_vertex.resize(dof_handler.get_triangulation().n_vertices());

    for (unsigned int ii = 0; ii < vertex_cell_count.size(); ii++)
      {
        vertex_cell_count[ii] = 0;
      }

    for (const auto &cell : dof_handler.cell_iterators_on_level(level) |
                              IteratorFilters::UserFlagSet())
      {
        for (const auto vertex : GeometryInfo<dim>::vertex_indices())
          {
            vertex_cell_count[cell->vertex_index(vertex)]++;
            cells_of_vertex[cell->vertex_index(vertex)].insert(cell->id());
          }
      }

    for (unsigned int ii = 0; ii < vertex_is_sheltered.size(); ii++)
      {
        if (vertex_cell_count[ii] >= shyness)
          vertex_is_sheltered[ii] = true;
        else
          vertex_is_sheltered[ii] = false;
      }

    for (const auto &cell : dof_handler.cell_iterators_on_level(level) |
                              IteratorFilters::UserFlagSet())
      {
        for (const auto &face : cell->face_iterators())
          {
            for (const auto vertex : GeometryInfo<dim - 1>::vertex_indices())
              {
                types::global_vertex_index global_vertex_index =
                  face->vertex_index(vertex);
                if (vertex_is_sheltered[global_vertex_index] == true)
                  {
                    if (vertex_on_patch[global_vertex_index] == false)
                      {
                        vertex_on_patch[global_vertex_index] = true;
                        index_from_global_vertex_to_patch.emplace(
                          global_vertex_index, i);
                        i++;
                      }
                  }
              }
          }
      }


    std::vector<types::global_cell_index> candidate_vertices;
    std::vector<unsigned int>             cells_shared;
    std::set<types::global_vertex_index>  modified_vertices;
    while (true)
      {
        modified_vertices.clear();
        for (const auto &cell : dof_handler.cell_iterators_on_level(level) |
                                  IteratorFilters::UserFlagSet())
          {
            for (const auto &face : cell->face_iterators())
              {
                for (const auto vertex : GeometryInfo<1>::vertex_indices())
                  {
                    types::global_vertex_index global_vertex_index =
                      face->vertex_index(vertex);
                    if (vertex_on_patch[global_vertex_index] == false)
                      {
                        types::global_vertex_index other_vertex =
                          face->vertex_index((vertex + 1) % 2);
                        if (vertex_on_patch[other_vertex])
                          {
                            if (modified_vertices.count(other_vertex) == 0)
                              {
                                bool patch_contains_all_cells = true;
                                for (dealii::CellId cell_index :
                                     cells_of_vertex[global_vertex_index])
                                  {
                                    if (cells_of_vertex[other_vertex].find(
                                          cell_index) ==
                                        cells_of_vertex[other_vertex].end())
                                      {
                                        patch_contains_all_cells = false;
                                        break;
                                      }
                                  }

                                if (patch_contains_all_cells)
                                  {
                                    vertex_on_patch[global_vertex_index] = true;
                                    modified_vertices.insert(
                                      global_vertex_index);
                                    index_from_global_vertex_to_patch.emplace(
                                      global_vertex_index,
                                      index_from_global_vertex_to_patch
                                        [other_vertex]);
                                  }
                              }
                          }
                      }
                  }
              }
          }
        if (modified_vertices.empty())
          break;
      }

    bool all_vertices_on_patch = true;
    for (const auto &cell : dof_handler.cell_iterators_on_level(level) |
                              IteratorFilters::UserFlagSet())
      for (const auto vertex : GeometryInfo<dim>::vertex_indices())
        if (vertex_on_patch[cell->vertex_index(vertex)] == false)
          {
            all_vertices_on_patch = false;
            break;
          }

    std::cout << "All vertices on patch  = " << all_vertices_on_patch
              << std::endl;

    if (all_vertices_on_patch == false)
      {
        throw std::runtime_error(" Not all vertices have a patch :( ");
      }

    std::vector<std::set<types::global_dof_index>> patches_indices(i);
    std::vector<types::global_dof_index>           object_dofs;
    const FiniteElement<dim>                      &fe = dof_handler.get_fe(0);
    types::fe_index                                fe_index;

    object_dofs.reserve(fe.n_dofs_per_cell());

    if (dim == 2)
      {
        unsigned int n_dof_quad = dof_handler.get_fe().n_dofs_per_quad();
        unsigned int n_dof_line   = dof_handler.get_fe().n_dofs_per_line();
        unsigned int n_dof_vertex = dof_handler.get_fe().n_dofs_per_vertex();
        for (const auto &cell : dof_handler.cell_iterators_on_level(level) |
                                  IteratorFilters::UserFlagSet())
          {
            fe_index = cell->active_fe_index();
            for (const auto &face : cell->face_iterators())
              {
                object_dofs.clear();
                for (unsigned int j = 0; j < n_dof_line; j++)
                  if (is_active)
                    object_dofs.push_back(face->dof_index(j));
                  else
                    object_dofs.push_back(face->mg_dof_index(level, j));
                for (const auto vertex :
                     GeometryInfo<dim - 1>::vertex_indices())
                  {
                    unsigned int vertex_global_index =
                      face->vertex_index(vertex);
                    unsigned int vertex_patch_index =
                      index_from_global_vertex_to_patch[vertex_global_index];
                    patches_indices[vertex_patch_index].insert(
                      object_dofs.begin(), object_dofs.end());
                  }
              }

            object_dofs.clear();
            for (unsigned int j = 0; j < n_dof_quad; j++)
              object_dofs.push_back(cell->mg_dof_index(level, j));
            for (const auto vertex : GeometryInfo<dim>::vertex_indices())
              {
                unsigned int vertex_global_index = cell->vertex_index(vertex);
                unsigned int vertex_patch_index =
                  index_from_global_vertex_to_patch[vertex_global_index];
                for (const auto dof : object_dofs)
                  patches_indices[vertex_patch_index].insert(dof);
                for (unsigned int j = 0; j < n_dof_vertex; j++)
                  {
                    patches_indices[vertex_patch_index].insert(
                      cell->mg_vertex_dof_index(level, vertex, j, fe_index));
                  }
              }
          }

        block_list.reinit(patches_indices.size(),
                          dof_handler.n_dofs(level),
                          dof_handler.get_fe().n_dofs_per_cell() *
                            std::pow(2, dim));
        for (i = 0; i < patches_indices.size(); i++)
          {
            for (const auto dof_index : patches_indices[i])
              {
                block_list.add(i, dof_index);
              }
          }
      }
    std::cout << "Displaying Block_list : \n";
    block_list.print(std::cout);
    std::cout << "\nEnd Block_list \n\n";
    block_list.compress();
  }

  template void
  make_shy_vertex_patches<2, 2>(SparsityPattern          &block_list,
                                 const DoFHandler<2, 2> &dof_handler,
                                 const unsigned int       level,
                                 unsigned int             shyness);

} // namespace Step85
