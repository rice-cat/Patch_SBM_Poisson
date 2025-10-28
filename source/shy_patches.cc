#include "shy_patches.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <fstream>
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

  template <int dim, int spacedim>
  void
  make_full_residual_vertex_patches(
    SparsityPattern                 &block_list,
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    const std::function<bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
      &cell_is_in_domain)
  {
    AssertDimension(dim, 2);

    using patch_index_type = unsigned int;
    using cell_iterator = typename DoFHandler<dim, spacedim>::cell_iterator;

    // Step 1: Build mapping from DoF indices to cells that use them
    std::map<types::global_dof_index, std::set<dealii::CellId>> cells_of_dof;
    
    // Step 2: Build mapping from vertices to cells that use them
    std::map<types::global_vertex_index, std::set<dealii::CellId>> cells_of_vertex;
    
    // Step 3: Build mapping from CellId to cell iterator for efficient lookup
    std::map<dealii::CellId, cell_iterator> cell_map;
    
    const FiniteElement<dim> &fe = dof_handler.get_fe(0);
    std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());

    // Iterate over all cells on the level that are in the domain
    for (const auto &cell : dof_handler.cell_iterators_on_level(level))
      {
        if (!cell_is_in_domain(cell))
          continue;
        
        // Store the cell iterator for later use
        cell_map[cell->id()] = cell;
          
        // Add this cell to the vertex mapping
        for (const auto vertex : GeometryInfo<dim>::vertex_indices())
          {
            cells_of_vertex[cell->vertex_index(vertex)].insert(cell->id());
          }
        
        // Add this cell to the DoF mapping
        cell->get_mg_dof_indices(local_dof_indices);
        for (const auto dof_index : local_dof_indices)
          {
            cells_of_dof[dof_index].insert(cell->id());
          }
      }

    // Step 4: Create patches - one patch per vertex
    std::map<types::global_vertex_index, patch_index_type> vertex_to_patch_index;
    patch_index_type patch_count = 0;
    
    // Assign a patch index to each vertex that has cells
    for (const auto &vertex_cells_pair : cells_of_vertex)
      {
        vertex_to_patch_index[vertex_cells_pair.first] = patch_count++;
      }

    // Step 5: For each vertex patch, collect DoFs where all using cells are in the patch
    std::vector<std::set<types::global_dof_index>> patches_dofs(patch_count);
    
    for (const auto &vertex_cells_pair : cells_of_vertex)
      {
        types::global_vertex_index vertex_index = vertex_cells_pair.first;
        const std::set<dealii::CellId> &vertex_cells = vertex_cells_pair.second;
        patch_index_type patch_idx = vertex_to_patch_index[vertex_index];
        
        // Collect all DoFs from cells touching this vertex
        std::set<types::global_dof_index> candidate_dofs;
        for (const auto &cell_id : vertex_cells)
          {
            // Use the cell_map for O(log n) lookup instead of O(n) search
            const auto &cell = cell_map[cell_id];
            cell->get_mg_dof_indices(local_dof_indices);
            for (const auto dof_index : local_dof_indices)
              {
                candidate_dofs.insert(dof_index);
              }
          }
        
        // For each candidate DoF, check if all cells using it are in the vertex patch
        for (const auto dof_index : candidate_dofs)
          {
            const std::set<dealii::CellId> &dof_cells = cells_of_dof[dof_index];
            
            // Check if all cells using this DoF are in the vertex patch
            bool all_cells_in_patch = true;
            for (const auto &cell_id : dof_cells)
              {
                if (vertex_cells.find(cell_id) == vertex_cells.end())
                  {
                    all_cells_in_patch = false;
                    break;
                  }
              }
            
            if (all_cells_in_patch)
              {
                patches_dofs[patch_idx].insert(dof_index);
              }
          }
      }

    // Step 6: Build the sparsity pattern
    block_list.reinit(patch_count,
                      dof_handler.n_dofs(level),
                      fe.n_dofs_per_cell() * std::pow(2, dim));
    
    for (patch_index_type i = 0; i < patch_count; i++)
      {
        for (const auto dof_index : patches_dofs[i])
          {
            block_list.add(i, dof_index);
          }
      }

    std::cout << "Displaying Block_list for full residual patches: \n";
    block_list.print(std::cout);
    std::cout << "\nEnd Block_list \n\n";
    block_list.compress();
  }

  template void
  make_full_residual_vertex_patches<2, 2>(
    SparsityPattern          &block_list,
    const DoFHandler<2, 2> &dof_handler,
    const unsigned int       level,
    const std::function<bool(const typename DoFHandler<2, 2>::cell_iterator &)>
      &cell_is_in_domain);

  template <int dim, int spacedim>
  void
  output_patches_vtk(const SparsityPattern               &block_list,
                     const DoFHandler<dim, spacedim>     &dof_handler,
                     const unsigned int                   level,
                     const std::string                   &filename)
  {
    // Get support points for DoFs on the given level
    const MappingQ1<dim, spacedim> mapping;
    std::map<types::global_dof_index, Point<spacedim>> dof_location_map;
    
    // Manually get support points for the multigrid level
    const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    
    for (const auto &cell : dof_handler.cell_iterators_on_level(level))
      {
        cell->get_mg_dof_indices(local_dof_indices);
        
        const std::vector<Point<dim>> &unit_support_points = 
          fe.get_unit_support_points();
          
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const Point<spacedim> support_point = 
              mapping.transform_unit_to_real_cell(cell, unit_support_points[i]);
            dof_location_map[local_dof_indices[i]] = support_point;
          }
      }

    std::ofstream vtk_file(filename);
    
    // Write VTK header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Shy Patches Visualization\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // Count total points (DoFs can appear in multiple patches)
    unsigned int total_points = 0;
    for (unsigned int patch = 0; patch < block_list.n_rows(); ++patch)
      {
        for (SparsityPattern::iterator it = block_list.begin(patch);
             it != block_list.end(patch); ++it)
          {
            ++total_points;
          }
      }

    // Write points
    vtk_file << "POINTS " << total_points << " double\n";
    for (unsigned int patch = 0; patch < block_list.n_rows(); ++patch)
      {
        for (SparsityPattern::iterator it = block_list.begin(patch);
             it != block_list.end(patch); ++it)
          {
            const types::global_dof_index dof_index = it->column();
            const Point<spacedim> &point = dof_location_map[dof_index];
            
            vtk_file << point[0] << " " << point[1];
            if (spacedim == 3)
              vtk_file << " " << point[2];
            else
              vtk_file << " 0.0";
            vtk_file << "\n";
          }
      }

    // Write cells (each point is a vertex cell)
    vtk_file << "\nCELLS " << total_points << " " << (total_points * 2) << "\n";
    for (unsigned int i = 0; i < total_points; ++i)
      vtk_file << "1 " << i << "\n";

    // Write cell types (VTK_VERTEX = 1)
    vtk_file << "\nCELL_TYPES " << total_points << "\n";
    for (unsigned int i = 0; i < total_points; ++i)
      vtk_file << "1\n";

    // Write point data
    vtk_file << "\nPOINT_DATA " << total_points << "\n";
    
    // DoF indices as scalar data
    vtk_file << "SCALARS dof_index int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (unsigned int patch = 0; patch < block_list.n_rows(); ++patch)
      {
        for (SparsityPattern::iterator it = block_list.begin(patch);
             it != block_list.end(patch); ++it)
          {
            vtk_file << it->column() << "\n";
          }
      }

    // Patch indices as scalar data
    vtk_file << "\nSCALARS patch_index int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (unsigned int patch = 0; patch < block_list.n_rows(); ++patch)
      {
        for (SparsityPattern::iterator it = block_list.begin(patch);
             it != block_list.end(patch); ++it)
          {
            vtk_file << patch << "\n";
          }
      }

    vtk_file.close();
    std::cout << "Patches output to VTK file: " << filename << std::endl;
  }

  template void
  output_patches_vtk<2, 2>(const SparsityPattern     &block_list,
                           const DoFHandler<2, 2>   &dof_handler,
                           const unsigned int         level,
                           const std::string         &filename);

} // namespace Step85
