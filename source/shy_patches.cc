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
    const std::function<
      bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
                &cell_is_in_domain,
    unsigned int shyness,
    unsigned int n_boundary_passes)
  {
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

    // Iterate over cells - use active_cell_iterators if level is invalid
    auto process_cells = [&](const auto &cell_range) {
      for (const auto &cell : cell_range)
        {
          if (!cell_is_in_domain(cell))
            continue;

          for (const auto vertex : GeometryInfo<dim>::vertex_indices())
            {
              vertex_cell_count[cell->vertex_index(vertex)]++;
              cells_of_vertex[cell->vertex_index(vertex)].insert(cell->id());
            }
        }
    };

    if (is_active)
      process_cells(dof_handler.active_cell_iterators());
    else
      process_cells(dof_handler.cell_iterators_on_level(level));

    for (unsigned int ii = 0; ii < vertex_is_sheltered.size(); ii++)
      {
        if (vertex_cell_count[ii] >= shyness)
          vertex_is_sheltered[ii] = true;
        else
          vertex_is_sheltered[ii] = false;
      }

    // Iterate over cells to find sheltered vertices
    auto find_sheltered_vertices = [&](const auto &cell_range) {
      for (const auto &cell : cell_range)
        {
          if (!cell_is_in_domain(cell))
            continue;

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
    };

    if (is_active)
      find_sheltered_vertices(dof_handler.active_cell_iterators());
    else
      find_sheltered_vertices(dof_handler.cell_iterators_on_level(level));


    std::vector<types::global_cell_index> candidate_vertices;
    std::vector<unsigned int>             cells_shared;
    std::set<types::global_vertex_index>  modified_vertices;
    while (true)
      {
        modified_vertices.clear();

        auto expand_patches = [&](const auto &cell_range) {
          for (const auto &cell : cell_range)
            {
              if (!cell_is_in_domain(cell))
                continue;

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
                                      vertex_on_patch[global_vertex_index] =
                                        true;
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
        };

        if (is_active)
          expand_patches(dof_handler.active_cell_iterators());
        else
          expand_patches(dof_handler.cell_iterators_on_level(level));

        if (modified_vertices.empty())
          break;
      }

    bool all_vertices_on_patch = true;

    auto check_all_vertices = [&](const auto &cell_range) {
      for (const auto &cell : cell_range)
        {
          if (!cell_is_in_domain(cell))
            continue;

          for (const auto vertex : GeometryInfo<dim>::vertex_indices())
            if (vertex_on_patch[cell->vertex_index(vertex)] == false)
              {
                all_vertices_on_patch = false;
                break;
              }
          if (!all_vertices_on_patch)
            break;
        }
    };

    if (is_active)
      check_all_vertices(dof_handler.active_cell_iterators());
    else
      check_all_vertices(dof_handler.cell_iterators_on_level(level));

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
        unsigned int n_dof_quad   = dof_handler.get_fe().n_dofs_per_quad();
        unsigned int n_dof_line   = dof_handler.get_fe().n_dofs_per_line();
        unsigned int n_dof_vertex = dof_handler.get_fe().n_dofs_per_vertex();

        auto collect_dofs = [&](const auto &cell_range) {
          for (const auto &cell : cell_range)
            {
              if (!cell_is_in_domain(cell))
                continue;

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
                if (is_active)
                  object_dofs.push_back(cell->dof_index(j));
                else
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
                      if (is_active)
                        patches_indices[vertex_patch_index].insert(
                          cell->vertex_dof_index(vertex, j, fe_index));
                      else
                        patches_indices[vertex_patch_index].insert(
                          cell->mg_vertex_dof_index(
                            level, vertex, j, fe_index));
                    }
                }
            }
        };


        if (is_active)
          collect_dofs(dof_handler.active_cell_iterators());
        else
          collect_dofs(dof_handler.cell_iterators_on_level(level));

        const bool force_outside_dofs_to_singleton_patches = false;
        if (force_outside_dofs_to_singleton_patches)
          {
            std::vector<bool> dof_assigned_to_patch(is_active ?
                                                      dof_handler.n_dofs() :
                                                      dof_handler.n_dofs(level),
                                                    false);

            for (i = 0; i < patches_indices.size(); i++)
              {
                for (const auto dof_index : patches_indices[i])
                  {
                    dof_assigned_to_patch[dof_index] = true;
                  }
              }

            for (size_t i = 0; i < dof_assigned_to_patch.size(); i++)
              {
                if (dof_assigned_to_patch[i] == false)
                  {
                    patches_indices.push_back(
                      std::set<types::global_dof_index>());
                    patches_indices.back().insert(i);
                  }
              }
          }

        // Identify boundary patches based on DoF count
        // Standard patches have (2*degree - 1)^dim DoFs
        const unsigned int fe_degree = fe.degree;
        const unsigned int base = 2 * fe_degree - 1;
        const unsigned int standard_dof_count = (dim == 2) ? base * base : base * base * base;
        
        std::vector<char> is_boundary_patch(patches_indices.size(), false);
        for (i = 0; i < patches_indices.size(); i++)
          {
            if (patches_indices[i].size() != standard_dof_count)
              {
                is_boundary_patch[i] = true;
              }
          }

        // Calculate total number of patch entries (including duplicates for boundary patches)
        // First pass includes all patches, subsequent passes only include boundary patches
        size_t num_boundary_patches = 0;
        for (i = 0; i < patches_indices.size(); i++)
          {
            if (is_boundary_patch[i])
              num_boundary_patches++;
          }
        size_t total_patch_entries = patches_indices.size() + 
                                      num_boundary_patches * (n_boundary_passes - 1);

        block_list.reinit(total_patch_entries,
                          is_active ? dof_handler.n_dofs() :
                                      dof_handler.n_dofs(level),
                          dof_handler.get_fe().n_dofs_per_cell() *
                            (dim == 2 ? 4 : 8));


        // Add patches to block_list with correct ordering:
        // First pass: all patches (A, B, C, D)
        // Subsequent passes: only boundary patches (B, C, B, C, ...)
        size_t row_index = 0;
        
        // First pass: add all patches
        for (i = 0; i < patches_indices.size(); i++)
          {
            for (const auto dof_index : patches_indices[i])
              {
                block_list.add(row_index, dof_index);
              }
            row_index++;
          }
        
        // Subsequent passes: add only boundary patches
        for (unsigned int pass = 1; pass < n_boundary_passes; pass++)
          {
            for (i = 0; i < patches_indices.size(); i++)
              {
                if (is_boundary_patch[i])
                  {
                    for (const auto dof_index : patches_indices[i])
                      {
                        block_list.add(row_index, dof_index);
                      }
                    row_index++;
                  }
              }
          }
      }
    block_list.compress();
  }

  template void
  make_shy_vertex_patches<DEAL_II_DIMENSION, DEAL_II_DIMENSION>(
    SparsityPattern        &block_list,
    const DoFHandler<DEAL_II_DIMENSION, DEAL_II_DIMENSION> &dof_handler,
    const unsigned int      level,
    const std::function<bool(const typename DoFHandler<DEAL_II_DIMENSION, DEAL_II_DIMENSION>::cell_iterator &)>
                &cell_is_in_domain,
    unsigned int shyness,
    unsigned int n_boundary_passes);

  template <int dim, int spacedim>
  void
  make_full_residual_vertex_patches(
    SparsityPattern                 &block_list,
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    const std::function<
      bool(const typename DoFHandler<dim, spacedim>::cell_iterator &)>
                &cell_is_in_domain,
    unsigned int shyness,
    unsigned int n_boundary_passes)
  {
    using patch_index_type = unsigned int;
    using cell_iterator    = typename DoFHandler<dim, spacedim>::cell_iterator;

    // Step 1: Build mapping from DoF indices to cells that use them
    std::map<types::global_dof_index, std::set<dealii::CellId>> cells_of_dof;

    // Step 2: Build mapping from vertices to cells that use them
    std::map<types::global_vertex_index, std::set<dealii::CellId>>
      cells_of_vertex;

    // Step 3: Build mapping from CellId to cell iterator for efficient lookup
    std::map<dealii::CellId, cell_iterator> cell_map;

    const FiniteElement<dim>            &fe = dof_handler.get_fe(0);
    std::vector<types::global_dof_index> local_dof_indices(
      fe.n_dofs_per_cell());

    bool is_active = level == numbers::invalid_unsigned_int;

    // Iterate over all cells on the level (or active cells) that are in the
    // domain
    auto process_cells = [&](const auto &cell_range) {
      for (const auto &cell : cell_range)
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
          if (is_active)
            cell->get_dof_indices(local_dof_indices);
          else
            cell->get_mg_dof_indices(local_dof_indices);
          for (const auto dof_index : local_dof_indices)
            {
              cells_of_dof[dof_index].insert(cell->id());
            }
        }
    };

    if (is_active)
      process_cells(dof_handler.active_cell_iterators());
    else
      process_cells(dof_handler.cell_iterators_on_level(level));

    // Step 4: Create patches - one patch per vertex (filtered by shyness)
    std::map<types::global_vertex_index, patch_index_type>
                     vertex_to_patch_index;
    patch_index_type patch_count = 0;

    // Assign a patch index to each vertex that has enough cells (shyness check)
    for (const auto &vertex_cells_pair : cells_of_vertex)
      {
        // If shyness is invalid_unsigned_int, accept all vertices
        // Otherwise, only accept vertices with at least 'shyness' cells
        if (shyness == numbers::invalid_unsigned_int ||
            vertex_cells_pair.second.size() >= shyness)
          {
            vertex_to_patch_index[vertex_cells_pair.first] = patch_count++;
          }
      }

    // Step 5: For each vertex patch, collect DoFs where all using cells are in
    // the patch
    std::vector<std::set<types::global_dof_index>> patches_dofs(patch_count);

    for (const auto &vertex_cells_pair : cells_of_vertex)
      {
        types::global_vertex_index      vertex_index = vertex_cells_pair.first;
        const std::set<dealii::CellId> &vertex_cells = vertex_cells_pair.second;

        // Skip vertices that didn't pass the shyness filter
        if (vertex_to_patch_index.find(vertex_index) ==
            vertex_to_patch_index.end())
          continue;

        patch_index_type patch_idx = vertex_to_patch_index[vertex_index];

        // Collect all DoFs from cells touching this vertex
        std::set<types::global_dof_index> candidate_dofs;
        for (const auto &cell_id : vertex_cells)
          {
            // Use the cell_map for O(log n) lookup instead of O(n) search
            const auto &cell = cell_map[cell_id];
            if (is_active)
              cell->get_dof_indices(local_dof_indices);
            else
              cell->get_mg_dof_indices(local_dof_indices);
            for (const auto dof_index : local_dof_indices)
              {
                candidate_dofs.insert(dof_index);
              }
          }

        // For each candidate DoF, check if all cells using it are in the vertex
        // patch
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

    // Step 6: Identify boundary patches based on DoF count
    // Standard patches have (2*degree - 1)^dim DoFs
    const unsigned int fe_degree = fe.degree;
    const unsigned int base = 2 * fe_degree - 1;
    const unsigned int standard_dof_count = (dim == 2) ? base * base : base * base * base;
    
    std::vector<char> is_boundary_patch(patch_count, false);
    for (patch_index_type i = 0; i < patch_count; i++)
      {
        if (patches_dofs[i].size() != standard_dof_count)
          {
            is_boundary_patch[i] = true;
          }
      }

    // Calculate total number of patch entries (including duplicates for boundary patches)
    // First pass includes all patches, subsequent passes only include boundary patches
    size_t num_boundary_patches = 0;
    for (patch_index_type i = 0; i < patch_count; i++)
      {
        if (is_boundary_patch[i])
          num_boundary_patches++;
      }
    size_t total_patch_entries = patch_count + 
                                  num_boundary_patches * (n_boundary_passes - 1);

    // Step 7: Build the sparsity pattern
    block_list.reinit(total_patch_entries,
                      is_active ? dof_handler.n_dofs() :
                                  dof_handler.n_dofs(level),
                      fe.n_dofs_per_cell() * (dim == 2 ? 4 : 8));

    // Add patches to block_list with correct ordering:
    // First pass: all patches (A, B, C, D)
    // Subsequent passes: only boundary patches (B, C, B, C, ...)
    size_t row_index = 0;
    
    // First pass: add all patches
    for (patch_index_type i = 0; i < patch_count; i++)
      {
        for (const auto dof_index : patches_dofs[i])
          {
            block_list.add(row_index, dof_index);
          }
        row_index++;
      }
    
    // Subsequent passes: add only boundary patches
    for (unsigned int pass = 1; pass < n_boundary_passes; pass++)
      {
        for (patch_index_type i = 0; i < patch_count; i++)
          {
            if (is_boundary_patch[i])
              {
                for (const auto dof_index : patches_dofs[i])
                  {
                    block_list.add(row_index, dof_index);
                  }
                row_index++;
              }
          }
      }

    block_list.compress();
  }

  template void
  make_full_residual_vertex_patches<DEAL_II_DIMENSION, DEAL_II_DIMENSION>(
    SparsityPattern        &block_list,
    const DoFHandler<DEAL_II_DIMENSION, DEAL_II_DIMENSION> &dof_handler,
    const unsigned int      level,
    const std::function<bool(const typename DoFHandler<DEAL_II_DIMENSION, DEAL_II_DIMENSION>::cell_iterator &)>
                &cell_is_in_domain,
    unsigned int shyness,
    unsigned int n_boundary_passes);

  template <int dim, int spacedim>
  void
  output_patches_vtk(const SparsityPattern           &block_list,
                     const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int               level,
                     const std::string               &filename)
  {
    // Get support points for DoFs on the given level
    const MappingQ1<dim, spacedim>                     mapping;
    std::map<types::global_dof_index, Point<spacedim>> dof_location_map;

    // Manually get support points for the multigrid level
    const FiniteElement<dim, spacedim> &fe            = dof_handler.get_fe();
    const unsigned int                  dofs_per_cell = fe.n_dofs_per_cell();

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
             it != block_list.end(patch);
             ++it)
          {
            ++total_points;
          }
      }

    // Write points
    vtk_file << "POINTS " << total_points << " double\n";
    for (unsigned int patch = 0; patch < block_list.n_rows(); ++patch)
      {
        for (SparsityPattern::iterator it = block_list.begin(patch);
             it != block_list.end(patch);
             ++it)
          {
            const types::global_dof_index dof_index = it->column();
            const Point<spacedim>        &point = dof_location_map[dof_index];

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
             it != block_list.end(patch);
             ++it)
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
             it != block_list.end(patch);
             ++it)
          {
            vtk_file << patch << "\n";
          }
      }

    vtk_file.close();
    std::cout << "Patches output to VTK file: " << filename << std::endl;
  }

  template void
  output_patches_vtk<DEAL_II_DIMENSION, DEAL_II_DIMENSION>(
    const SparsityPattern  &block_list,
    const DoFHandler<DEAL_II_DIMENSION, DEAL_II_DIMENSION> &dof_handler,
    const unsigned int      level,
    const std::string      &filename);

} // namespace Step85
