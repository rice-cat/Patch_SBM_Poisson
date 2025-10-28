#include "step-85.h"

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "shy_patches.h"

namespace Step85
{
  using namespace dealii;

  template <int dim>
  LaplaceSolver<dim>::LaplaceSolver(const Settings &settings)
    : fe_degree(settings.fe_degree)
    , rhs_function(0)
    , boundary_condition(0.0)
    , triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , fe_level_set(fe_degree)
    , level_set_dof_handler(triangulation)
    , level_set()
    , fe_poisson(fe_degree)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
    , settings(settings)
  {}

  template <int dim>
  void
  LaplaceSolver<dim>::make_grid()
  {
    std::cout << "Creating background mesh" << std::endl;

    GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
    triangulation.refine_global(2);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::setup_discrete_level_set()
  {
    std::cout << "Setting up discrete level set function" << std::endl;

    level_set_dof_handler.distribute_dofs(fe_level_set);
    level_set.reinit(level_set_dof_handler.n_dofs());

    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;
    VectorTools::interpolate(level_set_dof_handler,
                             signed_distance_sphere,
                             level_set);
  }

  enum ActiveFEIndex
  {
    lagrange = 0,
    nothing  = 1
  };

  template <int dim>
  void
  LaplaceSolver<dim>::distribute_dofs()
  {
    std::cout << "Distributing degrees of freedom" << std::endl;


    triangulation.clear_user_flags();
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const NonMatching::LocationToLevelSet cell_location =
          mesh_classifier.location_to_level_set(cell);

        // only langrage elements are inside
        if (cell_location != NonMatching::LocationToLevelSet::inside)
          {}
        else
          cell->set_user_flag();
      }

    dof_handler.distribute_dofs(fe_poisson);
    dof_handler.distribute_mg_dofs();
  }

  template <int dim>
  void
  LaplaceSolver<dim>::initialize_matrices()
  {
    std::cout << "Initializing matrices" << std::endl;

    const auto face_has_flux_coupling = [&](const auto &       cell,
                                            const unsigned int face_index) {
      return this->face_has_ghost_penalty(cell, face_index);
    };

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    const unsigned int n_components = 1; // fe_collection.n_components();
    Table<2, DoFTools::Coupling> cell_coupling(n_components, n_components);
    Table<2, DoFTools::Coupling> face_coupling(n_components, n_components);
    cell_coupling[0][0] = DoFTools::always;
    face_coupling[0][0] = DoFTools::always;

    const AffineConstraints<double> constraints;
    const bool                      keep_constrained_dofs = true;

    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         constraints,
                                         keep_constrained_dofs,
                                         cell_coupling,
                                         face_coupling,
                                         numbers::invalid_subdomain_id,
                                         face_has_flux_coupling);
    sparsity_pattern.copy_from(dsp);

    stiffness_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    rhs.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  bool
  LaplaceSolver<dim>::face_has_ghost_penalty(
    const typename Triangulation<dim>::active_cell_iterator &cell,
    const unsigned int                                       face_index) const
  {
    if (cell->at_boundary(face_index))
      return false;

    const NonMatching::LocationToLevelSet cell_location =
      mesh_classifier.location_to_level_set(cell);

    const NonMatching::LocationToLevelSet neighbor_location =
      mesh_classifier.location_to_level_set(cell->neighbor(face_index));

    if (cell_location == NonMatching::LocationToLevelSet::intersected &&
        neighbor_location != NonMatching::LocationToLevelSet::outside)
      return true;

    if (neighbor_location == NonMatching::LocationToLevelSet::intersected &&
        cell_location != NonMatching::LocationToLevelSet::outside)
      return true;

    return false;
  }

  template <int dim>
  void
  LaplaceSolver<dim>::assemble_system()
  {
    Step85::assemble_system(dof_handler,
                            fe_poisson,
                            mesh_classifier,
                            fe_degree,
                            rhs_function,
                            boundary_condition,
                            stiffness_matrix,
                            rhs,
                            active_dofs);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::setup_smoother()
  {
    unsigned int level = triangulation.n_levels() - 1;
    make_shy_vertex_patches(smoother_data.block_list, dof_handler, level);
    smoother_data.relaxation = 1.;
    smoother_data.inversion  = PreconditionBlockBase<double>::svd;
    auto smoother            = std::make_unique<SmootherType>();
    smoother->initialize(stiffness_matrix, smoother_data);
    patch_smoother = std::move(smoother);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::write_dof_locations(const std::string &filename)
  {
    std::map<types::global_dof_index, Point<dim>> dof_location_map =
      DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler);


    for (auto map_element = dof_location_map.begin();
         map_element != dof_location_map.end();)
      {
        if (active_dofs[map_element->first] == false)
          {
            auto temp = map_element;
            map_element++;
            dof_location_map.erase(temp);
          }
        else
          {
            map_element++;
          }
      }



    std::ofstream dof_location_file(filename);
    DoFTools::write_gnuplot_dof_support_point_info(dof_location_file,
                                                   dof_location_map);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::solve()
  {
    std::random_device rd;
    std::mt19937       e2(rd());
    e2.seed(85);
    std::uniform_real_distribution<double> dist(-1, 1);

    for (unsigned int i = 0; i < solution.size(); i++)
      {
        if (active_dofs[i])
          solution[i] = dist(e2);
        // if (i < 20)
        //   std::cout << solution[i] << std::endl;
      }

    std::cout << "Before : " << solution.norm_sqr() << std::endl;
    for (unsigned int i = 0; i < 5; i++)
      patch_smoother->step(solution, rhs);
    std::cout << "After  : " << solution.norm_sqr() << std::endl << std::endl;


    // SparseDirectUMFPACK sparse_direct;
    // sparse_direct.initialize(stiffness_matrix);
    // sparse_direct.vmult(solution, rhs);
    // return;
    // const unsigned int max_iterations = solution.size();
    // SolverControl      solver_control(max_iterations);
    // // SolverCG<>         solver(solver_control);
    // SolverGMRES<> solver(solver_control);
    // solver.solve(stiffness_matrix, solution, rhs, PreconditionIdentity());
  }

  template <int dim>
  void
  LaplaceSolver<dim>::output_results() const
  {
    std::cout << "Writing vtu file" << std::endl;

    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");
    // data_out.add_data_vector()


    data_out.build_patches();
    std::ofstream output("patch_smoother.vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  double
  AnalyticalSolution<dim>::value(const Point<dim> & point,
                                 const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;

    return 1. - 2. / dim * (point.norm_square() - 1.);
  }

  template <int dim>
  double
  LaplaceSolver<dim>::compute_L2_error() const
  {
    std::cout << "Computing L2 error" << std::endl;

    const QGauss<dim> quadrature_formula(fe_degree + 1);

    // NonMatching::RegionUpdateFlags region_update_flags;
    // region_update_flags.inside =
    //   update_values | update_JxW_values | update_quadrature_points;

    FEValues<dim> cell_fe_values(fe_poisson,
                                 quadrature_formula,
                                 update_values | update_JxW_values |
                                   update_quadrature_points);

    AnalyticalSolution<dim> analytical_solution;
    double                  error_L2_squared = 0;

    for (const auto &cell :
         dof_handler.active_cell_iterators() | IteratorFilters::UserFlagSet())
      {
        cell_fe_values.reinit(cell);


        std::vector<double> solution_values(cell_fe_values.n_quadrature_points);
        cell_fe_values.get_function_values(solution, solution_values);

        for (const unsigned int q : cell_fe_values.quadrature_point_indices())
          {
            const Point<dim> &point = cell_fe_values.quadrature_point(q);
            const double      error_at_point =
              solution_values.at(q) - analytical_solution.value(point);
            error_L2_squared +=
              std::pow(error_at_point, 2) * cell_fe_values.JxW(q);
          }
      }

    return std::sqrt(error_L2_squared);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::run()
  {
    ConvergenceTable   convergence_table;
    const unsigned int n_refinements = 1;

    make_grid();
    for (unsigned int cycle = 0; cycle <= n_refinements; cycle++)
      {
        std::cout << "Refinement cycle " << cycle << std::endl;
        triangulation.refine_global(1);
        setup_discrete_level_set();
        std::cout << "Classifying cells" << std::endl;
        mesh_classifier.reclassify();
        distribute_dofs();
        std::ofstream mesh_file("mesh.gnuplot");
        GridOut().write_gnuplot(triangulation, mesh_file);

        initialize_matrices();
        assemble_system();

        write_dof_locations("dof-locations-1.gnuplot");
        setup_smoother();
        solve(); //"SOLVE does smoother steps"
                 // if (cycle == 1)
        output_results();
        const double error_L2 = compute_L2_error();
        const double cell_side_length =
          triangulation.begin_active()->minimum_vertex_distance();

        convergence_table.add_value("Cycle", cycle);
        convergence_table.add_value("Mesh size", cell_side_length);
        convergence_table.add_value("L2-Error", error_L2);

        convergence_table.evaluate_convergence_rates(
          "L2-Error", ConvergenceTable::reduction_rate_log2);
        convergence_table.set_scientific("L2-Error", true);

        std::cout << std::endl;
        convergence_table.write_text(std::cout);
        std::cout << std::endl;
      }
  }

  // Standalone assembly function for reuse in smoother testing and MG solver
  template <int dim>
  void
  assemble_system(const DoFHandler<dim>                     &dof_handler,
                  const FE_Q<dim>                           &fe_poisson,
                  const NonMatching::MeshClassifier<dim>    &mesh_classifier,
                  const unsigned int                         fe_degree,
                  const Functions::ConstantFunction<dim>    &rhs_function,
                  const Functions::ConstantFunction<dim>    &boundary_condition,
                  SparseMatrix<double>                      &stiffness_matrix,
                  Vector<double>                            &rhs,
                  std::vector<bool>                         &active_dofs)
  {
    std::cout << "Assembling" << std::endl;

    active_dofs.clear();
    active_dofs.resize(dof_handler.n_dofs(), false);

    const unsigned int n_dofs_per_cell = fe_poisson.dofs_per_cell;
    FullMatrix<double> local_stiffness(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double>     local_rhs(n_dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

    const double nitsche_parameter = 5 * (fe_degree + 1) * fe_degree;

    const QGauss<dim - 1> face_quadrature(fe_degree + 1);

    FEFaceValues<dim> face_fe_values(fe_poisson,
                                     face_quadrature,
                                     update_values | update_gradients |
                                       update_normal_vectors |
                                       update_JxW_values |
                                       update_quadrature_points);

    const QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> cell_fe_values(fe_poisson,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_JxW_values |
                                   update_quadrature_points);


    for (const auto &cell :
         dof_handler.active_cell_iterators() | IteratorFilters::UserFlagSet())
      {
        local_stiffness = 0;
        local_rhs       = 0;

        const double cell_side_length = cell->minimum_vertex_distance();
        const double alpha            = nitsche_parameter / cell_side_length;
        cell_fe_values.reinit(cell);

        for (const unsigned int q : cell_fe_values.quadrature_point_indices())
          {
            const Point<dim> &point = cell_fe_values.quadrature_point(q);
            for (const unsigned int i : cell_fe_values.dof_indices())
              {
                for (const unsigned int j : cell_fe_values.dof_indices())
                  {
                    local_stiffness(i, j) += cell_fe_values.shape_grad(i, q) *
                                             cell_fe_values.shape_grad(j, q) *
                                             cell_fe_values.JxW(q);
                  }
                local_rhs(i) += rhs_function.value(point) *
                                cell_fe_values.shape_value(i, q) *
                                cell_fe_values.JxW(q);
              }
          }

        for (size_t face_index = 0; face_index < cell->n_faces(); face_index++)
          {
            bool face_on_boundary = false;
            if (cell->neighbor_index(face_index) != -1)
              {
                const NonMatching::LocationToLevelSet neighbor_location =
                  mesh_classifier.location_to_level_set(
                    cell->neighbor(face_index));
                if (neighbor_location ==
                    NonMatching::LocationToLevelSet::intersected)
                  face_on_boundary = true;
              }
            else
              {
                face_on_boundary = true;
              }
            if (face_on_boundary) // face is on the shifted boundary
              {
                cell->face(face_index)->set_user_flag();
                face_fe_values.reinit(cell, face_index);
                for (const unsigned int q :
                     face_fe_values.quadrature_point_indices())
                  {
                    const Point<dim> &point =
                      face_fe_values.quadrature_point(q);
                    const Point<dim> unit_point =
                      cell_fe_values.get_mapping().transform_real_to_unit_cell(
                        cell, point);
                    const Point<dim> closest_boundary_point =
                      point / point.norm();
                    const Point<dim> unit_boundary_point =
                      cell_fe_values.get_mapping().transform_real_to_unit_cell(
                        cell, closest_boundary_point);

                    std::vector<double> shifted_shape_values(n_dofs_per_cell);
                    std::vector<double> shape_values(n_dofs_per_cell);
                    std::vector<Tensor<1, dim>> shape_grad(n_dofs_per_cell);

                    const FiniteElement<dim> &fe_lagrange = fe_poisson;
                    for (size_t i = 0; i < n_dofs_per_cell; i++)
                      {
                        shifted_shape_values[i] =
                          fe_lagrange.shape_value(i, unit_boundary_point);
                        shape_values[i] =
                          fe_lagrange.shape_value(i, unit_point);
                        shape_grad[i] = fe_lagrange.shape_grad(i, unit_point);
                      }

                    for (size_t i = 0; i < n_dofs_per_cell; i++)
                      {
                        for (size_t j = 0; j < n_dofs_per_cell; j++)
                          {
                            local_stiffness(i, j) -=
                              face_fe_values.normal_vector(q) *
                              face_fe_values.shape_grad(i, q) *
                              shifted_shape_values[j] * face_fe_values.JxW(q);
                            local_stiffness(i, j) -=
                              face_fe_values.normal_vector(q) *
                              face_fe_values.shape_grad(j, q) *
                              shape_values[i] * face_fe_values.JxW(q);

                            local_stiffness(i, j) +=
                              alpha * shifted_shape_values[j] *
                              shifted_shape_values[i] * face_fe_values.JxW(q);
                          }
                        local_rhs[i] -=
                          boundary_condition.value(closest_boundary_point) *
                          face_fe_values.normal_vector(q) *
                          face_fe_values.shape_grad(i, q) *
                          face_fe_values.JxW(q);
                        local_rhs[i] +=
                          alpha *
                          boundary_condition.value(closest_boundary_point) *
                          shifted_shape_values[i] * face_fe_values.JxW(q);
                      }
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);
        for (const auto index : local_dof_indices)
          active_dofs[index] = true;

        stiffness_matrix.add(local_dof_indices, local_stiffness);
        rhs.add(local_dof_indices, local_rhs);
      }
  }

  // Explicit template instantiations
  template class LaplaceSolver<2>;
  template class AnalyticalSolution<2>;
  template void assemble_system<2>(const DoFHandler<2>                     &,
                                   const FE_Q<2>                           &,
                                   const NonMatching::MeshClassifier<2>    &,
                                   const unsigned int                       ,
                                   const Functions::ConstantFunction<2>    &,
                                   const Functions::ConstantFunction<2>    &,
                                   SparseMatrix<double>                    &,
                                   Vector<double>                          &,
                                   std::vector<bool>                       &);

} // namespace Step85
