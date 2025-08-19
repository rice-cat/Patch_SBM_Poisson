

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>

namespace Step85
{
  using namespace dealii;

  template <int dim>
  class LaplaceSolver
  {
  public:
    LaplaceSolver();

    void
    run();

  private:
    void
    make_grid();

    void
    setup_discrete_level_set();

    void
    distribute_dofs();

    void
    initialize_matrices();

    void
    assemble_system();

    void
    solve();

    void
    output_results() const;

    double
    compute_L2_error() const;

    bool
    face_has_ghost_penalty(
      const typename Triangulation<dim>::active_cell_iterator &cell,
      const unsigned int face_index) const;

    const unsigned int fe_degree;

    const Functions::ConstantFunction<dim> rhs_function;
    const Functions::ConstantFunction<dim> boundary_condition;

    Triangulation<dim> triangulation;

    const FE_Q<dim> fe_level_set;
    DoFHandler<dim> level_set_dof_handler;
    Vector<double>  level_set;

    hp::FECollection<dim> fe_collection;
    DoFHandler<dim>       dof_handler;
    Vector<double>        solution;

    NonMatching::MeshClassifier<dim> mesh_classifier;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    Vector<double>       rhs;
  };

  template <int dim>
  LaplaceSolver<dim>::LaplaceSolver()
    : fe_degree(2)
    , rhs_function(4.0)
    , boundary_condition(1.0)
    , fe_level_set(fe_degree)
    , level_set_dof_handler(triangulation)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
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

    fe_collection.push_back(FE_Q<dim>(fe_degree));
    fe_collection.push_back(FE_Nothing<dim>());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const NonMatching::LocationToLevelSet cell_location =
          mesh_classifier.location_to_level_set(cell);

        // This is another place where code is changed from step-85, now the
        // only langrage elements are inside
        if (cell_location != NonMatching::LocationToLevelSet::inside)
          cell->set_active_fe_index(ActiveFEIndex::nothing);
        else
          cell->set_active_fe_index(ActiveFEIndex::lagrange);
      }

    dof_handler.distribute_dofs(fe_collection);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::initialize_matrices()
  {
    std::cout << "Initializing matrices" << std::endl;

    const auto face_has_flux_coupling = [&](const auto        &cell,
                                            const unsigned int face_index) {
      return this->face_has_ghost_penalty(cell, face_index);
    };

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    const unsigned int           n_components = fe_collection.n_components();
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
    std::cout << "Assembling" << std::endl;

    const unsigned int n_dofs_per_cell = fe_collection[0].dofs_per_cell;
    FullMatrix<double> local_stiffness(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double>     local_rhs(n_dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

    // const double ghost_parameter   = 0.5;
    const double nitsche_parameter = 5 * (fe_degree + 1) * fe_degree;

    const QGauss<dim - 1> face_quadrature(fe_degree + 1);
    // FEInterfaceValues<dim> fe_interface_values(fe_collection[0],
    //                                            face_quadrature,
    //                                            update_gradients |
    //                                              update_JxW_values |
    //                                              update_normal_vectors);

    FEFaceValues<dim> face_fe_values(fe_collection[0],
                                     face_quadrature,
                                     update_values | update_gradients |
                                       update_normal_vectors |
                                       update_JxW_values |
                                       update_quadrature_points);

    const QGauss<1> quadrature_1D(fe_degree + 1);

    NonMatching::RegionUpdateFlags region_update_flags;
    region_update_flags.inside = update_values | update_gradients |
                                 update_JxW_values | update_quadrature_points;
    region_update_flags.surface = update_values | update_gradients |
                                  update_JxW_values | update_quadrature_points |
                                  update_normal_vectors;

    NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                      quadrature_1D,
                                                      region_update_flags,
                                                      mesh_classifier,
                                                      level_set_dof_handler,
                                                      level_set);

    // uuh nvm this part
    // const QGauss<2> quadrature_2D(fe_degree + 1);
    // FEValues<dim>   ordinary_fe_values(fe_collection[0],
    //                                  quadrature_2D,
    //                                  update_values | update_gradients |
    //                                    update_JxW_values |
    //                                    update_quadrature_points);

    for (const auto &cell :
         dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
        local_stiffness = 0;
        local_rhs       = 0;

        const double cell_side_length = cell->minimum_vertex_distance();
        const double alpha            = nitsche_parameter / cell_side_length;
        non_matching_fe_values.reinit(cell);

        const std_cxx17::optional<FEValues<dim>> &inside_fe_values =
          non_matching_fe_values.get_inside_fe_values();

        // if (inside_fe_values) // I dont think i need this anymore, as all
        //                       // lagrange elements are now "inside"
        for (const unsigned int q :
             inside_fe_values->quadrature_point_indices())
          {
            const Point<dim> &point = inside_fe_values->quadrature_point(q);
            for (const unsigned int i : inside_fe_values->dof_indices())
              {
                for (const unsigned int j : inside_fe_values->dof_indices())
                  {
                    local_stiffness(i, j) +=
                      inside_fe_values->shape_grad(i, q) *
                      inside_fe_values->shape_grad(j, q) *
                      inside_fe_values->JxW(q);
                  }
                local_rhs(i) += rhs_function.value(point) *
                                inside_fe_values->shape_value(i, q) *
                                inside_fe_values->JxW(q);
              }
          }

        for (size_t face_index = 0; face_index < cell->n_faces(); face_index++)
          {
            const NonMatching::LocationToLevelSet neighbor_location =
              mesh_classifier.location_to_level_set(cell->neighbor(
                face_index)); // Q: is face->index correct replacement for
                              // face_index   A: NO, cell->face_accessors() is
                              // functionally useless?
            if (neighbor_location ==
                NonMatching::LocationToLevelSet::
                  intersected) // face is on the shifted boundary
              {
                face_fe_values.reinit(cell, face_index);
                for (const unsigned int q :
                     face_fe_values.quadrature_point_indices())
                  {
                    const Point<dim> &point =
                      face_fe_values.quadrature_point(q);
                    const Point<dim> unit_point =
                      inside_fe_values->get_mapping()
                        .transform_real_to_unit_cell(cell, point);
                    const Point<dim> closest_boundary_point =
                      point / point.norm();
                    const Point<dim> unit_boundary_point =
                      inside_fe_values->get_mapping()
                        .transform_real_to_unit_cell(cell,
                                                     closest_boundary_point);

                    std::vector<double> shifted_shape_values(n_dofs_per_cell);
                    std::vector<double> shape_values(n_dofs_per_cell);
                    std::vector<Tensor<1, dim>> shape_grad(n_dofs_per_cell);

                    const FiniteElement<dim> &fe_lagrange = fe_collection[0];
                    for (size_t i = 0; i < n_dofs_per_cell; i++)
                      {
                        shifted_shape_values[i] =
                          fe_lagrange.shape_value(i, unit_boundary_point);
                        shape_values[i] =
                          fe_lagrange.shape_value(i, unit_point);
                        shape_grad[i] = fe_lagrange.shape_grad(i, unit_point);
                      }

                    // Q : does the dof indices on a face, include all dofs from
                    // the cell? or is it just the face endpoint dofs?
                    //   because if its only the face endpoint dofs we have a
                    //   problem

                    for (size_t i = 0; i < n_dofs_per_cell; i++)
                      {
                        for (size_t j = 0; j < n_dofs_per_cell; j++)
                          {
                            // Q: assuming j is for u and i is fo v, is this
                            // right ?
                            local_stiffness(i, j) -=
                              face_fe_values.normal_vector(q) *
                              face_fe_values.shape_grad(i, q) *
                              shifted_shape_values[j] * face_fe_values.JxW(q);
                            local_stiffness(i, j) -=
                              face_fe_values.normal_vector(q) *
                              face_fe_values.shape_grad(j, q) *
                              shape_values[i] * face_fe_values.JxW(q);

                            // this term we can shift both terms or one of them,
                            // as long as its consistent with the
                            // alpha/nitsche_parameter term in the rhs, i
                            // decided to shift both
                            local_stiffness(i, j) +=
                              alpha * shifted_shape_values[j] *
                              shape_values[i] * face_fe_values.JxW(q);
                          }
                        local_rhs[i] -=
                          boundary_condition.value(closest_boundary_point) *
                          face_fe_values.normal_vector(q) *
                          face_fe_values.shape_grad(i, q) *
                          face_fe_values.JxW(q);
                        local_rhs[i] +=
                          alpha *
                          boundary_condition.value(closest_boundary_point) *
                          shape_values[i] * face_fe_values.JxW(q);
                      }
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);

        stiffness_matrix.add(local_dof_indices, local_stiffness);
        rhs.add(local_dof_indices, local_rhs);

      }
  }

  template <int dim>
  void
  LaplaceSolver<dim>::solve()
  {
    std::cout << "Solving system" << std::endl;
    SparseDirectUMFPACK sparse_direct;
    sparse_direct.initialize(stiffness_matrix);
    sparse_direct.vmult(solution, rhs);
    return;
    const unsigned int max_iterations = solution.size();
    SolverControl      solver_control(max_iterations);
    // SolverCG<>         solver(solver_control);
    SolverGMRES<> solver(solver_control);
    solver.solve(stiffness_matrix, solution, rhs, PreconditionIdentity());
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

    data_out.set_cell_selection(
      [this](const typename Triangulation<dim>::cell_iterator &cell) {
        return cell->is_active() &&
               mesh_classifier.location_to_level_set(cell) ==
                 NonMatching::LocationToLevelSet::inside;
      });

    data_out.build_patches();
    std::ofstream output("step-85.vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  class AnalyticalSolution : public Function<dim>
  {
  public:
    double
    value(const Point<dim>  &point,
          const unsigned int component = 0) const override;
  };

  template <int dim>
  double
  AnalyticalSolution<dim>::value(const Point<dim>  &point,
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

    const QGauss<1> quadrature_1D(fe_degree + 1);

    NonMatching::RegionUpdateFlags region_update_flags;
    region_update_flags.inside =
      update_values | update_JxW_values | update_quadrature_points;

    NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                      quadrature_1D,
                                                      region_update_flags,
                                                      mesh_classifier,
                                                      level_set_dof_handler,
                                                      level_set);

    AnalyticalSolution<dim> analytical_solution;
    double                  error_L2_squared = 0;

    for (const auto &cell :
         dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
        non_matching_fe_values.reinit(cell);

        const std_cxx17::optional<FEValues<dim>> &fe_values =
          non_matching_fe_values.get_inside_fe_values();

        if (fe_values)
          {
            std::vector<double> solution_values(fe_values->n_quadrature_points);
            fe_values->get_function_values(solution, solution_values);

            for (const unsigned int q : fe_values->quadrature_point_indices())
              {
                const Point<dim> &point = fe_values->quadrature_point(q);
                const double      error_at_point =
                  solution_values.at(q) - analytical_solution.value(point);
                error_L2_squared +=
                  std::pow(error_at_point, 2) * fe_values->JxW(q);
              }
          }
      }

    return std::sqrt(error_L2_squared);
  }

  template <int dim>
  void
  LaplaceSolver<dim>::run()
  {
    ConvergenceTable   convergence_table;
    const unsigned int n_refinements = 5;

    make_grid();
    for (unsigned int cycle = 0; cycle <= n_refinements; cycle++)
      {
        std::cout << "Refinement cycle " << cycle << std::endl;
        triangulation.refine_global(1);
        setup_discrete_level_set();
        std::cout << "Classifying cells" << std::endl;
        mesh_classifier.reclassify();
        distribute_dofs();
        initialize_matrices();
        assemble_system();
        solve();
        if (cycle == 1)
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

} // namespace Step85

int
main()
{
  const int dim = 2;

  Step85::LaplaceSolver<dim> laplace_solver;
  laplace_solver.run();
}
