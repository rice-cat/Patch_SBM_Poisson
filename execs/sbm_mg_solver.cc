#include "step-85.h"
#include "shy_patches.h"

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_base.h>

#include <fstream>
#include <iostream>

namespace Step85
{
  using namespace dealii;

  struct MGParameters
  {
    double       omega               = 1.0;
    unsigned int shyness             = 3;
    bool         multiplicative      = true;
    std::string  smoother_type       = "shy_patches";
    unsigned int fe_degree           = 2;
    unsigned int n_refinements       = 2;
    unsigned int max_iterations      = 100;
    double       solver_tolerance    = 1e-10;

    void
    declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Multigrid parameters");
      {
        prm.declare_entry("omega",
                          "1.0",
                          Patterns::Double(0.0),
                          "Relaxation parameter for smoother");
        prm.declare_entry("shyness",
                          "3",
                          Patterns::Integer(0),
                          "Shyness parameter for patch construction");
        prm.declare_entry("multiplicative",
                          "true",
                          Patterns::Bool(),
                          "Use multiplicative (true) or additive (false) smoother");
        prm.declare_entry("smoother_type",
                          "shy_patches",
                          Patterns::Selection("shy_patches|full_residual"),
                          "Type of patch smoother to use");
        prm.declare_entry("fe_degree",
                          "2",
                          Patterns::Integer(1),
                          "Finite element polynomial degree");
        prm.declare_entry("n_refinements",
                          "2",
                          Patterns::Integer(0),
                          "Number of global refinements");
        prm.declare_entry("max_iterations",
                          "100",
                          Patterns::Integer(1),
                          "Maximum solver iterations");
        prm.declare_entry("solver_tolerance",
                          "1e-10",
                          Patterns::Double(0.0),
                          "Solver tolerance");
      }
      prm.leave_subsection();
    }

    void
    parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Multigrid parameters");
      {
        omega           = prm.get_double("omega");
        shyness         = prm.get_integer("shyness");
        multiplicative  = prm.get_bool("multiplicative");
        smoother_type   = prm.get("smoother_type");
        fe_degree       = prm.get_integer("fe_degree");
        n_refinements   = prm.get_integer("n_refinements");
        max_iterations  = prm.get_integer("max_iterations");
        solver_tolerance = prm.get_double("solver_tolerance");
      }
      prm.leave_subsection();
    }
  };

  template <int dim>
  class MGSolver
  {
  public:
    MGSolver(const MGParameters &params);

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
    setup_multigrid();

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

    const FE_Q<dim>   fe_poisson;
    DoFHandler<dim>   dof_handler;
    Vector<double>    solution;
    std::vector<bool> active_dofs;

    NonMatching::MeshClassifier<dim> mesh_classifier;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    Vector<double>       rhs;

    // Multigrid level objects
    MGLevelObject<SparsityPattern>      mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double>> mg_matrices;
    MGLevelObject<Vector<double>>       mg_rhs;

    const MGParameters mg_params;
  };

  template <int dim>
  MGSolver<dim>::MGSolver(const MGParameters &params)
    : fe_degree(params.fe_degree)
    , rhs_function(0)
    , boundary_condition(0.0)
    , triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , fe_level_set(fe_degree)
    , level_set_dof_handler(triangulation)
    , level_set()
    , fe_poisson(fe_degree)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
    , mg_params(params)
  {}

  template <int dim>
  void
  MGSolver<dim>::make_grid()
  {
    std::cout << "Creating background mesh" << std::endl;
    GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
    triangulation.refine_global(mg_params.n_refinements);
  }

  template <int dim>
  void
  MGSolver<dim>::setup_discrete_level_set()
  {
    std::cout << "Setting up discrete level set function" << std::endl;
    level_set_dof_handler.distribute_dofs(fe_level_set);
    level_set.reinit(level_set_dof_handler.n_dofs());

    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;
    VectorTools::interpolate(level_set_dof_handler,
                             signed_distance_sphere,
                             level_set);
  }

  template <int dim>
  void
  MGSolver<dim>::distribute_dofs()
  {
    std::cout << "Distributing degrees of freedom" << std::endl;

    triangulation.clear_user_flags();
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const NonMatching::LocationToLevelSet cell_location =
          mesh_classifier.location_to_level_set(cell);

        if (cell_location == NonMatching::LocationToLevelSet::inside)
          cell->set_user_flag();
      }

    dof_handler.distribute_dofs(fe_poisson);
    dof_handler.distribute_mg_dofs();
  }

  template <int dim>
  void
  MGSolver<dim>::initialize_matrices()
  {
    std::cout << "Initializing matrices" << std::endl;

    const auto face_has_flux_coupling = [&](const auto &       cell,
                                            const unsigned int face_index) {
      return this->face_has_ghost_penalty(cell, face_index);
    };

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    const unsigned int n_components = 1;
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
  MGSolver<dim>::face_has_ghost_penalty(
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
  MGSolver<dim>::assemble_system()
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
  MGSolver<dim>::setup_multigrid()
  {
    std::cout << "Setting up multigrid" << std::endl;
    
    const unsigned int n_levels = triangulation.n_levels();
    
    // Resize MG level objects for all levels
    mg_sparsity_patterns.resize(0, n_levels - 1);
    mg_matrices.resize(0, n_levels - 1);
    mg_rhs.resize(0, n_levels - 1);
    
    std::cout << "  Number of MG levels: " << n_levels << std::endl;
    std::cout << "  Using " << mg_params.smoother_type << " smoother" << std::endl;
    std::cout << "  Omega: " << mg_params.omega << std::endl;
    std::cout << "  Shyness: " << mg_params.shyness << std::endl;
    std::cout << "  Mode: " << (mg_params.multiplicative ? "multiplicative" : "additive") << std::endl;
    
    // Initialize sparsity patterns and matrices for each level
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
                                   dof_handler.n_dofs(level));
        MGTools::make_sparsity_pattern(dof_handler, dsp, level);
        mg_sparsity_patterns[level].copy_from(dsp);
        mg_matrices[level].reinit(mg_sparsity_patterns[level]);
        mg_rhs[level].reinit(dof_handler.n_dofs(level));
      }
    
    std::cout << "  MG matrices initialized for " << n_levels << " levels" << std::endl;
  }

  template <int dim>
  void
  MGSolver<dim>::solve()
  {
    std::cout << "Solving with MG preconditioner" << std::endl;

    // For now, use simple CG solver
    // TODO: Replace with MG-preconditioned solver
    SolverControl solver_control(mg_params.max_iterations, 
                                   mg_params.solver_tolerance);
    SolverCG<> solver(solver_control);

    try
      {
        solver.solve(stiffness_matrix, solution, rhs, PreconditionIdentity());
        std::cout << "  Converged in " << solver_control.last_step() 
                  << " iterations" << std::endl;
      }
    catch (std::exception &e)
      {
        std::cout << "  Solver failed: " << e.what() << std::endl;
      }
  }

  template <int dim>
  void
  MGSolver<dim>::output_results() const
  {
    std::cout << "Writing vtu file" << std::endl;

    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");

    data_out.build_patches();
    std::ofstream output("mg_solver.vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  double
  MGSolver<dim>::compute_L2_error() const
  {
    std::cout << "Computing L2 error" << std::endl;

    const QGauss<dim> quadrature_formula(fe_degree + 1);

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
  MGSolver<dim>::run()
  {
    make_grid();
    setup_discrete_level_set();
    std::cout << "Classifying cells" << std::endl;
    mesh_classifier.reclassify();
    distribute_dofs();

    initialize_matrices();
    assemble_system();
    setup_multigrid();
    solve();
    output_results();

    const double error_L2 = compute_L2_error();
    std::cout << "\nL2 Error: " << error_L2 << std::endl;
  }

} // namespace Step85

int
main(int argc, char **argv)
{
  const int dim = 2;

  Step85::MGParameters params;
  dealii::ParameterHandler prm;
  params.declare_parameters(prm);

  // Parse command line parameter file if provided
  if (argc > 1)
    {
      try
        {
          prm.parse_input(argv[1]);
        }
      catch (std::exception &e)
        {
          std::cout << "Failed to parse parameter file: " << e.what() << std::endl;
          return 1;
        }
    }
  else
    {
      std::cout << "Using default parameters (no parameter file provided)" << std::endl;
    }

  params.parse_parameters(prm);

  Step85::MGSolver<dim> mg_solver(params);
  mg_solver.run();

  return 0;
}
