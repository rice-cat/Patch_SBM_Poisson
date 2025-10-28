#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <fstream>
#include <iostream>
#include <memory>

#include "shy_patches.h"
#include "step-85.h"

namespace Step85
{
  using namespace dealii;

  struct MGParameters
  {
    double       omega            = 1.0;
    unsigned int shyness          = 3;
    bool         multiplicative   = true;
    std::string  smoother_type    = "shy_patches";
    unsigned int fe_degree        = 2;
    unsigned int n_refinements    = 2;
    unsigned int max_iterations   = 100;
    double       solver_tolerance = 1e-10;

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
        prm.declare_entry(
          "multiplicative",
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
        omega            = prm.get_double("omega");
        shyness          = prm.get_integer("shyness");
        multiplicative   = prm.get_bool("multiplicative");
        smoother_type    = prm.get("smoother_type");
        fe_degree        = prm.get_integer("fe_degree");
        n_refinements    = prm.get_integer("n_refinements");
        max_iterations   = prm.get_integer("max_iterations");
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

    const unsigned int fe_degree;

    const Functions::ConstantFunction<dim> rhs_function;
    const Functions::ConstantFunction<dim> boundary_condition;

    // Multigrid level objects - separate triangulation and DoFHandler per level
    MGLevelObject<Triangulation<dim>> mg_triangulations;
    FE_Q<dim>                         fe_level_set;
    MGLevelObject<DoFHandler<dim>>    mg_level_set_dof_handlers;
    MGLevelObject<Vector<double>>     mg_level_sets;

    MGLevelObject<std::unique_ptr<FE_Q<dim>>> mg_fe_poisson;
    MGLevelObject<DoFHandler<dim>>            mg_dof_handlers;
    MGLevelObject<std::unique_ptr<NonMatching::MeshClassifier<dim>>>
      mg_mesh_classifiers;

    // Fine level solution and RHS
    Vector<double>    solution;
    Vector<double>    rhs;
    std::vector<bool> active_dofs;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;

    // Multigrid matrices
    MGLevelObject<SparsityPattern>      mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double>> mg_matrices;

    const MGParameters mg_params;
  };

  template <int dim>
  MGSolver<dim>::MGSolver(const MGParameters &params)
    : fe_degree(params.fe_degree)
    , rhs_function(0)
    , boundary_condition(0.0)
    , fe_level_set(params.fe_degree + 1)
    , mg_params(params)
  {}

  template <int dim>
  void
  MGSolver<dim>::make_grid()
  {
    std::cout << "Creating background meshes for MG levels" << std::endl;

    const unsigned int n_levels = mg_params.n_refinements + 1;

    // Resize all level objects
    mg_triangulations.resize(0, n_levels - 1);
    mg_level_set_dof_handlers.resize(0, n_levels - 1);
    mg_level_sets.resize(0, n_levels - 1);
    mg_fe_poisson.resize(0, n_levels - 1);
    mg_dof_handlers.resize(0, n_levels - 1);
    mg_mesh_classifiers.resize(0, n_levels - 1);

    // Create triangulations for each level (h-coarsening)
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        mg_triangulations[level].set_mesh_smoothing(
          Triangulation<dim>::limit_level_difference_at_vertices);
        GridGenerator::hyper_cube(mg_triangulations[level], -1.21, 1.21);
        mg_triangulations[level].refine_global(level);

        // Initialize FE objects for this level
        mg_fe_poisson[level] = std::make_unique<FE_Q<dim>>(fe_degree);

        std::cout << "  Level " << level << ": "
                  << mg_triangulations[level].n_active_cells() << " cells"
                  << std::endl;
      }
  }

  template <int dim>
  void
  MGSolver<dim>::setup_discrete_level_set()
  {
    std::cout << "Setting up discrete level set functions for all levels"
              << std::endl;

    const unsigned int n_levels = mg_triangulations.n_levels();
    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;

    for (unsigned int level = 0; level < n_levels; ++level)
      {
        mg_level_set_dof_handlers[level].reinit(mg_triangulations[level]);
        mg_level_set_dof_handlers[level].distribute_dofs(fe_level_set);
        mg_level_sets[level].reinit(mg_level_set_dof_handlers[level].n_dofs());

        VectorTools::interpolate(mg_level_set_dof_handlers[level],
                                 signed_distance_sphere,
                                 mg_level_sets[level]);

        // Initialize mesh classifier for this level (construct with dof handler
        // and level set vector; MeshClassifier has no default ctor)
        mg_mesh_classifiers[level] =
          std::make_unique<NonMatching::MeshClassifier<dim>>(
            mg_level_set_dof_handlers[level], mg_level_sets[level]);

        std::cout << "  Level " << level << ": "
                  << mg_level_set_dof_handlers[level].n_dofs()
                  << " level set DoFs" << std::endl;
      }
  }

  template <int dim>
  void
  MGSolver<dim>::distribute_dofs()
  {
    std::cout << "Distributing degrees of freedom for all levels" << std::endl;

    const unsigned int n_levels = mg_triangulations.n_levels();

    for (unsigned int level = 0; level < n_levels; ++level)
      {
        mg_triangulations[level].clear_user_flags();

        // Reinitialize DoFHandler for this level's triangulation
        mg_dof_handlers[level].reinit(mg_triangulations[level]);

        // Set user flags for cells inside the domain
        for (const auto &cell : mg_dof_handlers[level].active_cell_iterators())
          {
            const NonMatching::LocationToLevelSet cell_location =
              mg_mesh_classifiers[level]->location_to_level_set(cell);

            if (cell_location == NonMatching::LocationToLevelSet::inside)
              cell->set_user_flag();
          }

        mg_dof_handlers[level].distribute_dofs(*mg_fe_poisson[level]);

        std::cout << "  Level " << level << ": "
                  << mg_dof_handlers[level].n_dofs() << " DoFs" << std::endl;
      }
  }

  template <int dim>
  void
  MGSolver<dim>::initialize_matrices()
  {
    std::cout << "Initializing matrices for finest level" << std::endl;


    // Initialize fine level matrix and vectors
    for (unsigned int level = 0; level < mg_dof_handlers.n_levels(); ++level)
      {
        DynamicSparsityPattern dsp(mg_dof_handlers[level].n_dofs(),
                                   mg_dof_handlers[level].n_dofs());

        DoFTools::make_sparsity_pattern(mg_dof_handlers[level], dsp);
        sparsity_pattern.copy_from(dsp);

        stiffness_matrix.reinit(sparsity_pattern);
      }

    const unsigned int finest_level = mg_dof_handlers.n_levels() - 1;

    solution.reinit(mg_dof_handlers[finest_level].n_dofs());
    rhs.reinit(mg_dof_handlers[finest_level].n_dofs());
  }

  template <int dim>
  void
  MGSolver<dim>::assemble_system()
  {
    const unsigned int finest_level = mg_dof_handlers.n_levels() - 1;

    Step85::assemble_system(mg_dof_handlers[finest_level],
                            *mg_fe_poisson[finest_level],
                            *mg_mesh_classifiers[finest_level],
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

    const unsigned int n_levels = mg_dof_handlers.n_levels();

    // Resize MG level objects for all levels
    mg_sparsity_patterns.resize(0, n_levels - 1);
    mg_matrices.resize(0, n_levels - 1);

    std::cout << "  Number of MG levels: " << n_levels << std::endl;
    std::cout << "  Using " << mg_params.smoother_type << " smoother"
              << std::endl;
    std::cout << "  Omega: " << mg_params.omega << std::endl;
    std::cout << "  Shyness: " << mg_params.shyness << std::endl;
    std::cout << "  Mode: "
              << (mg_params.multiplicative ? "multiplicative" : "additive")
              << std::endl;

    // Initialize sparsity patterns and matrices for each level
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        DynamicSparsityPattern dsp(mg_dof_handlers[level].n_dofs(),
                                   mg_dof_handlers[level].n_dofs());
        DoFTools::make_sparsity_pattern(mg_dof_handlers[level], dsp);
        mg_sparsity_patterns[level].copy_from(dsp);
        mg_matrices[level].reinit(mg_sparsity_patterns[level]);

        std::cout << "  Level " << level << ": "
                  << mg_dof_handlers[level].n_dofs() << " x "
                  << mg_dof_handlers[level].n_dofs() << " matrix" << std::endl;
      }

    std::cout << "  MG matrices initialized for " << n_levels << " levels"
              << std::endl;
  }

  template <int dim>
  void
  MGSolver<dim>::solve()
  {
    // ADDED, not working. Make it work
    std::cout << "Solving system with multigrid" << std::endl;
    mg::Matrix<VectorType> mg_matrix(system_matrices);

    using SmootherType =
      Patch::PreconditionPatch<SparseMatrix<double>, VectorType, dim>;

    typename SmootherType::AdditionalData additional_data_example;

    MGSmootherPrecondition<SparseMatrixType, SmootherType, VectorType>
                                                         mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
      0, max_hp_level, additional_data_example);
    for (unsigned int level = 0; level <= max_hp_level; ++level)
      {
        smoother_data[level].dof_handler = &dof_handlers[level];
        smoother_data[level].relaxation  = dim == 2 ? 0.20 : 1. / 9.;
      }

    mg_smoother.initialize(system_matrices, smoother_data);
    mg_smoother.set_steps(parameters.n_smoothing_steps);

    CoarseDirectSolver coarse_direct;
    coarse_direct.initialize(system_matrices[0]);

    IterationNumberControl  coarse_control(10000, 1e-12);
    SolverGMRES<VectorType> coarse_solver(coarse_control);
    SmootherType            coarse_preconditioner;
    coarse_preconditioner.initialize(system_matrices[0], smoother_data[0]);

    MGCoarseGridIterativeSolver<VectorType,
                                SolverGMRES<VectorType>,
                                SparseMatrixType,
                                SmootherType>
      mg_coarse_solver;
    mg_coarse_solver.initialize(coarse_solver,
                                system_matrices[0],
                                coarse_preconditioner);

    AffineConstraints<double> constrains_dummy;
    constrains_dummy.close();

    using TwoLevelTransfer = MGTwoLevelTransfer<dim, VectorType>;
    // using TwoLevelTransfer = SBM::MGTransfer<dim, VectorType>;
    MGLevelObject<TwoLevelTransfer> transfers;
    AffineConstraints<double>       constraints_dummy;
    transfers.resize(0, max_hp_level);
    for (unsigned int level = 0; level < max_hp_level; ++level)
      transfers[level + 1].reinit(dof_handlers[level + 1],
                                  dof_handlers[level],
                                  constraints_dummy,
                                  constraints_dummy);


    using MGTransferType = MGTransferGlobalCoarsening<dim, VectorType>;
    MGTransferType mg_transfer(transfers, [&](const auto l, auto &vec) {
      // if (l == max_level + 1)
      //   vec.reinit(dof_handlers[l].n_dofs());
      // else
      vec.reinit(dof_handlers[l].n_dofs());
    });


    SBM::Multigrid<VectorType> mg_new(mg_matrix,
                                      coarse_direct,
                                      mg_transfer,
                                      mg_smoother,
                                      mg_smoother,
                                      0,
                                      max_hp_level,
                                      SBM::Multigrid<VectorType>::v_cycle);

    Multigrid<VectorType> mg(mg_matrix,
                             coarse_direct,
                             mg_transfer,
                             mg_smoother,
                             mg_smoother,
                             0,
                             max_hp_level,
                             Multigrid<VectorType>::v_cycle);

    MGLevelObject<SBM::PrePostTransfer<dim, VectorType>> mg_pre_post_transfer;
    mg_pre_post_transfer.resize(0, max_hp_level);
    for (unsigned int level = 0; level <= max_hp_level; ++level)
      {
        mg_pre_post_transfer[level].reinit(dof_handlers[level]);
      }

    // mg_new.set_pre_transfers(mg_pre_post_transfer);

    using PreconditionerType =
      SBM::PreconditionMG<dim, VectorType, MGTransferType>;
    PreconditionerType preconditioner(dof_handlers[max_hp_level],
                                      mg_new,
                                      mg_transfer);

    SolverControl solver_control(parameters.n_iterations, parameters.tolerance);

    SolverGMRES<VectorType>(solver_control)
      .solve(system_matrices[max_hp_level], solution, rhs, preconditioner);
    std::cout << "Solved in " << solver_control.last_step()
              << " iterations, final residual" << solver_control.last_value()
              << std::endl;
  }

  template <int dim>
  void
  MGSolver<dim>::output_results() const
  {
    std::cout << "Writing vtu file" << std::endl;

    const unsigned int finest_level = mg_dof_handlers.n_levels() - 1;

    DataOut<dim> data_out;
    data_out.add_data_vector(mg_dof_handlers[finest_level],
                             solution,
                             "solution");
    data_out.add_data_vector(mg_level_set_dof_handlers[finest_level],
                             mg_level_sets[finest_level],
                             "level_set");

    data_out.build_patches();
    std::ofstream output("mg_solver.vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  double
  MGSolver<dim>::compute_L2_error() const
  {
    std::cout << "Computing L2 error" << std::endl;

    const unsigned int finest_level = mg_dof_handlers.n_levels() - 1;

    const QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> cell_fe_values(*mg_fe_poisson[finest_level],
                                 quadrature_formula,
                                 update_values | update_JxW_values |
                                   update_quadrature_points);

    AnalyticalSolution<dim> analytical_solution;
    double                  error_L2_squared = 0;

    for (const auto &cell :
         mg_dof_handlers[finest_level].active_cell_iterators() |
           IteratorFilters::UserFlagSet())
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

    std::cout << "Classifying cells for all levels" << std::endl;
    const unsigned int n_levels = mg_triangulations.n_levels();
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        if (mg_mesh_classifiers[level])
          mg_mesh_classifiers[level]->reclassify();
      }

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

  Step85::MGParameters     params;
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
          std::cout << "Failed to parse parameter file: " << e.what()
                    << std::endl;
          return 1;
        }
    }
  else
    {
      std::cout << "Using default parameters (no parameter file provided)"
                << std::endl;
    }

  params.parse_parameters(prm);

  Step85::MGSolver<dim> mg_solver(params);
  mg_solver.run();

  return 0;
}
