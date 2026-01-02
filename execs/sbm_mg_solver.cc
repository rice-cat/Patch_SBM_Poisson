#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <fstream>
#include <iostream>
#include <memory>

#include "shy_patches.h"
#include "step-85.h"

namespace Step85
{
  using namespace dealii;

  // Proxy class to wrap block smoothers for use with distributed vectors
  // This allows RelaxationBlockSOR (which only works with Vector<double>)
  // to be used with LinearAlgebra::distributed::Vector<double>
  template <typename SmootherType>
  class BlockSmootherProxy
  {
  public:
    using VectorType     = LinearAlgebra::distributed::Vector<double>;
    using AdditionalData = typename SmootherType::AdditionalData;

    BlockSmootherProxy() = default;

    void
    initialize(const SparseMatrix<double> &matrix, const AdditionalData &data)
    {
      smoother.initialize(matrix, data);
    }

    void
    vmult(VectorType &dst, const VectorType &src) const
    {
      // Copy distributed vector to regular vector
      Vector<double> dst_serial(dst.size());
      Vector<double> src_serial(src.size());

      for (unsigned int i = 0; i < src.size(); ++i)
        src_serial(i) = src(i);

      // Apply the smoother
      smoother.vmult(dst_serial, src_serial);

      // Copy back to distributed vector
      for (unsigned int i = 0; i < dst.size(); ++i)
        dst(i) = dst_serial(i);
    }

    void
    Tvmult(VectorType &dst, const VectorType &src) const
    {
      // Copy distributed vector to regular vector
      Vector<double> dst_serial(dst.size());
      Vector<double> src_serial(src.size());

      for (unsigned int i = 0; i < src.size(); ++i)
        src_serial(i) = src(i);

      // Apply the smoother
      smoother.Tvmult(dst_serial, src_serial);

      // Copy back to distributed vector
      for (unsigned int i = 0; i < dst.size(); ++i)
        dst(i) = dst_serial(i);
    }

    void
    clear()
    {
      smoother.clear();
    }

  private:
    SmootherType smoother;
  };

  struct MGParameters
  {
    double       omega             = 1.0;
    unsigned int shyness           = 3;
    std::string  smoother_type     = "shy_patches";
    unsigned int fe_degree         = 2;
    double       cell_threshold    = 1;
    unsigned int n_refinements     = 1;
    unsigned int base_refinements  = 3;
    unsigned int max_iterations    = 100;
    unsigned int n_smoothing_steps = 2;
    double       solver_tolerance  = 1e-10;
    bool         use_p_multigrid   = false;

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
        // Note: multiplicative/additive mode removed - only multiplicative
        // smoothing is supported in this executable.
        prm.declare_entry("smoother_type",
                          "shy_patches",
                          Patterns::Selection("shy_patches|full_residual"),
                          "Type of patch smoother to use");
        prm.declare_entry("fe_degree",
                          "2",
                          Patterns::Integer(1),
                          "Finite element polynomial degree");
        prm.declare_entry(
          "cell_threshold",
          "1.0",
          Patterns::Double(0.0, 1.0),
          "Threshold for cell inclusion in domain when using level set");
        prm.declare_entry(
          "n_refinements",
          "2",
          Patterns::Integer(0),
          "Number of global refinements (levels = n_refinements + 1)");
        prm.declare_entry(
          "base_refinements",
          "0",
          Patterns::Integer(0),
          "Number of base global refinements applied to the base mesh before creating MG levels");
        prm.declare_entry("max_iterations",
                          "100",
                          Patterns::Integer(1),
                          "Maximum solver iterations");
        prm.declare_entry("n_smoothing_steps",
                          "2",
                          Patterns::Integer(1),
                          "Number of smoothing steps");
        prm.declare_entry("solver_tolerance",
                          "1e-10",
                          Patterns::Double(0.0),
                          "Solver tolerance");
        prm.declare_entry("use_p_multigrid",
                          "false",
                          Patterns::Bool(),
                          "Whether to use p-multigrid instead of h-multigrid");
      }
      prm.leave_subsection();
    }

    void
    parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Multigrid parameters");
      {
        omega             = prm.get_double("omega");
        shyness           = prm.get_integer("shyness");
        cell_threshold    = prm.get_double("cell_threshold");
        smoother_type     = prm.get("smoother_type");
        fe_degree         = prm.get_integer("fe_degree");
        n_refinements     = prm.get_integer("n_refinements");
        base_refinements  = prm.get_integer("base_refinements");
        max_iterations    = prm.get_integer("max_iterations");
        n_smoothing_steps = prm.get_integer("n_smoothing_steps");
        solver_tolerance  = prm.get_double("solver_tolerance");
        use_p_multigrid   = prm.get_bool("use_p_multigrid");
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
    solve_amg();

    void
    output_results() const;

    double
    compute_L2_error() const;

    const unsigned int fe_degree;

    const Functions::ConstantFunction<dim> rhs_function;
    const Functions::ConstantFunction<dim> boundary_condition;

    // Type aliases for MG
    using VectorType       = LinearAlgebra::distributed::Vector<double>;
    using SparseMatrixType = SparseMatrix<double>;

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
    VectorType        solution;
    VectorType        rhs;
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
    , rhs_function(1)
    , boundary_condition(0.0)
    , fe_level_set(params.fe_degree + 1)
    , mg_params(params)
  {}

  template <int dim>
  void
  MGSolver<dim>::make_grid()
  {
    std::cout << "Creating background meshes for MG levels" << std::endl;

    const unsigned int n_levels = mg_params.use_p_multigrid ?
                                    mg_params.fe_degree :
                                    mg_params.n_refinements + 1;

    // Resize all level objects
    mg_triangulations.resize(0, n_levels - 1);
    mg_level_set_dof_handlers.resize(0, n_levels - 1);
    mg_level_sets.resize(0, n_levels - 1);
    mg_fe_poisson.resize(0, n_levels - 1);
    mg_dof_handlers.resize(0, n_levels - 1);
    mg_mesh_classifiers.resize(0, n_levels - 1);

    // Create triangulations for each level
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        mg_triangulations[level].set_mesh_smoothing(
          Triangulation<dim>::limit_level_difference_at_vertices);
        GridGenerator::hyper_cube(mg_triangulations[level], -1.21, 1.21);

        if (mg_params.use_p_multigrid)
          {
            // For p-multigrid, all levels have the same mesh refinement
            mg_triangulations[level].refine_global(mg_params.base_refinements +
                                                   mg_params.n_refinements);
            // FE degree varies from 1 to fe_degree
            mg_fe_poisson[level] = std::make_unique<FE_Q<dim>>(level + 1);
          }
        else
          {
            // For h-multigrid, mesh refinement varies
            mg_triangulations[level].refine_global(mg_params.base_refinements +
                                                   level);
            // FE degree is constant
            mg_fe_poisson[level] = std::make_unique<FE_Q<dim>>(fe_degree);
          }

        std::cout << "  Level " << level << ": "
                  << mg_triangulations[level].n_active_cells() << " cells, "
                  << "FE degree " << mg_fe_poisson[level]->degree << std::endl;
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
    std::cout
      << "Distributing degrees of freedom for all levels with cell threshold "
      << mg_params.cell_threshold << std::endl;

    const unsigned int n_levels = mg_triangulations.n_levels();


    const QGauss<1>       quadrature_1D(2);
    hp::FECollection<dim> collection_dummy;
    collection_dummy.push_back(fe_level_set);

    NonMatching::RegionUpdateFlags region_update_flags;
    region_update_flags.inside = update_JxW_values;

    for (unsigned int level = 0; level < n_levels; ++level)
      {
        mg_triangulations[level].clear_user_flags();

        // Reinitialize DoFHandler for this level's triangulation
        mg_dof_handlers[level].reinit(mg_triangulations[level]);


        std::vector<types::global_dof_index> level_set_dof_indices(
          fe_level_set.dofs_per_cell);

        NonMatching::FEValues<dim> non_matching_fe_values(
          collection_dummy,
          quadrature_1D,
          region_update_flags,
          *mg_mesh_classifiers[level],
          mg_level_set_dof_handlers[level],
          mg_level_sets[level]);

        // Set user flags for cells inside the domain
        for (const auto &cell : mg_dof_handlers[level].active_cell_iterators())
          {
            const NonMatching::LocationToLevelSet cell_location =
              mg_mesh_classifiers[level]->location_to_level_set(cell);

            if (cell_location == NonMatching::LocationToLevelSet::inside)
              cell->set_user_flag();

            if (cell_location == NonMatching::LocationToLevelSet::intersected)
              {
                const typename DoFHandler<dim>::active_cell_iterator
                  lvl_set_cell(&mg_triangulations[level],
                               cell->level(),
                               cell->index(),
                               &mg_level_set_dof_handlers[level]);
                non_matching_fe_values.reinit(lvl_set_cell);
                const Quadrature<dim> inside_quad =
                  non_matching_fe_values.get_inside_fe_values()
                    ->get_quadrature();
                double inside_ratio = 0;
                for (const double q_weight : inside_quad.get_weights())
                  {
                    inside_ratio += q_weight;
                  }
                if (inside_ratio > mg_params.cell_threshold)
                  cell->set_user_flag();
              }
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
    // Note: only multiplicative smoothing is supported by this executable.

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
  MGSolver<dim>::assemble_system()
  {
    const unsigned int finest_level = mg_dof_handlers.n_levels() - 1;

    Step85::assemble_system(mg_dof_handlers[finest_level],
                            *mg_fe_poisson[finest_level],
                            *mg_mesh_classifiers[finest_level],
                            mg_fe_poisson[finest_level]->degree,
                            rhs_function,
                            boundary_condition,
                            stiffness_matrix,
                            rhs,
                            true,
                            active_dofs);

    std::vector<bool> active_dofs_on_level;
    for (unsigned int level = 0; level <= finest_level; ++level)
      {
        Step85::assemble_system(mg_dof_handlers[level],
                                *mg_fe_poisson[level],
                                *mg_mesh_classifiers[level],
                                mg_fe_poisson[level]->degree,
                                rhs_function,
                                boundary_condition,
                                mg_matrices[level],
                                rhs,
                                false,
                                active_dofs_on_level);
        if (level == 0)
          // Ensure inactive rows are well-posed by putting a unit diagonal on
          // them.
          for (types::global_dof_index d = 0;
               d < static_cast<types::global_dof_index>(
                     active_dofs_on_level.size());
               ++d)
            if (!active_dofs_on_level[d])
              mg_matrices[level].set(d, d, 1.0);
      }
  }



  class CoarseDirectSolver : public ::dealii::MGCoarseGridBase<
                               LinearAlgebra::distributed::Vector<double>>
  {
    using VectorType = LinearAlgebra::distributed::Vector<double>;

  public:
    void
    operator()(const unsigned int level,
               VectorType &       dst,
               const VectorType & src) const override
    {
      (void)level;
      Vector<double> tmp_src;
      tmp_src.reinit(src.size());
      for (unsigned int i = 0; i < src.size(); ++i)
        tmp_src(i) = src(i);
      Vector<double> tmp_dst;
      tmp_dst.reinit(dst.size());
      solver_direct.vmult(tmp_dst, tmp_src);
      for (unsigned int i = 0; i < dst.size(); ++i)
        dst(i) = tmp_dst(i);
    }

    void
    initialize(const SparseMatrix<double> &matrix)
    {
      solver_direct.initialize(matrix);
    }

    SparseDirectUMFPACK solver_direct;
  };


  template <int dim>
  void
  MGSolver<dim>::solve()
  {
    std::cout << "Solving system with multigrid" << std::endl;

    const unsigned int min_level = mg_matrices.min_level();
    const unsigned int max_level = mg_matrices.max_level();

    // Create MG matrix wrapper
    mg::Matrix<VectorType> mg_matrix(mg_matrices);

    // Setup coarse solver
    CoarseDirectSolver mg_coarse_solver;
    mg_coarse_solver.initialize(mg_matrices[min_level]);

    // Setup smoother using RelaxationBlock with proxy wrapper
    using BaseSmootherType =
      RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>;
    using SmootherType               = BlockSmootherProxy<BaseSmootherType>;
    using SmootherAdditionalDataType = BaseSmootherType::AdditionalData;

    MGSmootherPrecondition<SparseMatrixType, SmootherType, VectorType>
      mg_smoother;

    MGLevelObject<SmootherAdditionalDataType> smoother_data(min_level,
                                                            max_level);

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        // Create a predicate function that uses user flags to determine if a
        // cell is in the domain
        auto cell_is_in_domain =
          [](const typename DoFHandler<dim>::cell_iterator &cell) -> bool {
          return cell->user_flag_set();
        };

        // Create patches for each level
        if (mg_params.smoother_type == "shy_patches")
          make_shy_vertex_patches(
            smoother_data[level].block_list,
            mg_dof_handlers[level],
            numbers::invalid_unsigned_int, // level within the DoFHandler (use
                                           // numbers::invalid_unsigned_int for
                                           // active)
            cell_is_in_domain,
            mg_params.shyness,
            mg_params.n_smoothing_steps);
        else if (mg_params.smoother_type == "full_residual")
          make_full_residual_vertex_patches(
            smoother_data[level].block_list,
            mg_dof_handlers[level],
            numbers::invalid_unsigned_int, // level within the DoFHandler (use
                                           // numbers::invalid_unsigned_int for
                                           // active)
            cell_is_in_domain,
            mg_params.shyness,
            mg_params.n_smoothing_steps);
        else if (true)
          AssertThrow(false,
                      ExcMessage("Unknown smoother type specified in "
                                 "parameters."));

        smoother_data[level].relaxation = mg_params.omega;
        smoother_data[level].inversion  = PreconditionBlockBase<double>::svd;
      }

    mg_smoother.initialize(mg_matrices, smoother_data);
    mg_smoother.set_steps(1); // Number of smoothing steps

    // Measure single smoother sweep on the finest level
    {
      VectorType dst(rhs);
      VectorType src(rhs);
      Timer      timer;
      mg_smoother.smooth(max_level, dst, src);
      const double smoother_time = timer.wall_time();
      std::cout << "SUMMARY: smoother_sweep_time " << smoother_time << " s"
                << std::endl;
    }

    using TwoLevelTransfer = MGTwoLevelTransfer<dim, VectorType>;
    MGLevelObject<TwoLevelTransfer> transfers;
    AffineConstraints<double>       constraints_dummy;
    transfers.resize(0, max_level);
    for (unsigned int level = 0; level < max_level; ++level)
      transfers[level + 1].reinit(mg_dof_handlers[level + 1],
                                  mg_dof_handlers[level],
                                  constraints_dummy,
                                  constraints_dummy);


    using MGTransferType = MGTransferGlobalCoarsening<dim, VectorType>;
    MGTransferType mg_transfer(transfers, [&](const auto l, auto &vec) {
      // if (l == max_level + 1)
      //   vec.reinit(dof_handlers[l].n_dofs());
      // else
      vec.reinit(mg_dof_handlers[l].n_dofs());
    });

    // Create multigrid object
    Multigrid<VectorType> mg(mg_matrix,
                             mg_coarse_solver,
                             mg_transfer,
                             mg_smoother,
                             mg_smoother,
                             min_level,
                             max_level);

    // Create MG preconditioner
    PreconditionMG<dim, VectorType, MGTransferType> preconditioner(
      mg_dof_handlers[max_level], mg, mg_transfer);

    // Solve with GMRES
    IterationNumberControl  solver_control(mg_params.max_iterations,
                                          mg_params.solver_tolerance);
    SolverGMRES<VectorType> solver(solver_control);

    try
      {
        Timer timer;
        solver.solve(stiffness_matrix, solution, rhs, preconditioner);
        const double solve_time = timer.wall_time();
        std::cout << "SUMMARY: gmg_solve_time " << solve_time << " s"
                  << std::endl;
        std::cout << "SUMMARY: gmg_iterations " << solver_control.last_step()
                  << std::endl;

        std::cout << "  Solved in " << solver_control.last_step()
                  << " iterations, final residual "
                  << solver_control.last_value() << std::endl;
      }
    catch (std::exception &e)
      {
        std::cout << "  Solver failed: " << e.what() << std::endl;
      }
  }



  template <int dim>
  void
  MGSolver<dim>::solve_amg()
  {
    std::cout << "Solving system with AMG" << std::endl;

    Timer timer;

    TrilinosWrappers::SparseMatrix trilinos_matrix;
    trilinos_matrix.reinit(stiffness_matrix);

    IndexSet index_set(solution.size());
    index_set.add_range(0, solution.size());

    TrilinosWrappers::MPI::Vector trilinos_solution(index_set, MPI_COMM_SELF);
    TrilinosWrappers::MPI::Vector trilinos_rhs(index_set, MPI_COMM_SELF);
    for (unsigned int i = 0; i < rhs.size(); ++i)
      trilinos_rhs(i) = rhs(i);

    TrilinosWrappers::PreconditionAMG                 preconditioner;
    TrilinosWrappers::PreconditionAMG::AdditionalData data;
    data.elliptic = true;
    preconditioner.initialize(trilinos_matrix, data);

    IterationNumberControl solver_control(mg_params.max_iterations,
                                          mg_params.solver_tolerance);
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

    try
      {
        timer.restart();
        solver.solve(trilinos_matrix,
                     trilinos_solution,
                     trilinos_rhs,
                     preconditioner);
        const double solve_time = timer.wall_time();

        for (unsigned int i = 0; i < solution.size(); ++i)
          solution(i) = trilinos_solution(i);

        std::cout << "SUMMARY: amg_solve_time " << solve_time << " s"
                  << std::endl;
        std::cout << "SUMMARY: amg_iterations " << solver_control.last_step()
                  << std::endl;

        std::cout << "  Solved in " << solver_control.last_step()
                  << " iterations, final residual "
                  << solver_control.last_value() << std::endl;
      }
    catch (std::exception &e)
      {
        std::cout << "  AMG Solver failed: " << e.what() << std::endl;
      }
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
    setup_multigrid();
    assemble_system();
    solve();
    // only solve with AMG for linear elements, as a comparison
    // for higher order elements the AMG preconditioner does not converge
    // anyway.
    if (fe_degree == 1)
      solve_amg();

    output_results();

    const double error_L2 = compute_L2_error();
    std::cout << "\nL2 Error: " << error_L2 << std::endl;
  }


} // namespace Step85

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  constexpr int                            dim = DEAL_II_DIMENSION;

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
