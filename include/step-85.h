#ifndef STEP_85_H
#define STEP_85_H

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/relaxation_block.h>
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

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace Step85
{
  using namespace dealii;

  struct Settings
  {
    enum DoFRenumberingStrategy
    {
      none,
      downstream,
      upstream,
      random
    };

    double                 epsilon;
    unsigned int           fe_degree;
    std::string            smoother_type;
    unsigned int           smoothing_steps;
    DoFRenumberingStrategy dof_renumbering;
    bool                   with_streamline_diffusion;
    bool                   output;
  };

  template <int dim>
  class LaplaceSolver
  {
  public:
    LaplaceSolver(const Settings &settings);

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
    write_dof_locations(const std::string &filename);

    void
    initialize_matrices();

    void
    assemble_system();

    void
    setup_smoother();

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

    hp::MappingCollection<dim> mapping_collection;

    NonMatching::MeshClassifier<dim> mesh_classifier;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    Vector<double>       rhs;

    using SmootherType =
      RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>;
    using SmootherAdditionalDataType = SmootherType::AdditionalData;

    std::unique_ptr<SmootherType> patch_smoother;
    SmootherAdditionalDataType    smoother_data;

    const Settings settings;
  };

  template <int dim>
  class AnalyticalSolution : public Function<dim>
  {
  public:
    double
    value(const Point<dim>  &point,
          const unsigned int component = 0) const override;
  };

  // Standalone assembly function for reuse in smoother testing and MG solver
  template <int dim, typename VectorType = Vector<double>>
  void
  assemble_system(const DoFHandler<dim>                     &dof_handler,
                  const FE_Q<dim>                           &fe_poisson,
                  const NonMatching::MeshClassifier<dim>    &mesh_classifier,
                  const unsigned int                         fe_degree,
                  const Functions::ConstantFunction<dim>    &rhs_function,
                  const Functions::ConstantFunction<dim>    &boundary_condition,
                  SparseMatrix<double>                      &stiffness_matrix,
                  VectorType                                &rhs,
                  std::vector<bool>                         &active_dofs);

} // namespace Step85

#endif // STEP_85_H
