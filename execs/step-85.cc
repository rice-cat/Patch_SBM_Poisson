#include "step-85.h"

int
main()
{
  const int        dim = 2;
  Step85::Settings settings;
  settings.fe_degree       = 2;
  settings.smoothing_steps = 5;
  Step85::LaplaceSolver<dim> laplace_solver(settings);
  laplace_solver.run();
}
