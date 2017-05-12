#ifndef SNOWSIM_SEMIIMPLICITMETHOD_H
#define SNOWSIM_SEMIIMPLICITMETHOD_H

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include "grid.hpp"
#include "materialPoints.hpp"

using namespace Eigen;

namespace SnowSimulator {

class SemiImplicitMethod {
public:
  SemiImplicitMethod(Grid *grid);

  // Implicit integration methods

  void solve(double timestep, double beta = 0.5,
             double residualThreshold = 1e-20, double maxIterations = 100);
  Vector3f computePotentialHessian(GridNode *node, Vector3f deltaU);
  Matrix3f computeAp(MaterialPoint *mp, Vector3f &deltaU, float mu,
                     float lambda);
  Matrix3f computeDeltaCofactorFe(Matrix3f &Fe, Matrix3f &dFe);
  Matrix3f computeDeltaRe(Matrix3f &Re, Matrix3f &Se, Matrix3f &deltaFe);

  // Utility methods

  float frobeniusInnerProduct(Matrix3f &x, Matrix3f &y) const;
  Matrix3f cofactor(Matrix3f &x) const;

private:
  Grid *m_grid;
  std::shared_ptr<spdlog::logger> logger;
};

} // namespace SnowSimulator

#endif // SNOWSIM_SEMIIMPLICITMETHOD_H
