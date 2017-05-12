#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include "grid.hpp"
#include "materialPoints.hpp"
#include "semiImplicitMethod.hpp"

using namespace Eigen;
using namespace SnowSimulator;

typedef Array<float, 1, Dynamic> RowArrayXf;

SemiImplicitMethod::SemiImplicitMethod(Grid *grid) : m_grid(grid) {
  logger = spdlog::get("snowsim");

  int numNodes = m_grid->nodes().size();

  r = Matrix3Xd::Zero(3, numNodes);
  s = Matrix3Xd::Zero(3, numNodes);
}

/**
 * Uses the conjugate residual method to solve for the new velocities on each
 * GridNode.
 */
void SemiImplicitMethod::solve(double timestep, double beta,
                               double residualThreshold, double maxIterations) {
  auto nodes = m_grid->nodes();
  int numNodes = nodes.size();

  // Initialize our linear system

  r.setZero();
  s.setZero();

  for (int i = 0; i < numNodes; i++) {
    auto node = nodes[i];

    if (node->mass() > 0) {
      double scale = beta * timestep * timestep / node->mass();

      r.col(i) = scale * computePotentialHessian(node, node->velocityStar());
      s.col(i) = r.col(i) + scale * computePotentialHessian(node, r.col(i));
    }
  }

  Matrix3Xf p(r);
  Matrix3Xf q(s);

  RowVectorXf gamma = r.cwiseProduct(s).colwise().sum();

  double residual = residualThreshold + 1;
  int iteration = 0;

  while (iteration < maxIterations && residual > residualThreshold) {
    RowVectorXf alpha =
        gamma.array() / (q.colwise().squaredNorm().array() + 0.000001);
    Matrix3Xf residuals = p.array().rowwise() * alpha.array();

    r -= residuals;

    for (int i = 0; i < numNodes; i++) {
      auto node = nodes[i];

      // alpha(i) = gamma(i) /

      node->nextVelocity() += residuals.col(i);

      if (node->mass() > 0) {
        double scale = beta * timestep * timestep / node->mass();
        s.col(i) = r.col(i) + scale * computePotentialHessian(node, r.col(i));
      }
    }

    RowVectorXf b =
        r.cwiseProduct(s).colwise().sum().array() / (gamma.array() + 0.000001);
    gamma = b.cwiseProduct(gamma);

    p.array().rowwise() *= b.array();
    q.array().rowwise() *= b.array();

    iteration++;
    residual = residuals.colwise().squaredNorm().sum();
  }

  logger->info("Conjugate residual method completed in {} iteration(s) with "
               "final residual {}",
               iteration + 1, residual);
}

Vector3f SemiImplicitMethod::computePotentialHessian(GridNode *node,
                                                     Vector3f deltaU) {
  Vector3f hessian = Vector3f::Zero();

  for (auto &cell : node->surroundingCells()) {
    for (auto &mp : cell->particles()) {
      hessian -= mp->volume() * computeAp(mp, deltaU, mp->mu(), mp->lambda()) *
                 mp->elasticDeformation().transpose() *
                 node->gradBasisFunction(mp->position());
    }
  }

  return hessian;
}

Matrix3f SemiImplicitMethod::computeAp(MaterialPoint *mp, Vector3f &deltaU,
                                       float mu, float lambda) {
  Matrix3f Fe = mp->elasticDeformation();
  float Je = mp->Je();

  Matrix3f deltaFe = Matrix3f::Zero();

  // Compute deltaFe

  m_grid->forEachNeighbor(mp, [&mp, &deltaFe, deltaU, Fe](GridNode *node) {
    deltaFe +=
        deltaU * node->gradBasisFunction(mp->position()).transpose() * Fe;
  });

  // double Je = mp->.determinant();  // this is already cached

  // This can be rewritten as Je(Fe)^{-T} but taking the cofactor is faster...?

  Matrix3f cofactorFe = cofactor(Fe);
  float innerProductFe = frobeniusInnerProduct(cofactorFe, deltaFe);

  Matrix3f deltaCofactorFe = computeDeltaCofactorFe(Fe, deltaFe);
  Matrix3f deltaRe = computeDeltaRe(mp->Re(), mp->Se(), deltaFe);

  return 2 * mu * (deltaFe - deltaRe) + lambda * cofactorFe * innerProductFe +
         lambda * (Je - 1) * deltaCofactorFe;
}

// Utility methods

/**
 * Computes the d(J_{E_p}F^{-T}_{E_p}) term. Performs the double inner loop.
 */
Matrix3f SemiImplicitMethod::computeDeltaCofactorFe(Matrix3f &Fe,
                                                    Matrix3f &dFe) {
  Matrix3f deltaCofactorFe;

  // Outer loop

  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     deltaCofactor(i, j) = FeInvT.cwiseProduct(FeInvT).cwiseProduct(dFe);
  //   }
  // }

  // Unpacking this is faster.

  deltaCofactorFe(0, 0) = Fe(1, 1) * dFe(2, 2) + Fe(2, 2) * dFe(1, 1) -
                          Fe(1, 2) * dFe(2, 1) - Fe(2, 1) * dFe(1, 2);
  deltaCofactorFe(0, 1) = Fe(1, 2) * dFe(2, 0) + Fe(2, 0) * dFe(1, 2) -
                          Fe(1, 0) * dFe(2, 2) - Fe(2, 2) * dFe(1, 0);
  deltaCofactorFe(0, 2) = Fe(1, 0) * dFe(2, 1) + Fe(2, 1) * dFe(1, 0) -
                          Fe(1, 1) * dFe(2, 0) - Fe(2, 0) * dFe(1, 1);
  deltaCofactorFe(1, 0) = Fe(0, 2) * dFe(2, 1) + Fe(2, 1) * dFe(0, 2) -
                          Fe(0, 1) * dFe(2, 2) - Fe(2, 2) * dFe(0, 1);
  deltaCofactorFe(1, 1) = Fe(0, 0) * dFe(2, 2) + Fe(2, 2) * dFe(0, 0) -
                          Fe(0, 2) * dFe(2, 0) - Fe(2, 0) * dFe(0, 2);
  deltaCofactorFe(1, 2) = Fe(0, 1) * dFe(2, 0) + Fe(2, 0) * dFe(0, 1) -
                          Fe(0, 0) * dFe(2, 1) - Fe(2, 1) * dFe(0, 0);
  deltaCofactorFe(2, 0) = Fe(0, 1) * dFe(1, 2) + Fe(1, 2) * dFe(0, 1) -
                          Fe(0, 2) * dFe(1, 1) - Fe(1, 1) * dFe(0, 2);
  deltaCofactorFe(2, 1) = Fe(0, 2) * dFe(1, 0) + Fe(1, 0) * dFe(0, 2) -
                          Fe(0, 0) * dFe(1, 2) - Fe(1, 2) * dFe(0, 0);
  deltaCofactorFe(2, 2) = Fe(0, 0) * dFe(1, 1) + Fe(1, 1) * dFe(0, 0) -
                          Fe(0, 1) * dFe(1, 0) - Fe(1, 0) * dFe(0, 1);

  return deltaCofactorFe;
}

Matrix3f SemiImplicitMethod::computeDeltaRe(Matrix3f &Re, Matrix3f &Se,
                                            Matrix3f &deltaFe) {
  Matrix3f V = Re.transpose() * deltaFe - deltaFe.transpose() * Re;

  // Set up the linear equation Ax=b and solve it as described in the paper

  Matrix3f A;
  A << Se(0, 0) + Se(1, 1), Se(2, 1), -Se(0, 2), Se(1, 2), Se(0, 0) + Se(2, 2),
      Se(0, 1), -Se(0, 2), Se(1, 0), Se(1, 1) + Se(2, 2);

  Vector3f b;
  b << V(0, 1), V(0, 2), V(1, 2);

  Vector3f x = A.inverse() * b;
  Matrix3f U;

  U << 0, x(0), x(1), -x(0), 0, x(2), -x(1), -x(2), 0;

  return Re * U;
}

/**
 * Computes the Frobenius inner product of two matrices as described in
 * the paper (a = B : C). This is NOT the extended notation to deal with
 * nested matrices!
 */
float SemiImplicitMethod::frobeniusInnerProduct(Matrix3f &x,
                                                Matrix3f &y) const {
  return x.cwiseProduct(y).sum();
}

/**
 * Computes the cofactor of a given 3x3 matrix.
 */
Matrix3f SemiImplicitMethod::cofactor(Matrix3f &x) const {
  Matrix3f cofactor;

  cofactor(0, 0) = x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1);
  cofactor(0, 1) = x(1, 0) * x(2, 2) - x(1, 2) * x(2, 0);
  cofactor(0, 2) = x(1, 0) * x(2, 2) - x(1, 2) * x(2, 0);

  cofactor(1, 0) = x(0, 1) * x(2, 2) - x(0, 2) * x(2, 1);
  cofactor(1, 1) = x(0, 0) * x(2, 2) - x(0, 2) * x(2, 0);
  cofactor(1, 2) = x(0, 0) * x(2, 2) - x(1, 2) * x(2, 0);

  cofactor(2, 0) = x(0, 1) * x(1, 2) - x(0, 2) * x(1, 1);
  cofactor(2, 1) = x(0, 0) * x(1, 2) - x(0, 2) * x(1, 0);
  cofactor(2, 2) = x(0, 0) * x(1, 1) - x(0, 1) * x(1, 0);

  return cofactor;
}
