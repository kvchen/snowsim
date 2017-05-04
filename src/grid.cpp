#include <Eigen/Dense>
#include <math.h>

#include "grid.hpp"

using namespace SnowSimulator;

// int ParticleGrid::index(int i, int j, int k) {
//   return i * (m_dim.y() * m_dim.z()) + j * (m_dim.z()) + k;
// }

/**
 * The cubic B-spline formulation of the grid basis function used for
 * determining how particle properties are transferred to the grid. Both this
 * function and the gradient of this function below are called for each pair
 * of particle and grid node.
 *
 * TODO(kvchen): Look into parallelizing this.
 */
float ParticleGrid::basisFunction(float x) {
  float absx = fabs(x);

  if (absx >= 2) {
    return 0;
  } else if (absx < 1) {
    return (0.5f * absx - 1.f) * absx * absx + (2.0f / 3);
  } else {
    return (((-1.0f / 6) * absx + 1) * absx - 2) * absx + (4.0f / 3);
  }
}

float ParticleGrid::gradBasisFunction(float x) {
  if (x < 0) {
    // TODO(kvchen): This recursive call is probably inefficient. Should we
    // bother listing out the additional cases? Benchmark after confirmed
    // working.
    return -gradBasisFunction(-x);
  } else if (x < 1) {
    return (1.5f * x - 2) * x;
  } else if (x < 2) {
    return (-0.5f * x + 2) * x - 2;
  } else {
    return 0;
  }
}

/**
 * Computes the weight w_{ip} used to rasterize a particle to the grid. This is
 * used in steps 1, 2, 7, and 8 of the MPM procedure.
 */
float ParticleGrid::transferWeight(Eigen::Vector3i gridIdx,
                                   Eigen::Vector3f particlePos) {
  float invSpacing = 1.0f / m_spacing;
  Eigen::Vector3f fractionalCellOffset =
      invSpacing * particlePos - gridIdx.cast<float>();

  return basisFunction(fractionalCellOffset.x()) *
         basisFunction(fractionalCellOffset.y()) *
         basisFunction(fractionalCellOffset.z());
}

float ParticleGrid::gradTransferWeight(Eigen::Vector3i gridIdx,
                                       Eigen::Vector3f particlePos) {
  return 1.0f;
}
