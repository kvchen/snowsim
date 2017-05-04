#include <Eigen/Dense>
#include <boost/math/special_functions/sign.hpp>
#include <math.h>

#include "grid.hpp"

using namespace Eigen;
using namespace SnowSimulator;

GridCell::GridCell(GridNode *gridNode) : m_gridNode(gridNode) {}

GridNode::GridNode(Vector3i idx, Grid *grid) : m_idx(idx), m_grid(grid) {
  m_surroundingCells.push_back(new GridCell(this));
}

/**
 * The cubic B-spline formulation of the grid basis function used for
 * determining how particle properties are transferred to the grid. Both this
 * function and the gradient of this function below are called for each pair
 * of particle and grid node. This is referred to as N in the paper.
 *
 * TODO(kvchen): Look into parallelizing this.
 */
inline float GridNode::basisFunction(float x) const {
  float absx = fabs(x);

  if (absx >= 2) {
    return 0;
  } else if (absx < 1) {
    return (0.5f * absx - 1.f) * absx * absx + (2.0f / 3);
  } else {
    return (((-1.0f / 6) * absx + 1) * absx - 2) * absx + (4.0f / 3);
  }
}

inline float GridNode::gradBasisFunction(float x) const {
  float absx = fabs(x);
  int sign = copysign(1, x);

  if (absx < 1) {
    return sign * (1.5 * absx - 2) * absx;
  } else if (absx < 2) {
    return sign * ((-0.5f * x + 2) * x - 2);
  } else {
    return 0;
  }
}

inline Vector3i Grid::idxToVector(int idx) const {
  return Vector3i(idx % m_dim.x(), (idx / m_dim.x()) % m_dim.y(),
                  idx / (m_dim.x() * m_dim.y()));
}

Grid::Grid(Vector3f origin, Vector3i dimensions, float spacing)
    : m_origin(origin), m_dim(dimensions), m_spacing(spacing) {

  // Allocate space for each of the gridNodes

  int num_nodes = m_dim.x() * m_dim.y() * m_dim.z();
  for (int i = 0; i < num_nodes; i++) {
    m_gridNodes.push_back(new GridNode(idxToVector(i), this));
  }

  //

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < 64; j++) {
      Vector3i offset(j % 4, (j / 4) % 4, j / 16);
      // idxToVector(i) + offset
    }
  }

  // Vector3i m_bbox_min =
}

/**
 * Computes the weight w_{ip} used to rasterize a particle to the grid. This is
 * used in steps 1, 2, 7, and 8 of the MPM procedure.
 */
float Grid::transferWeight(Eigen::Vector3i gridIdx,
                           Eigen::Vector3f particlePos) {
  float invSpacing = 1.0f / m_spacing;
  Eigen::Vector3f fractionalCellOffset =
      invSpacing * particlePos - gridIdx.cast<float>();

  return 1.0f;
  // return basisFunction(fractionalCellOffset.x()) *
  //        basisFunction(fractionalCellOffset.y()) *
  //        basisFunction(fractionalCellOffset.z());
}

float Grid::gradTransferWeight(Eigen::Vector3i gridIdx,
                               Eigen::Vector3f particlePos) {
  return 1.0f;
}
