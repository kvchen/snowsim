#include <Eigen/Dense>
#include <boost/math/special_functions/sign.hpp>
#include <iostream>
#include <math.h>

#include "grid.hpp"
#include "spdlog/spdlog.h"

using namespace Eigen;
using namespace SnowSimulator;

// ============================================================================
// GRIDCELL METHODS
// ============================================================================

// GridCell::GridCell(GridNode *gridNode) : m_gridNode(gridNode) {}
GridCell::GridCell() {}

// ============================================================================
// GRIDNODE METHODS
// ============================================================================

GridNode::GridNode(Vector3i idx, Grid *grid) : m_idx(idx), m_grid(grid) {
  // m_surroundingCells.push_back(new GridCell(this));
}

/**
 * The cubic B-spline formulation of the grid basis function used for
 * determining how particle properties are transferred to the grid. Both this
 * function and the gradient of this function below are called for each pair
 * of particle and grid node. This is referred to as N in the paper.
 *
 * TODO(kvchen): Look into parallelizing this.
 */
float GridNode::cubicBSpline(float x) const {
  float absx = fabs(x);

  if (absx >= 2) {
    return 0;
  } else if (absx < 1) {
    return (0.5f * absx - 1.f) * absx * absx + (2.0f / 3);
  } else {
    return (((-1.0f / 6) * absx + 1) * absx - 2) * absx + (4.0f / 3);
  }
}

float GridNode::gradCubicBSpline(float x) const {
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

/**
 * Computes the weight necessary for rasterizing the properties of a
 * materialPoint to the grid. This is referred to as w in the paper.
 */
float GridNode::basisFunction(Vector3f particlePos) const {
  float invSpacing = 1.0f / m_grid->m_spacing;
  Vector3f offset = invSpacing * particlePos - m_idx.cast<float>();

  return cubicBSpline(offset.x()) * cubicBSpline(offset.y()) *
         cubicBSpline(offset.z());
}

Vector3f GridNode::gradBasisFunction(Vector3f particlePos) const {
  float invSpacing = 1.0f / m_grid->m_spacing;
  Vector3f offset = invSpacing * particlePos - m_idx.cast<float>();

  float gx = gradCubicBSpline(offset.x()) * cubicBSpline(offset.y()) *
             cubicBSpline(offset.z());
  float gy = cubicBSpline(offset.x()) * gradCubicBSpline(offset.y()) *
             cubicBSpline(offset.z());
  float gz = cubicBSpline(offset.x()) * cubicBSpline(offset.y()) *
             gradCubicBSpline(offset.z());
  return invSpacing * Vector3f(gx, gy, gz);
}

// ============================================================================
// GRID METHODS
// ============================================================================

Grid::Grid(Vector3f origin, Vector3i dimensions, float spacing)
    : m_origin(origin), m_dim(dimensions), m_spacing(spacing) {
  auto logger = spdlog::get("snowsim");
  int num_nodes = m_dim.x() * m_dim.y() * m_dim.z();

  logger->info("Creating grid of size ({}, {}, {})", m_dim.x(), m_dim.y(),
               m_dim.z());
  logger->info("Total node count: {}", num_nodes);

  // Allocate space for each of the gridNodes

  logger->info("Allocating grid nodes and cells...");

  for (int i = 0; i < num_nodes; i++) {
    m_gridNodes.push_back(new GridNode(idxToVector(i), this));
    m_gridCells.push_back(new GridCell());
  }

  // Make sure each GridNode has a reference to the surrounding GridCells
  // within the support radius of the basis function.

  logger->info("Linking grid nodes to neighboring cells...");

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < 64; j++) {
      Vector3i offset(j % 4 - 2, (j / 4) % 4 - 2, j / 16 - 2);
      Vector3i offsetIdx = idxToVector(i) + offset;

      if ((0 <= offsetIdx.x() && offsetIdx.x() < m_dim.x()) &&
          (0 <= offsetIdx.y() && offsetIdx.y() < m_dim.y()) &&
          (0 <= offsetIdx.z() && offsetIdx.z() < m_dim.z())) {

        // std::cout << "------------" << std::endl;
        // std::cout << offsetIdx << std::endl;
        // std::cout << idxToVector(i) << std::endl;

        GridNode *current = m_gridNodes[i];
        GridCell *neighbor = m_gridCells[vectorToIdx(offsetIdx)];
        current->m_neighbors.push_back(neighbor);
      }
    }
  }
}

// void Grid::rasterizeMaterialPoints(MaterialPoints &materialPoints) {
//   // materialPoints
//   // Move each material point into its respective GridCell.
// }

inline Vector3i Grid::idxToVector(int idx) {
  return Vector3i(idx % m_dim.x(), (idx / m_dim.x()) % m_dim.y(),
                  idx / (m_dim.x() * m_dim.y()));
}

// TODO(kvchen): FIX THIS THIS IS INCORRECT
inline int Grid::vectorToIdx(Vector3i idx) {
  return idx.x() + m_dim.x() * (idx.y() + m_dim.y() * idx.z());
}
