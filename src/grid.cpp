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

// GridCell::GridCell(GridNode *node) : m_node(node) {}
GridCell::GridCell() {}

void GridCell::addMaterialPoint(MaterialPoint *materialPoint) {
  m_materialPoints.push_back(materialPoint);
  materialPoint->m_cell = this;
}

void GridCell::clear() {
  // for (auto &mp : m_materialPoints) {
  //   mp->m_cell = nullptr;
  // }
  m_materialPoints.clear();
}

// ============================================================================
// GRIDNODE METHODS
// ============================================================================

GridNode::GridNode(Vector3i idx, Grid *grid)
    : m_idx(idx), m_grid(grid), m_mass(0), m_velocity(Vector3f::Zero()) {
  // m_surroundingCells.push_back(new GridCell(this));
}

/**
 * TODO(kvchen): Look into parallelizing this with OpenMP.
 */
void GridNode::rasterizeMaterialPoints() {
  // Rasterize mass to the grid
  for (auto const &cell : m_neighborCells) {
    for (auto const &mp : cell->m_materialPoints) {
      m_mass += mp->m_mass * basisFunction(mp->m_position);
    }
  }

  // Rasterize velocity to the grid
  if (m_mass == 0) {
    m_velocity.setZero();
  } else {
    for (auto const &cell : m_neighborCells) {
      for (auto const &mp : cell->m_materialPoints) {
        m_velocity +=
            mp->m_velocity * mp->m_mass * basisFunction(mp->m_position);
      }
    }

    m_velocity /= m_mass;
  }
}

/**
 * The cubic B-spline formulation of the grid basis function used for
 * determining how particle properties are transferred to the grid. Both
 * this function and the gradient of this function below are called for each
 * pair of particle and grid node. This is referred to as N in the paper.
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

  logger->info("Creating grid of size ({}, {}, {})", m_dim.x(), m_dim.y(),
               m_dim.z());

  int num_nodes = m_dim.prod();
  logger->info("Total node count: {}", num_nodes);

  // Allocate space for each of the gridNodes

  logger->info("Allocating grid nodes and cells...");

  for (int i = 0; i < num_nodes; i++) {
    GridNode *node = new GridNode(idxToVector(i), this);
    m_gridNodes.push_back(node);
    m_gridCells.push_back(new GridCell());
  }

  // Make sure each GridNode has a reference to the 4^3 surrounding GridCells
  // within the support radius of the basis function.
  // Each GridCell also needs to know about the 6^3 surrounding GridNodes.

  logger->info("Linking grid nodes to neighboring cells...");

  for (int i = 0; i < num_nodes; i++) {
    Vector3i idx = idxToVector(i);
    GridNode *current = m_gridNodes[i];

    // GridCell indexing is such that the corresponding GridNode is at the
    // top-left corner of the GridCell.

    Vector3i minIdx = (idx.array() - 2).max(0);
    Vector3i maxIdx = (idx.array() + 2).min(m_dim.array());

    // Vector3i offsetDim = maxIdx - maxIdx;
    // Vector3i numNeighboringCells = offsetDim.prod();

    for (int ox = minIdx.x(); ox < maxIdx.x(); ox++) {
      for (int oy = minIdx.y(); oy < maxIdx.y(); oy++) {
        for (int oz = minIdx.z(); oz < maxIdx.z(); oz++) {
          int cellIdx = vectorToIdx(Vector3i(ox, oy, oz));
          GridCell *neighborCell = m_gridCells[cellIdx];
          current->m_neighborCells.push_back(neighborCell);
        }
      }
    }
  }
}

/**
 * Uses the grid basis function to computes the weighting factor for each pair
 * of particle and grid cell.
 */
// void Grid::computeWeights() {}

// void Grid::rasterizeMaterialPoints(MaterialPoints &materialPoints) {
//   // materialPoints
//   // Move each material point into its respective GridCell.
// }

// inline GridNode *Grid::getNode(Vector3i idx) {
//   return m_gridNodes[vectorToIdx(idx)];
// }

void Grid::setInitialVolumesAndDensities(MaterialPoints &materialPoints) {
  for (auto &mp : materialPoints.m_materialPoints) {
    for (auto &node : getNearbyNodes(mp, 2.0)) {
      mp->m_density += node->m_mass * node->basisFunction(mp->m_position);
    }

    if (mp->m_density != 0) {
      mp->m_volume = mp->m_mass / mp->m_density;
    }
  }
}

std::vector<GridNode *> Grid::getNearbyNodes(MaterialPoint *particle,
                                             double radius) {
  Vector3f min =
      (particle->m_position.array() / m_spacing - radius).ceil().max(0);

  Vector3f max = (particle->m_position.array() / m_spacing + radius)
                     .floor()
                     .max(m_dim.array().cast<float>());
  std::vector<GridNode *> nearbyNodes;

  for (int x = min.x(); x < max.x(); x++) {
    for (int y = min.y(); y < max.y(); y++) {
      for (int z = min.z(); z < max.z(); z++) {
        int idx = vectorToIdx(Vector3i(x, y, z));
        nearbyNodes.push_back(m_gridNodes[idx]);
      }
    }
  }

  return nearbyNodes;
}

inline Vector3i Grid::idxToVector(int idx) {
  return Vector3i(idx % m_dim.x(), (idx / m_dim.x()) % m_dim.y(),
                  idx / (m_dim.x() * m_dim.y()));
}

inline int Grid::vectorToIdx(Vector3i idx) {
  return idx.x() + m_dim.x() * (idx.y() + m_dim.y() * idx.z());
}
