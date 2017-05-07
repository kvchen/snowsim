#include <Eigen/Dense>
#include <assert.h>
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

GridCell::GridCell() {}

void GridCell::addMaterialPoint(MaterialPoint *materialPoint) {
  m_materialPoints.push_back(materialPoint);
  materialPoint->m_cell = this;
}

void GridCell::clear() { m_materialPoints.clear(); }

// ============================================================================
// GRIDNODE METHODS
// ============================================================================

GridNode::GridNode(Vector3i idx, Grid *grid)
    : m_idx(idx), m_grid(grid), m_mass(0), m_velocity(Vector3f::Zero()),
      m_nextVelocity(Vector3f::Zero()), m_force(Vector3f::Zero()) {}

void GridNode::zeroForce() { m_force = Vector3f::Zero(); }

void GridNode::addForce(Vector3f force) { m_force += force; }

Vector3f GridNode::getVelocity() { return m_velocity; }

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
  // int sign = copysign(1, x);
  int sign = x > 0 ? 1 : -1;

  if (absx < 1) {
    return sign * (1.5 * absx - 2) * absx;
  } else if (absx < 2) {
    return sign * ((-0.5f * absx + 2) * absx - 2);
  } else {
    return 0;
  }
}

/**
 * Computes the weight necessary for rasterizing the properties of a
 * materialPoint to the grid. This is referred to as w in the paper.
 */
float GridNode::basisFunction(Vector3f particlePos) const {
  Vector3f offset = particlePos / m_grid->m_spacing - m_idx.cast<float>();

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
  return Vector3f(gx, gy, gz);
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
    Vector3i maxIdx = (idx.array() + 1).min(m_dim.array() - 1);

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

std::vector<GridNode *> Grid::getAllNodes() { return m_gridNodes; }
