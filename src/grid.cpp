#include <Eigen/Dense>
#include <math.h>

#include "grid.hpp"

using namespace Eigen;
using namespace SnowSimulator;

GridNode::GridNode(Vector3i idx, Grid *grid) : m_idx(idx), m_grid(grid) {
  m_surroundingCells.push_back(new GridCell(this));
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
    for (int j = 0; j < 27; j++) {
      // offset = Vector3i();
      // m_gridNodes[]
    }
  }

  // Vector3i m_bbox_min =
}

/**
 * The cubic B-spline formulation of the grid basis function used for
 * determining how particle properties are transferred to the grid. Both this
 * function and the gradient of this function below are called for each pair
 * of particle and grid node.
 *
 * TODO(kvchen): Look into parallelizing this.
 */
float Grid::basisFunction(float x) {
  float absx = fabs(x);

  if (absx >= 2) {
    return 0;
  } else if (absx < 1) {
    return (0.5f * absx - 1.f) * absx * absx + (2.0f / 3);
  } else {
    return (((-1.0f / 6) * absx + 1) * absx - 2) * absx + (4.0f / 3);
  }
}

float Grid::gradBasisFunction(float x) {
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
float Grid::transferWeight(Eigen::Vector3i gridIdx,
                           Eigen::Vector3f particlePos) {
  float invSpacing = 1.0f / m_spacing;
  Eigen::Vector3f fractionalCellOffset =
      invSpacing * particlePos - gridIdx.cast<float>();

  return basisFunction(fractionalCellOffset.x()) *
         basisFunction(fractionalCellOffset.y()) *
         basisFunction(fractionalCellOffset.z());
}

float Grid::gradTransferWeight(Eigen::Vector3i gridIdx,
                               Eigen::Vector3f particlePos) {
  return 1.0f;
}

// TODO(elbertlin168): We'll definitely want to parallelize this.
// Also might be better to associate mpoints with nodes rather than cells,
// since we need to calculate the gradient of w_{ip} between each mpoint/node 
// pair (unless I'm doing something stupid looping through the points).
// Instead might try looping through every particle, then transfering to
// nearby nodes
// Need to define params LAMBDA0, HARDENING
// LAMBDA0 = Young's * Poisson's / (1 + Poisson's) / (1 - 2 * Poisson's)
void Grid::computeForces() {
  for (std::vector<GridNode *>::iterator i1 = m_gridNodes.begin(); i1 != m_gridNodes.end(); ++i1) {
    GridNode *curNode = *i1;
    for (std::vector<GridCell *>::iterator i2 = curNode->m_surroundingCells.begin(); i2 != curNode->m_surroundingCells.end(); ++i2) {
      GridCell *curCell = *i2;
      for (std::vector<MaterialPoint *>::iterator i3 = curCell->m_materialPoints.begin(); i3 != curCell->m_materialPoints.end(); ++i3) {
        MaterialPoint *curPoint = *i3;
        float Jp = curPoint->defPlastic.determinant();
        float lambda = LAMBDA0 * exp(HARDENING * (1 - Jp));
      }
    }
  }
}
