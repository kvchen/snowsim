#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>

#include "materialPoints.hpp"

using namespace Eigen;

namespace SnowSimulator {

// Forward declarations

class GridCell;
class GridNode;
class Grid;

/**
 * A simple class that binds a collection of particles to a particular
 * GridNode.
 */
class GridCell {
public:
  GridCell(GridNode *gridNode) : m_gridNode(gridNode) {}

  GridNode *m_gridNode;
  std::vector<MaterialPoint *> m_materialPoints;
};

class GridNode {
public:
  GridNode(Vector3i idx, Grid *grid) : m_idx(idx), m_grid(grid) {}
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Vector3i m_idx;
  Grid *m_grid;
};

class Grid {
public:
  Grid(MaterialPoints &materialPoints, float spacing)
      : m_materialPoints(materialPoints), m_spacing(spacing) {}
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Rasterization methods

  void rasterizeParticlesToGrid();
  void computeParticleVolumesAndDensities();

  // Helper functions

  float basisFunction(float x);
  float gradBasisFunction(float x);
  float transferWeight(Vector3i gridIdx, Vector3f particlePos);
  float gradTransferWeight(Vector3i gridIdx, Vector3f particlePos);

  // int index(int i, int j, int k);

private:
  MaterialPoints &m_materialPoints;

  // Grid sizing
  Vector3i m_bbox_min;
  Vector3i m_bbox_max;

  float m_spacing; // h in the paper

  // Eigen::VectorXf m_masses;
  // Eigen::VectorXf m_velocities;
  // Eigen::VectorXf m_prevVelocities;
};

} // namespace SnowSimulator

#endif // GRID_H
