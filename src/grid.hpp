#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>

#include "materialPoints.hpp"

using namespace Eigen;

namespace SnowSimulator {

/**
 * A simple class that binds a collection of particles to a particular
 * GridNode.
 */
class GridCell {
public:
  GridCell(GridNode *gridNode) : m_gridNode(gridNode) {}

  GridNode *m_gridNode;
  std::vector<Particle *> m_particles;
}

class GridNode {
public:
  GridNode(Vector3i idx, Grid *grid) : m_idx(idx), m_grid(grid) {}

private:
  Vector3i m_idx;
  Grid *m_grid;
}

class Grid {
public:
  Grid(MaterialPointData &materialPoints, float spacing) : m_spacing(spacing) {}

  // Rasterization methods

  void rasterizeParticlesToGrid();
  void computeParticleVolumesAndDensities();

  // Helper functions

  float basisFunction(float x);
  float gradBasisFunction(float x);
  float transferWeight(Eigen::Vector3i gridIdx, Eigen::Vector3f particlePos);
  float gradTransferWeight(Eigen::Vector3i gridIdx,
                           Eigen::Vector3f particlePos);

  // int index(int i, int j, int k);

private:
  // Grid sizing

  Eigen::Vector3i m_bbox_min;
  Eigen::Vector3i m_bbox_max;
  float m_spacing; // h in the paper

  // Eigen::VectorXf m_masses;
  // Eigen::VectorXf m_velocities;
  // Eigen::VectorXf m_prevVelocities;
};

} // namespace SnowSimulator

#endif // GRID_H
