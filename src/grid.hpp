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
  GridCell(GridNode *gridNode);

  GridNode *m_gridNode;
  std::vector<MaterialPoint *> m_materialPoints;
};

class GridNode {
public:
  GridNode(Vector3i idx, Grid *grid);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  std::vector<GridCell *> m_surroundingCells;

  inline float basisFunction(float x) const;
  inline float gradBasisFunction(float x) const;

private:
  Vector3i m_idx;
  Grid *m_grid;
};

class Grid {
public:
  Grid(Vector3f origin, Vector3i dimensions, float spacing);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Rasterization methods

  void rasterizeParticlesToGrid();
  void computeParticleVolumesAndDensities();

  // Helper functions

  float transferWeight(Vector3i gridIdx, Vector3f particlePos);
  float gradTransferWeight(Vector3i gridIdx, Vector3f particlePos);

  // int index(int i, int j, int k);

private:
  inline Vector3i idxToVector(int idx) const;

  // MaterialPoints &m_materialPoints;

  // Grid sizing and location

  Vector3f m_origin;
  Vector3i m_dim;
  float m_spacing; // h in the paper

  // Grid data

  std::vector<GridNode *> m_gridNodes;
};

} // namespace SnowSimulator

#endif // GRID_H
