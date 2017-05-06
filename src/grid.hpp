#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>

#include "materialPoints.hpp"

using namespace Eigen;

namespace SnowSimulator {

// Forward declarations

class MaterialPoint;
class MaterialPoints;

class GridCell;
class GridNode;
class Grid;

/**
 * A simple class that binds a collection of particles to a particular
 * GridNode.
 */
class GridCell {
public:
  GridCell();
  // GridCell(GridNode *gridNode);
  void addMaterialPoint(MaterialPoint *materialPoint);
  void clear();

  GridNode *m_node; // Pointer to the immediately adjacent GridNode
  std::vector<MaterialPoint *> m_materialPoints;
};

class GridNode {
public:
  GridNode(Vector3i idx, Grid *grid);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  void rasterizeMaterialPoints();

  float basisFunction(Vector3f particlePos) const;
  Vector3f gradBasisFunction(Vector3f particlePos) const;

  std::vector<GridNode *> m_neighbors;
  std::vector<GridCell *> m_neighborCells;

private:
  float cubicBSpline(float x) const;
  float gradCubicBSpline(float x) const;

  Vector3i m_idx;
  Grid *m_grid;

  // Physical properties
  double m_mass;
  Vector3f m_velocity;
};

class Grid {
public:
  Grid(Vector3f origin, Vector3i dimensions, float spacing);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Rasterization methods

  void rasterizeParticlesToGrid();
  void computeParticleVolumesAndDensities(MaterialPoints &materialPoints);
  std::vector<GridNode *> getNearbyNodes(MaterialPoint *particle,
                                         double radius = 2.0);

  float m_spacing; // h in the paper

private:
  // Utility methods
  inline Vector3i idxToVector(int idx);
  inline int vectorToIdx(Vector3i idx);

  // MaterialPoints &m_materialPoints;

  // Grid sizing and location

  Vector3f m_origin;
  Vector3i m_dim;

  // Grid data
  // the x dimension varies the fastest, followed by y and then z.

  std::vector<GridNode *> m_gridNodes;
  std::vector<GridCell *> m_gridCells;
};

} // namespace SnowSimulator

#endif // GRID_H
