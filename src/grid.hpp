#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>

#include "materialPoints.hpp"
#include "snowModel.hpp"
#include "collisionObject.hpp"

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
  void zeroForce();
  void addForce(Vector3f force);
  void explicitUpdateVelocity(double timestep);
  void semiImplicitUpdateVelocity(double beta);
  Vector3f getVelocity();
  Vector3f getVelocityChange();
  void detectCollision(CollisionObject *co, double timestep);

  std::vector<GridNode *> m_neighbors;
  std::vector<GridCell *> m_neighborCells;

  // Physical properties
  Vector3f m_velocity;
  double m_mass;

private:
  float cubicBSpline(float x) const;
  float gradCubicBSpline(float x) const;

  Vector3i m_idx;
  Grid *m_grid;

  // Physical properties
  Vector3f m_velocityChange;
  Vector3f m_force;
};

class Grid {
public:
  Grid(Vector3f origin, Vector3i dimensions, float spacing);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Rasterization methods

  void rasterizeParticlesToGrid();
  void setInitialVolumesAndDensities(MaterialPoints &materialPoints);
  void computeGridForces(MaterialPoints &materialPoints,
                         struct SnowModel snowModel);
  std::vector<GridNode *> getNearbyNodes(MaterialPoint *particle,
                                         double radius = 2.0);

  // Grid sizing and location

  Vector3f m_origin;
  Vector3i m_dim;
  float m_spacing; // h in the paper

private:
  // Utility methods
  inline Vector3i idxToVector(int idx);
  inline int vectorToIdx(Vector3i idx);

  // Grid data
  // the x dimension varies the fastest, followed by y and then z.

  std::vector<GridNode *> m_gridNodes;
  std::vector<GridCell *> m_gridCells;
};

} // namespace SnowSimulator

#endif // GRID_H
