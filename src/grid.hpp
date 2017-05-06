#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>

#include "collisionObject.hpp"
#include "materialPoints.hpp"
#include "snowModel.hpp"

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
  void explicitUpdateVelocity(double delta_t);
  void semiImplicitUpdateVelocity(double beta);
  Vector3f getVelocity();
  Vector3f getVelocityChange();
  void detectCollision(CollisionObject *co, double delta_t);

  std::vector<GridNode *> m_neighbors;
  std::vector<GridCell *> m_neighborCells;

  // Physical properties
  Vector3f m_velocity;
  Vector3f m_nextVelocity;

  double m_mass;

private:
  float cubicBSpline(float x) const;
  float gradCubicBSpline(float x) const;

  Vector3i m_idx;
  Grid *m_grid;

  // Physical properties
  // Vector3f m_velocityChange;
  Vector3f m_force;
};

class Grid {
public:
  Grid(Vector3f origin, Vector3i dimensions, float spacing);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Rasterization methods

  void rasterizeMaterialPoints(MaterialPoints &materialPoints);
  void setInitialVolumesAndDensities(MaterialPoints &materialPoints);
  void computeGridForces(MaterialPoints &materialPoints,
                         struct SnowModel snowModel);

  // Helper functions

  // template <typename Lambda>
  // void forEachActiveCell(Lambda &&f) {
  //   for (auto &cell : m_gridCells) {
  //
  //   }
  // }

  template <typename Lambda>
  void forEachNeighbor(MaterialPoint *particle, Lambda &&f) {
    Array3i particleIdx =
        (particle->m_position.array() / m_spacing).floor().cast<int>();
    Array3i min = (particleIdx - 1).max(0);
    Array3i max = (particleIdx + 2).min(m_dim.array() - 1);

    for (int x = min.x(); x < max.x(); x++) {
      for (int y = min.y(); y < max.y(); y++) {
        for (int z = min.z(); z < max.z(); z++) {
          int neighborIdx = vectorToIdx(Vector3i(x, y, z));
          f(m_gridNodes[neighborIdx]);
        }
      }
    }
  }

  std::vector<GridNode *> getAllNodes();

  // Grid sizing and location

  Vector3f m_origin;
  Vector3i m_dim;
  float m_spacing; // h in the paper

private:
  // Utility methods

  inline Vector3i idxToVector(int idx) {
    return Vector3i(idx % m_dim.x(), (idx / m_dim.x()) % m_dim.y(),
                    idx / (m_dim.x() * m_dim.y()));
  }

  inline int vectorToIdx(const Vector3i &idx) {
    return idx.x() + m_dim.x() * (idx.y() + m_dim.y() * idx.z());
  }

  // Grid data
  // the x dimension varies the fastest, followed by y and then z.

  std::vector<GridNode *> m_gridNodes;
  std::vector<GridCell *> m_gridCells;

  // std::shared_ptr<spdlog::logger> logger;
};

} // namespace SnowSimulator

#endif // GRID_H
