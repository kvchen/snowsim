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

  std::vector<MaterialPoint *> &particles() { return m_materialPoints; }

  GridNode *m_node; // Pointer to the immediately adjacent GridNode
  std::vector<MaterialPoint *> m_materialPoints;
};

class GridNode {
public:
  GridNode(Vector3i idx, Grid *grid);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Vector3f position() const;

  float basisFunction(Vector3f particlePos) const;
  Vector3f gradBasisFunction(Vector3f particlePos) const;

  void zeroForce();
  void addForce(Vector3f force);

  Vector3f &velocity() { return m_velocity; };
  Vector3f &nextVelocity() { return m_nextVelocity; };
  Vector3f &force() { return m_force; };

  std::vector<GridNode *> m_neighbors;
  std::vector<GridCell *> m_neighborCells;

  // Physical properties

  Vector3f m_velocity;
  Vector3f m_nextVelocity;
  Vector3f m_velocityStar;

  double m_mass;
  Vector3f m_force;

private:
  float cubicBSpline(float x) const;
  float gradCubicBSpline(float x) const;

  Vector3i m_idx;
  Grid *m_grid;
};

class Grid {
public:
  Grid(Vector3f origin, Vector3i dimensions, float spacing);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  std::vector<GridNode *> getAllNodes();

  // Template functions

  template <typename Lambda>
  void forEachNeighbor(MaterialPoint *particle, Lambda &&f) {
    Array3i particleIdx =
        (particle->m_position.array() / m_spacing).floor().cast<int>();
    Array3i min = (particleIdx - 2).max(0);
    Array3i max = (particleIdx + 3).min(m_dim.array() - 1);

    for (int x = min.x(); x < max.x(); x++) {
      for (int y = min.y(); y < max.y(); y++) {
        for (int z = min.z(); z < max.z(); z++) {
          int neighborIdx = vectorToIdx(Vector3i(x, y, z));
          f(m_gridNodes[neighborIdx]);
        }
      }
    }
  }

  inline Vector3i idxToVector(int idx) {
    return Vector3i(idx % m_dim.x(), (idx / m_dim.x()) % m_dim.y(),
                    idx / (m_dim.x() * m_dim.y()));
  }

  inline int vectorToIdx(const Vector3i &idx) {
    return idx.x() + m_dim.x() * (idx.y() + m_dim.y() * idx.z());
  }

  int getParticleIdx(MaterialPoint *mp) {
    Vector3i idx = (mp->position().array() / m_spacing).cast<int>();
    return vectorToIdx(idx);
  }

  // Accessor methods

  Vector3f &origin() { return m_origin; }
  Vector3i &dim() { return m_dim; }

  float spacing() { return m_spacing; }

  std::vector<GridNode *> &nodes() { return m_gridNodes; }
  std::vector<GridCell *> &cells() { return m_gridCells; }

  // Grid sizing and location

  Vector3f m_origin;
  Vector3i m_dim;
  float m_spacing; // h in the paper

  // Grid data

  std::vector<GridNode *> m_gridNodes;
  std::vector<GridCell *> m_gridCells;
};

} // namespace SnowSimulator

#endif // GRID_H
