#ifndef MATERIALPOINTS_H
#define MATERIALPOINTS_H

#include <vector>

#include "forceField.hpp"
#include "grid.hpp"

using namespace Eigen;

namespace SnowSimulator {

class GridCell;

class MaterialPoint;
class MaterialPoints;

class MaterialPoint {
public:
  MaterialPoint(Vector3f &position, Vector3f &velocity)
      : m_position(position), m_velocity(velocity) {
    m_defElastic = Matrix3f::Identity();
    m_defPlastic = Matrix3f::Identity();
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Vector3f m_position;
  Vector3f m_velocity;
  Matrix3f m_defElastic;
  Matrix3f m_defPlastic;

  GridCell *m_cell;

  double m_mass;
  double m_volume;
  double m_density;
};

/**
 * Stores data about all the particles used in the simulation.
 */
class MaterialPoints {
public:
  MaterialPoints() {};

  std::vector<MaterialPoint *> m_materialPoints;
};

} // namespace SnowSimulator

#endif // MATERIALPOINTS_H
