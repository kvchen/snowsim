#ifndef MATERIALPOINTS_H
#define MATERIALPOINTS_H

#include <iostream>
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

  double m_mass = 1.0;
  double m_volume;
  double m_density;

  friend std::ostream &operator<<(std::ostream &os, const MaterialPoint &mp) {
    IOFormat InlineFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "",
                          "(", ")");
    os << "Particle(" << std::endl;
    os << "  position: " << mp.m_position.format(InlineFormat) << std::endl;
    os << "  velocity: " << mp.m_velocity.format(InlineFormat) << std::endl;
    os << "  defElastic: " << mp.m_defElastic.format(InlineFormat) << std::endl;
    os << "  defPlastic: " << mp.m_defPlastic.format(InlineFormat) << std::endl;
    os << "  mass: " << mp.m_mass << std::endl;
    os << "  volume: " << mp.m_volume << std::endl;
    os << "  density: " << mp.m_density << std::endl;
    os << ")" << std::endl;

    return os;
  };
};

/**
 * Stores data about all the particles used in the simulation.
 */
class MaterialPoints {
public:
  MaterialPoints(){};

  std::vector<MaterialPoint *> m_materialPoints;
};

} // namespace SnowSimulator

#endif // MATERIALPOINTS_H
