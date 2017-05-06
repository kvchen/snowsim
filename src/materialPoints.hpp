#ifndef MATERIALPOINTS_H
#define MATERIALPOINTS_H

#include <boost/functional/hash.hpp>
#include <vector>

// #include "forceField.hpp"
#include "grid.hpp"

using namespace Eigen;

namespace std {

// template <> struct hash<SnowSimulator::MaterialPoint> {
//   std::size_t operator()(const MaterialPoint &mp) const {
//     using boost::hash_value;
//     using boost::hash_combine;
//
//     std::size_t seed = 0;
//
//     hash_combine(seed, hash_value(mp.m_mass));
//     hash_combine(seed, hash_value(mp.m_volume));
//     hash_combine(seed, hash_value(mp.m_density));
//     hash_combine(seed, hash_value(mp.m_density));
//     hash_combine(seed, hash_value(mp.m_density));
//     hash_combine(seed, hash_value(mp.m_density));
//
//     // Return the result.
//     return seed;
//   }
// };

} // namespace std

namespace SnowSimulator {

// Forward declarations

class GridCell;
class MaterialPoint;
class MaterialPoints;

class MaterialPoint {
public:
  MaterialPoint(Vector3f &position, Vector3f &velocity)
      : m_position(position), m_velocity(velocity) {
    m_deformationGradient = Matrix3f::Identity();
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Used for hashing this material point

  bool operator==(const MaterialPoint &other) const {
    return (m_mass == other.m_mass && m_volume == other.m_volume &&
            m_density == other.m_density && m_position == other.m_position &&
            m_velocity == other.m_velocity &&
            m_deformationGradient == other.m_deformationGradient);
  }

  GridCell *m_cell;

  Vector3f m_position;
  Vector3f m_velocity;
  Matrix3f m_deformationGradient;

  double m_mass;
  double m_volume;
  double m_density;
};

/**
 * Stores data about all the particles used in the simulation.
 */
class MaterialPoints {
public:
  MaterialPoints();

  std::vector<MaterialPoint *> m_materialPoints;
};

} // namespace SnowSimulator

#endif // MATERIALPOINTS_H
