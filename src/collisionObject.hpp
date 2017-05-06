#ifndef COLLISIONOBJECT_H
#define COLLISIONOBJECT_H

#include <Eigen/Dense>

namespace SnowSimulator {

class CollisionObject {
public:
  virtual float phi(const Eigen::Vector3f &x) const = 0;
  virtual Eigen::Vector3f normal(const Eigen::Vector3f &x) {return Eigen::Vector3f::Zero();};

  Eigen::Vector3f m_velocity;
  double m_friction;
};

} // namespace SnowSimulator

#endif // COLLISIONOBJECT_H
