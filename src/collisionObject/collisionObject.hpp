#ifndef SNOWSIM_COLLISIONOBJECT_H
#define SNOWSIM_COLLISIONOBJECT_H

#include <Eigen/Dense>

using namespace Eigen;

namespace SnowSimulator {

class CollisionObject {
public:
  virtual float phi(const Vector3f &x) const = 0;

  virtual Vector3f normal(const Vector3f &x) const {
    return Eigen::Vector3f::Zero();
  };

  Vector3f velocity() { return m_velocity; }

  double friction() { return m_friction; }

private:
  Vector3f m_velocity;
  double m_friction;
};

} // namespace SnowSimulator

#endif // SNOWSIM_COLLISIONOBJECT_H
