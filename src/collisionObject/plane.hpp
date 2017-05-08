#ifndef SNOWSIM_PLANE_COLLISIONOBJECT_H
#define SNOWSIM_PLANE_COLLISIONOBJECT_H

#include <Eigen/Dense>

#include "collisionObject.hpp"

using namespace Eigen;

namespace SnowSimulator {

class Plane : public CollisionObject {
public:
  Plane(const Vector3f &point, const Vector3f &normal, double friction)
      : m_point(point), m_normal(normal.normalized()), m_friction(friction) {
    m_velocity.setZero();
  }

  float phi(const Eigen::Vector3f &x) const {
    return m_normal.dot(x - m_point);
  }

  Vector3f normal(const Eigen::Vector3f &x) const { return m_normal; }

  Vector3f velocity() { return m_velocity; }

  double friction() { return m_friction; }

private:
  Vector3f m_point;
  Vector3f m_normal;

  Vector3f m_velocity;
  double m_friction;
};

} // namespace SnowSimulator

#endif // SNOWSIM_PLANE_COLLISIONOBJECT_H
