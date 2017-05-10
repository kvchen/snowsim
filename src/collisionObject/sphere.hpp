#ifndef SNOWSIM_SPHERE_COLLISIONOBJECT_H
#define SNOWSIM_SPHERE_COLLISIONOBJECT_H

#include <Eigen/Dense>

#include "collisionObject.hpp"

using namespace Eigen;

namespace SnowSimulator {

class Sphere : public CollisionObject {
public:
  Sphere(const Vector3f &center, const Vector3f &radius, double friction)
      : m_point(point), m_normal(normal.normalized()), m_friction(friction) {
    m_velocity.setZero();
  }

  double phi(const Eigen::Vector3f &x) { return m_normal.dot(x - m_point); }

  Vector3f normal(const Eigen::Vector3f &x) { return m_normal; }

  Vector3f velocity() { return m_velocity; }

  double friction() { return m_friction; }

protected:
  Vector3f m_center;
  Vector3f m_radius;

  Vector3f m_velocity;
  double m_friction;
};

} // namespace SnowSimulator

#endif // SNOWSIM_SPHERE_COLLISIONOBJECT_H
