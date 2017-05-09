#ifndef SNOWSIM_COLLISIONOBJECT_H
#define SNOWSIM_COLLISIONOBJECT_H

#include <Eigen/Dense>

using namespace Eigen;

namespace SnowSimulator {

class CollisionObject {
public:
  virtual double phi(const Vector3f &x) = 0;

  virtual Vector3f normal(const Vector3f &x) = 0;

  virtual Vector3f velocity() = 0;

  virtual double friction() = 0;
};

} // namespace SnowSimulator

#endif // SNOWSIM_COLLISIONOBJECT_H
