#ifndef COLLISIONOBJECT_H
#define COLLISIONOBJECT_H

#include <Eigen/Dense>

namespace SnowSimulator {

class CollisionObject {
public:
  virtual float phi(const Eigen::Vector3f &x) const = 0;

private:
};

} // namespace SnowSimulator

#endif // COLLISIONOBJECT_H
