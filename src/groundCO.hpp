#ifndef GROUNDCO_H
#define GROUNDCO_H

#include <Eigen/Dense>

namespace SnowSimulator {

class GroundCO: public CollisionObject {
public:
  GroundCO(double friction) : m_friction(friction) {};
  float phi(const Eigen::Vector3f &x) const {if (x[1] <= 0) return -1; else return 1;};
  Eigen::Vector3f normal(const Eigen::Vector3f &x) {return Eigen::Vector3f(0, 1, 0);};

  Eigen::Vector3f m_velocity = Eigen::Vector3f::Zero();
  double m_friction;
};

} // namespace SnowSimulator

#endif // GROUNDCO_H
