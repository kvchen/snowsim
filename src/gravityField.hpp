#ifndef GRAVITYFIELD_H
#define GRAVITYFIELD_H

#include <Eigen/Dense>

#include "forceField.hpp"

namespace SnowSimulator {

class GravityField : public ForceField {
public:
  GravityField(const Eigen::Vector3f &gravity);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::Vector3f force(const Eigen::Vector3f &pos, float mass) const;
  Eigen::Matrix3f gradForce(const Eigen::Vector3f &pos, float mass) const;

private:
  Eigen::Vector3f m_gravity;
};

} // namespace SnowSimulator

#endif // GRAVITYFIELD_H
