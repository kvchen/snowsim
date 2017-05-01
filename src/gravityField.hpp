#ifndef GRAVITYFIELD_H
#define GRAVITYFIELD_H

#include "forceField.hpp"

namespace SnowSimulator {

class GravityField : public ForceField {
public:
  GravityField(const Eigen::Vector3f &gravity);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  virtual Eigen::Vector3f force(const Eigen::Vector3f &pos, float mass) const;
  virtual Eigen::Matrix3f dForce(const Eigen::Vector3f &pos, float mass) const;

private:
  Eigen::Vector3f m_gravity;
};

} // namespace SnowSimulator

#endif // GRAVITYFIELD_H
