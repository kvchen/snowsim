#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <Eigen/Dense>

namespace SnowSimulator {

class ForceField {
public:
  virtual Eigen::Vector3f force(const Eigen::Vector3f &pos, float mass) const;
  virtual Eigen::Matrix3f dForce(const Eigen::Vector3f &pos, float mass) const;

private:
};

} // namespace SnowSimulator

#endif // FORCEFIELD_H
