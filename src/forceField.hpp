#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <Eigen/Dense>

namespace SnowSimulator {

class ForceField {
public:
  virtual Eigen::Vector3f force(const Eigen::Vector3f &pos, float mass) const;
  virtual Eigen::Matrix3f gradForce(const Eigen::Vector3f &pos,
                                    float mass) const;
};

} // namespace SnowSimulator

#endif // FORCEFIELD_H
