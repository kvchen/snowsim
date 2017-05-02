#ifndef MATERIALPOINTS_H
#define MATERIALPOINTS_H

#include "forceField.hpp"

namespace SnowSimulator {

/**
 * Stores data about all the particles used in the simulation. This should
 * probably be a struct instead.
 */
class MaterialPoints {

public:
  MaterialPoints();
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::VectorXf m_mass;
  Eigen::VectorXf m_volumes;

  //
  Eigen::VectorXF m_weights;

  // (3 x num_particles) matrices of particle properties. Each particle
  // contains a number of
  Eigen::MatrixXf m_positions;
  Eigen::MatrixXf m_velocities;
  Eigen::MatrixXf m_deformationGradients;
};

} // namespace SnowSimulator

#endif // MATERIALPOINTS_H
