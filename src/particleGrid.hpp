#ifndef PARTICLEGRID_H
#define PARTICLEGRID_H

#include <Eigen/Dense>

namespace SnowSimulator {

class ParticleGrid {
public:
  ParticleGrid();
  ~ParticleGrid(){};

  float basisFunction(float x);
  float dBasisFunction(float x);

  int index(int i, int j, int k);

private:
  // Grid sizing

  Eigen::Vector3i m_bbox_min;
  Eigen::Vector3i m_bbox_max;
  float m_spacing; // h in the paper

  // Eigen::VectorXf m_masses;
  // Eigen::VectorXf m_velocities;
  // Eigen::VectorXf m_prevVelocities;
};

} // namespace SnowSimulator

#endif // PARTICLEGRID_H
