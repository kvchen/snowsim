#ifndef PARTICLEGRID_H
#define PARTICLEGRID_H

#include <nanogui/nanogui.h>

class ParticleGrid {
public:
  ParticleGrid();

  ~ParticleGrid(){};

  int index(int i, int j, int k);

protected:
  Eigen::Vector3i m_dim;
  float m_spacing; // h in the paper

  Eigen::VectorXf m_masses;
  Eigen::VectorXf m_velocities;
  Eigen::VectorXf m_prevVelocities;
};

#endif // PARTICLEGRID_H
