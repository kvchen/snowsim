#include "particleGrid.hpp"

int ParticleGrid::index(int i, int j, int k) {
  return i * (m_dim.y() * m_dim.z()) + j * (m_dim.z()) + k;
}
