#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

// #include "forceField.hpp"
#include "materialPoints.hpp"

namespace SnowSimulator {

class Simulator {
public:
  Simulator();
  void advance(double timestep);

private:
  // Object storing the particle data
  const MaterialPoints &m_materialPoints;

  // vector of all force fields
  // const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
