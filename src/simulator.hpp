#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

#include "forceField.hpp"

namespace SnowSimulator {

class Simulator {
public:
  Simulator();
  void advance(double timestep);

private:
  const MaterialPoints &m_materialPoints;

  // vector of all force fields
  const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
