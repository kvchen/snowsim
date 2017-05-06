#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

// #include "forceField.hpp"
#include "materialPoints.hpp"
#include "snowModel.hpp"
#include "grid.hpp"

namespace SnowSimulator {

class Simulator {
public:
  Simulator(MaterialPoints &materialPoints, Grid *grid);
  // void advance(double timestep);
  void updateDeformationGradient(double timestep, SnowModel snowModel);
  void updateParticleVelocities(double timestep, float alpha = 0.95);

private:
  // Object storing the particle data
  const MaterialPoints &m_materialPoints;
  Grid *m_grid;

  // vector of all force fields
  // const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
