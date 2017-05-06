#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

// #include "forceField.hpp"
#include "collisionObject.hpp"
#include "grid.hpp"
#include "materialPoints.hpp"
#include "snowModel.hpp"

namespace SnowSimulator {

class Simulator {
public:
  Simulator(MaterialPoints &materialPoints, Grid *grid,
            std::vector<CollisionObject *> colliders);

  void advance(double timestep, SnowModel snowModel);
  void rasterizeParticlesToGrid();

  void firstStep();

  void updateDeformationGradient(double timestep, SnowModel snowModel);
  void updateParticleVelocities(double timestep, float alpha = 0.95);
  void detectParticleCollisions(double timestep);
  void updateParticlePositions(double timestep);

  std::vector<CollisionObject *> m_colliders;

private:
  // Object storing the particle data
  MaterialPoints &m_materialPoints;
  Grid *m_grid;

  // auto logger = spdlog::get("snowsim");

  // vector of all force fields
  // const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
