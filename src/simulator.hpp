#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

// #include "forceField.hpp"
#include "collisionObject.hpp"
#include "grid.hpp"
#include "materialPoints.hpp"
#include "snowModel.hpp"
#include "spdlog/spdlog.h"

namespace SnowSimulator {

class Simulator {
public:
  Simulator(MaterialPoints &materialPoints, Grid *grid,
            std::vector<CollisionObject *> colliders);

  void advance(double delta_t, SnowModel snowModel);
  void rasterizeParticlesToGrid();

  void firstStep();

  void updateDeformationGradient(double delta_t, SnowModel snowModel);
  void updateParticleVelocities(double delta_t, float alpha = 0.95);
  void detectParticleCollisions(double delta_t);
  void updateParticlePositions(double delta_t);

  std::vector<CollisionObject *> m_colliders;

private:
  // Object storing the particle data
  MaterialPoints &m_materialPoints;
  Grid *m_grid;

  size_t stepCount;

  std::shared_ptr<spdlog::logger> logger;

  // vector of all force fields
  // const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
