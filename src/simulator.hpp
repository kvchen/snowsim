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
  Simulator(MaterialPoints &materialPoints, Grid *grid, SnowModel &snowModel,
            std::vector<CollisionObject *> colliders);

  size_t stepCount() { return m_stepCount; }
  std::string exportVolumeData();

  void advance(double delta_t);
  void rasterizeParticlesToGrid();
  void setParticleVolumesAndDensities();
  void computeGridForces();
  void updateGridVelocities(double delta_t);

  void detectGridCollisions(double delta_t);
  void explicitIntegration();

  void updateDeformationGradient(double delta_t);
  void updateParticleVelocities(double delta_t, float alpha = 0.95);
  void detectParticleCollisions(double delta_t);
  void updateParticlePositions(double delta_t);

  std::vector<CollisionObject *> m_colliders;

private:
  // Object storing the particle data
  MaterialPoints &m_materialPoints;
  Grid *m_grid;
  SnowModel &m_snowModel;

  size_t m_stepCount;

  std::shared_ptr<spdlog::logger> logger;

  // vector of all force fields
  // const std::vector<ForceField> &m_forceFields;
};

} // namespace SnowSimulator

#endif // SIMULATOR_H
