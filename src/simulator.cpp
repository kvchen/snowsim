#include <Eigen/Dense>

#include "simulator.hpp"
#include "spdlog/spdlog.h"

using namespace SnowSimulator;

Simulator::Simulator(MaterialPoints &materialPoints, Grid *grid,
                     std::vector<CollisionObject *> colliders)
    : m_materialPoints(materialPoints), m_grid(grid), m_colliders(colliders) {}

/**
 * The simulation proceeds as follows:
 *
 * Rasterize particle data to the grid (rasterizeToGrid)
 * Compute particle volumes and densities. This is done only once!
 *   (computeParticleProperties)
 * Compute grid forces (computeGridForces)
 * Update velocities on grid to v_i' (updateGridVelocities)
 * Grid-based body collisions (computeBodyCollisions)
 * Solve linear system for semi-implicit integration (solveLinearSystem)
 * Update deformation gradient (updateDeformationGradient)
 * Update particle velocities (updateParticleVelocities)
 * Particle-based body collisions (computeParticleBodyCollisions)
 * Update particle positions (updateParticlePostions)
 **/
void Simulator::advance(double timestep, SnowModel snowModel) {
  auto logger = spdlog::get("snowsim");
  logger->info("Advancing by {}", timestep);

  logger->info("Rasterizing material point data to grid");
  m_grid->rasterizeMaterialPoints(m_materialPoints);

  logger->info("Computing grid forces");
  m_grid->computeGridForces(m_materialPoints, snowModel);

  for (auto &node : m_grid->getAllNodes()) {
    node->explicitUpdateVelocity(timestep);
    for (auto &co : m_colliders) {
      node->detectCollision(co, timestep);
    }
  }

  // solveLinearSystem();
  updateDeformationGradient(timestep, snowModel);
  updateParticleVelocities(timestep);
  detectParticleCollisions(timestep);
  updateParticlePositions(timestep);
}

void Simulator::rasterizeParticlesToGrid() {}

// void Grid::rasterizeMaterialPoints(MaterialPoints &materialPoints) {
//   auto logger = spdlog::get("snowsim");
//   logger->info("Clearing cells...");
//   for (auto &cell : m_gridCells) {
//     cell->clear();
//   }
//   logger->info("Adding material points to cells...");
//   for (auto &mp : materialPoints.m_materialPoints) {
//     Vector3f idx = (mp->m_position.array() / m_spacing).floor();
//     int i = vectorToIdx(idx.cast<int>());
//     m_gridCells[i]->addMaterialPoint(mp);
//   }
//   logger->info("Rasterizing material points to grid nodes...");
//   for (auto &node : m_gridNodes) {
//     node->rasterizeMaterialPoints();
//   }
//   logger->info("Finished rasterizing");
// }

void Simulator::firstStep() {
  m_grid->rasterizeMaterialPoints(m_materialPoints);
  m_grid->setInitialVolumesAndDensities(m_materialPoints);
}

//
// Simulator::computeWeights() {}
//
// void Simulator::rasterizeToGrid() {
//   std::vector<GridCell *> cells = m_grid->getGridCells();
//   for (auto &cell : cells) {
//     cell->clear();
//   }
//   for (auto &mp : m_materialPoints.m_materialPoints) {
//     Vector3f idx = (mp->m_position / m_grid->m_spacing).floor();
//     GridCell *cell = cells[m_grid->vectorToIdx(idx)];
//     cell->addMaterialPoint(mp);
//   }
//   for (auto &node : m_grid->getAllNodes()) {
//     node->rasterizeMaterialPoints();
//   }
// }
//
// Simulator::computeParticleProperties() {}
//
// Simulator::computeGridForces() {}
//
// Simulator::updateGridVelocities() {}
//
// Simulator::computeBodyCollisions() {}
//
// Simulator::solveLinearSystem() {}
//
void Simulator::updateDeformationGradient(double timestep,
                                          SnowModel snowModel) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Matrix3f gradVelocity = Matrix3f::Zero();

    m_grid->forEachNeighbor(mp, [mp, &gradVelocity](GridNode *node) {
      gradVelocity += node->getVelocity() *
                      node->gradBasisFunction(mp->m_position).transpose();
    });

    // for (auto &node : m_grid->getNearbyNodes(mp)) {
    //   gradVelocity += node->getVelocity() *
    //                   node->gradBasisFunction(mp->m_position).transpose();
    // }

    mp->m_defElastic =
        (Matrix3f::Identity() + timestep * gradVelocity) * mp->m_defElastic;
    Matrix3f defUpdate = mp->m_defElastic * mp->m_defPlastic;

    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    Vector3f sigma = svd.singularValues();
    for (int i = 0; i < 3; i++) {
      if (sigma[i] > 1 - snowModel.criticalCompression)
        sigma[i] = 1 - snowModel.criticalCompression;
      if (sigma[i] < 1 + snowModel.criticalStretch)
        sigma[i] = 1 + snowModel.criticalStretch;
    }

    mp->m_defElastic =
        svd.matrixU() * sigma.asDiagonal() * svd.matrixV().transpose();
    mp->m_defPlastic = svd.matrixV() * sigma.asDiagonal().inverse() *
                       svd.matrixU().transpose() * defUpdate;
  }
}

/**
 * Update the particle velocities according to PIC and FLIP.
 */
void Simulator::updateParticleVelocities(double timestep, float alpha) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Vector3f velocityPIC = Vector3f::Zero();
    Vector3f velocityFLIP = mp->m_velocity;

    // for (auto &node : m_grid->getNearbyNodes(mp)) {
    //   float weight = node->basisFunction(mp->m_position);
    //   velocityPIC += node->getVelocity() * weight;
    //   velocityFLIP += node->getVelocityChange() * weight;
    // }

    m_grid->forEachNeighbor(
        mp, [mp, &velocityPIC, &velocityFLIP](GridNode *node) {
          float weight = node->basisFunction(mp->m_position);
          velocityPIC += node->getVelocity() * weight;
          velocityFLIP += node->getVelocityChange() * weight;
        });

    mp->m_velocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
  }
}

void Simulator::detectParticleCollisions(double timestep) {
  for (auto &co : m_colliders) {
    for (auto &mp : m_materialPoints.m_materialPoints) {
      Vector3f position = mp->m_position + timestep * mp->m_velocity;
      if (co->phi(position) <= 0) {
        Vector3f normal = co->normal(position);
        Vector3f relVelocity = mp->m_velocity - co->m_velocity;
        double magnitude = relVelocity.transpose() * normal;
        if (magnitude < 0) {
          Vector3f tangent = relVelocity - normal * magnitude;
          if (tangent.norm() <= -co->m_friction * magnitude) {
            relVelocity.setZero();
          } else {
            relVelocity =
                tangent + co->m_friction * magnitude * tangent / tangent.norm();
          }
        }
        mp->m_velocity = relVelocity + co->m_velocity;
      }
    }
  }
}

/**
 * Uses backwards Euler integration to update the particle position for each
 * timestep. At this point, the velocities should already have been updated
 * to the next timestep.
 */
void Simulator::updateParticlePositions(double timestep) {
  auto logger = spdlog::get("snowsim");
  for (auto &mp : m_materialPoints.m_materialPoints) {
    mp->m_position += timestep * mp->m_velocity;
    if (!mp->m_velocity.isZero())
      logger->info("Velocity ({}, {}, {}) should be 0", mp->m_velocity.x(),
                   mp->m_velocity.y(), mp->m_velocity.z());
  }
}
