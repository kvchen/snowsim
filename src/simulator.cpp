#include <Eigen/Dense>
#include <iostream>

#include "simulator.hpp"
#include "spdlog/spdlog.h"

using namespace SnowSimulator;

Simulator::Simulator(MaterialPoints &materialPoints, Grid *grid,
                     std::vector<CollisionObject *> colliders)
    : m_materialPoints(materialPoints), m_grid(grid), m_colliders(colliders) {
  logger = spdlog::get("snowsim");
}

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
void Simulator::advance(double delta_t, SnowModel snowModel) {
  logger->info("Advancing by {}", delta_t);

  logger->info("Rasterizing material point data to grid");
  m_grid->rasterizeMaterialPoints(m_materialPoints);

  logger->info("Computing grid forces");
  m_grid->computeGridForces(m_materialPoints, snowModel);

  for (auto &node : m_grid->getAllNodes()) {
    node->explicitUpdateVelocity(delta_t);
    // for (auto &co : m_colliders) {
    //   node->detectCollision(co, delta_t);
    // }
  }

  // solveLinearSystem();
  updateDeformationGradient(delta_t, snowModel);
  updateParticleVelocities(delta_t);
  detectParticleCollisions(delta_t);
  updateParticlePositions(delta_t);

  stepCount++;
}

// void Simulator::updateParticlePlacement() {
//   logger->info("Clearing cells");
//   for (auto &cell : m_grid->m_gridCells) {
//     cell->clear();
//   }
//
//   logger->info("Adding material points to cells...");
//   for (auto &mp : m_materialPoints.m_materialPoints) {
//     Vector3f idx = (mp->m_position.array() / m_spacing).floor();
//     int i = vectorToIdx(idx.cast<int>());
//     m_gridCells[i]->addMaterialPoint(mp);
//   }
// }
//
// void Simulator::rasterizeParticlesToGrid() {
//   logger->info("Rasterizing material points to grid nodes...");
//   for (auto &node : m_gridNodes) {
//     node->rasterizeMaterialPoints();
//   }
//
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
void Simulator::updateDeformationGradient(double delta_t, SnowModel snowModel) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Matrix3f gradVelocity = Matrix3f::Identity();
    m_grid->forEachNeighbor(mp, [mp, &gradVelocity](GridNode *node) {
      gradVelocity += node->getVelocity() *
                      node->gradBasisFunction(mp->m_position).transpose();
    });

    Matrix3f defUpdate = mp->m_defElastic * mp->m_defPlastic;

    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    // Vector3f sv = svd.singularValues();

    Matrix3f sigma = svd.singularValues()
                         .cwiseMax(1 - snowModel.criticalCompression)
                         .cwiseMin(1 + snowModel.criticalStretch)
                         .asDiagonal();

    mp->m_defElastic = svd.matrixU() * sigma * svd.matrixV().transpose();
    mp->m_defPlastic =
        svd.matrixV() * sigma.inverse() * svd.matrixU().transpose() * defUpdate;
  }
}

/**
 * Update the particle velocities according to PIC and FLIP.
 */
void Simulator::updateParticleVelocities(double delta_t, float alpha) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Vector3f velocityPIC = Vector3f::Zero();
    Vector3f velocityFLIP = mp->m_velocity;

    m_grid->forEachNeighbor(
        mp, [mp, &velocityPIC, &velocityFLIP](GridNode *node) {
          float weight = node->basisFunction(mp->m_position);
          velocityPIC += node->m_velocity * weight;
          velocityFLIP += (node->m_nextVelocity - node->m_velocity) * weight;
        });

    mp->m_velocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
    mp->m_velocity += Vector3f(0, -9.8, 0);
  }
}

void Simulator::detectParticleCollisions(double delta_t) {
  for (auto &co : m_colliders) {
    for (auto &mp : m_materialPoints.m_materialPoints) {
      Vector3f position = mp->m_position + delta_t * mp->m_velocity;
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
 * delta_t. At this point, the velocities should already have been updated
 * to the next delta_t.
 */
void Simulator::updateParticlePositions(double delta_t) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    mp->m_position += delta_t * mp->m_velocity;
    mp->m_position = mp->m_position.array().min(
        (m_grid->m_dim.array() - 1).cast<float>() * m_grid->m_spacing);
    mp->m_position = mp->m_position.cwiseMax(0);

    // assert(!mp->m_velocity.array().isNaN().any());

    // if (!mp->m_velocity.isZero())
    // logger->info("Velocity ({}, {}, {}) should be 0", mp->m_velocity.x(),
    //              mp->m_velocity.y(), mp->m_velocity.z());
  }
}
