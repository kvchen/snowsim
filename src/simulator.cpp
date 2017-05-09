#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "simulator.hpp"
#include "spdlog/spdlog.h"

using namespace SnowSimulator;

Simulator::Simulator(MaterialPoints &materialPoints, Grid *grid,
                     SnowModel &snowModel,
                     std::vector<CollisionObject *> colliders)
    : m_materialPoints(materialPoints), m_grid(grid), m_snowModel(snowModel),
      m_colliders(colliders) {
  logger = spdlog::get("snowsim");
}

/**
 * The simulation proceeds as follows:
 *
 * 1. Rasterize particle data to the grid (rasterizeParticlesToGrid)
 * 2. Compute particle volumes and densities (setParticleVolumesAndDensities)
 *    NOTE: This is done only on the first timestep!
 * 3. Compute grid forces (computeGridForces)
 * 4. Update velocities on grid to v_i' (updateGridVelocities)
 * 5. Grid-based body collisions (detectGridCollisions)
 * 6. Solve linear system for semi-implicit integration (solveLinearSystem)
 * 7. Update deformation gradient (updateDeformationGradient)
 * 8. Update particle velocities (updateParticleVelocities)
 * 9. Particle-based body collisions (detectParticleCollisions)
 * 10. Update particle positions (updateParticlePostions)
 **/
void Simulator::advance(double delta_t) {
  logger->info("[Beginning step {} -> advancing by {}]", stepCount, delta_t);

  rasterizeParticlesToGrid();

  // On the first timestep, we need to compute particle volumes and densities

  if (stepCount == 0) {
    setParticleVolumesAndDensities();
  }

  computeGridForces();
  updateGridVelocities(delta_t);

  detectGridCollisions(delta_t);

  explicitIntegration();
  // solveLinearSystem();

  updateDeformationGradient(delta_t);
  updateParticleVelocities(delta_t);
  detectParticleCollisions(delta_t);
  updateParticlePositions(delta_t);

  stepCount++;
}

void Simulator::rasterizeParticlesToGrid() {
  // Reset all node properties

  for (auto &node : m_grid->nodes()) {
    node->m_mass = 0;
    node->m_velocity.setZero();
    node->m_nextVelocity.setZero();
    node->m_velocityStar.setZero();
    node->m_force.setZero();
  }

  // Place all particles in their updated cells

  // TODO(kvchen): Comment this out once we fix this bug

  // for (auto &mp : m_materialPoints.particles()) {
  //   Vector3f idx = (mp->m_position.array() / m_grid->m_spacing).floor();
  //   int i = m_grid->vectorToIdx(idx.cast<int>());
  //
  //   if (i < 0 || i >= m_grid->m_gridCells.size()) {
  //     logger->error("Particle index {} exceeds grid boundaries", i);
  //     std::cout << *mp << std::endl;
  //
  //     logger->error("Other particles in the same cell as invalid particle:");
  //     for (auto &other : mp->cell()->particles()) {
  //       if (other == mp) {
  //         continue;
  //       }
  //       logger->error("Distance to invalid particle: {}",
  //                     (other->position() - mp->position()).norm());
  //       std::cout << *other << std::endl;
  //     }
  //   }
  // }

  for (auto &cell : m_grid->cells()) {
    cell->clear();
  }

  for (auto &mp : m_materialPoints.particles()) {
    int idx = m_grid->getParticleIdx(mp);
    m_grid->cells()[idx]->addMaterialPoint(mp);
  }

  for (auto &mp : m_materialPoints.particles()) {
    m_grid->forEachNeighbor(mp, [&mp](GridNode *node) {
      double w = node->basisFunction(mp->position());
      node->m_mass += mp->mass() * w;
      node->m_velocity += mp->velocity() * mp->mass() * w;
    });
  }

  for (auto &node : m_grid->nodes()) {
    if (node->m_mass > 0) {
      node->m_velocity /= node->m_mass;
    }
  }
}

void Simulator::setParticleVolumesAndDensities() {
  logger->info("Setting initial particle volumes and densities");

  for (auto const &mp : m_materialPoints.particles()) {
    mp->m_volume = 0;
    mp->m_density = 0;

    m_grid->forEachNeighbor(mp, [&](GridNode *node) {
      mp->m_density += node->m_mass * node->basisFunction(mp->m_position) /
                       pow(m_grid->m_spacing, 3);
    });

    if (mp->m_density != 0) {
      mp->m_volume = mp->m_mass / mp->m_density;
    }
  }
}

void Simulator::computeGridForces() {
  for (auto &mp : m_materialPoints.particles()) {
    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);

    double Jp = mp->m_defPlastic.determinant();
    double Je = mp->m_defElastic.determinant();

    Matrix3f Re = svd.matrixU() * svd.matrixV().transpose();

    double epsilon =
        exp(fmin(m_snowModel.hardeningCoefficient * (1 - Jp), 1e3));

    // Compute the Cauchy stress
    // mu * 2 * (fe - re)_f

    Matrix3f stress =
        epsilon *
        (2 * m_snowModel.initialMu * (mp->m_defElastic - Re) *
             mp->m_defElastic.transpose() +
         Matrix3f::Identity() * (m_snowModel.initialLambda * (Je - 1) * Je));

    m_grid->forEachNeighbor(mp, [mp, &stress](GridNode *node) {
      node->m_force -=
          mp->m_volume * stress * node->gradBasisFunction(mp->m_position);
    });
  }
}

void Simulator::updateGridVelocities(double delta_t) {
  for (auto &node : m_grid->getAllNodes()) {
    if (node->m_mass > 0) {
      node->m_nextVelocity =
          node->m_velocity + delta_t * node->m_force / node->m_mass;
    }

    // Add gravitational forces

    node->m_nextVelocity += Vector3f(0, -9.8, 0) * delta_t;
  }
}

void Simulator::detectGridCollisions(double delta_t) {
  for (auto &node : m_grid->getAllNodes()) {
    node->m_velocityStar = node->m_nextVelocity;

    for (auto &collider : m_colliders) {
      Vector3f position = node->position();
      if (collider->phi(position) > 0) {
        continue;
      }

      Vector3f normal = collider->normal(position);
      Vector3f relVelocity = node->m_velocityStar - collider->velocity();
      double magnitude = relVelocity.dot(normal);

      if (magnitude >= 0) {
        continue;
      }

      Vector3f tangentialVelocity = relVelocity - magnitude * normal;

      // Only apply dynamic friction if the tangential velocity is large
      // compared to the normal

      double normComponent = collider->friction() * magnitude;
      double tangentNorm = tangentialVelocity.norm();

      if (tangentNorm <= -normComponent) {
        relVelocity.setZero();
      } else {
        relVelocity = tangentialVelocity * (1 + normComponent / tangentNorm);
      }

      node->m_velocityStar = relVelocity + collider->velocity();
    }
  }
}

void Simulator::explicitIntegration() {
  for (auto &node : m_grid->getAllNodes()) {
    node->m_nextVelocity = node->m_velocityStar;
    // std::cout << node->m_nextVelocity << std::endl;
  }
}

void Simulator::updateDeformationGradient(double delta_t) {
  for (auto &mp : m_materialPoints.particles()) {
    Matrix3f gradVelocity = Matrix3f::Identity();
    m_grid->forEachNeighbor(mp, [mp, delta_t, &gradVelocity](GridNode *node) {
      gradVelocity += delta_t * node->nextVelocity() *
                      node->gradBasisFunction(mp->m_position).transpose();
    });

    Matrix3f defNext = gradVelocity * mp->m_defElastic * mp->m_defPlastic;
    mp->m_defElastic = gradVelocity * mp->m_defElastic;

    // Push deformations exceeding critical deformation thresholds into Fp

    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    Matrix3f sigma = svd.singularValues()
                         .cwiseMax(1 - m_snowModel.criticalCompression)
                         .cwiseMin(1 + m_snowModel.criticalStretch)
                         .asDiagonal();

    mp->m_defPlastic =
        svd.matrixV() * sigma.inverse() * svd.matrixU().transpose() * defNext;
    mp->m_defElastic = svd.matrixU() * sigma * svd.matrixV().transpose();
    if (!defNext.isApprox(mp->m_defElastic * mp->m_defPlastic)) {
      logger->warn("Deformation update is incorrect");
    }
  }
}

void Simulator::updateParticleVelocities(double delta_t, float alpha) {
  for (auto &mp : m_materialPoints.particles()) {
    Vector3f velocityPIC = Vector3f::Zero();
    Vector3f velocityFLIP = mp->m_velocity;

    m_grid->forEachNeighbor(
        mp, [mp, &velocityPIC, &velocityFLIP](GridNode *node) {
          float weight = node->basisFunction(mp->m_position);

          velocityPIC += node->m_nextVelocity * weight;
          velocityFLIP += (node->m_nextVelocity - node->m_velocity) * weight;
        });

    mp->m_velocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
  }
}

void Simulator::detectParticleCollisions(double delta_t) {
  for (auto &mp : m_materialPoints.particles()) {
    for (auto &collider : m_colliders) {
      Vector3f position = mp->m_position;
      if (collider->phi(position) > 0) {
        continue;
      }

      // std::cout << mp->position() << std::endl;

      Vector3f normal = collider->normal(position);
      Vector3f relVelocity = mp->velocity() - collider->velocity();

      // std::cout << collider->velocity() << std::endl;

      double magnitude = relVelocity.dot(normal);
      if (magnitude >= 0) {
        continue;
      }

      Vector3f tangentialVelocity = relVelocity - magnitude * normal;

      // Only apply dynamic friction if the tangential velocity is large
      // compared to the normal

      double normComponent = collider->friction() * magnitude;
      double tangentNorm = tangentialVelocity.norm();

      // if (tangentNorm <= -normComponent) {
      relVelocity.setZero();
      // } else {
      //   std::cout << "GOT HERE" << std::endl;
      //   relVelocity = tangentialVelocity * (1 + normComponent /
      //   tangentNorm);
      // }

      // std::cout << "BEFORE: " << mp->m_velocity << std::endl;
      // mp->m_velocity = relVelocity + collider->velocity();
      // std::cout << "AFTER: " << mp->m_velocity << std::endl;
    }
  }
}

/**
 * Uses backwards Euler integration to update the particle position for each
 * delta_t. At this point, the velocities should already have been updated
 * to the next delta_t.
 */
void Simulator::updateParticlePositions(double delta_t) {
  for (auto &mp : m_materialPoints.particles()) {
    mp->m_position += delta_t * mp->m_velocity;

    // mp->m_position = mp->m_position.array().min(
    //     (m_grid->m_dim.array()).cast<float>() * m_grid->m_spacing);
    // mp->m_position = mp->m_position.cwiseMax(0);
  }
}
