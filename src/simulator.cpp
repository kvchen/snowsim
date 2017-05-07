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

  // computeGridForces();
  updateGridVelocities(delta_t);

  // detectGridCollisions(delta_t);

  // explicitIntegration();

  // solveLinearSystem();
  // updateDeformationGradient(delta_t);
  updateParticleVelocities(delta_t);
  // detectParticleCollisions(delta_t);
  updateParticlePositions(delta_t);

  stepCount++;
}

void Simulator::rasterizeParticlesToGrid() {
  // Reset all node properties

  for (auto &node : m_grid->m_gridNodes) {
    node->m_mass = 0;
    node->m_velocity.setZero();
    node->m_nextVelocity.setZero();
    node->m_force.setZero();
  }

  // Place all particles in their updated cells

  for (auto &cell : m_grid->m_gridCells) {
    cell->clear();
  }

  for (auto &mp : m_materialPoints.m_materialPoints) {
    Vector3f idx = (mp->m_position.array() / m_grid->m_spacing).floor();
    int i = m_grid->vectorToIdx(idx.cast<int>());

    // TODO(kvchen): Comment this out once we fix this bug

    if (i < 0 || i >= m_grid->m_gridCells.size()) {
      logger->error("Particle index {} exceeds grid boundaries", i);
      std::cout << *mp << std::endl;
    }

    m_grid->m_gridCells[i]->addMaterialPoint(mp);
  }

  for (auto &mp : m_materialPoints.m_materialPoints) {
    m_grid->forEachNeighbor(mp, [&mp](GridNode *node) {
      double w = node->basisFunction(mp->m_position);
      node->m_mass += mp->m_mass * w;
      node->m_velocity += mp->m_velocity * mp->m_mass * w;
    });
  }

  for (auto &node : m_grid->m_gridNodes) {
    if (node->m_mass > 0) {
      node->m_velocity /= node->m_mass;
    }
  }
}

void Simulator::setParticleVolumesAndDensities() {
  logger->info("Setting initial particle volumes and densities");

  for (auto const &mp : m_materialPoints.m_materialPoints) {
    mp->m_volume = 0;
    mp->m_density = 0;

    m_grid->forEachNeighbor(mp, [&mp](GridNode *node) {
      mp->m_density += node->m_mass * node->basisFunction(mp->m_position);
    });

    // std::cout << mp->m_mass << std::endl;

    // assert(!std::isnan(mp->m_density));

    if (mp->m_density != 0) {
      mp->m_volume = mp->m_mass / mp->m_density;
    }

    // assert(!std::isnan(mp->m_volume));
  }
}

void Simulator::computeGridForces() {
  for (auto &node : m_grid->m_gridNodes) {
    node->zeroForce();
    // node->addForce(Vector3f(0, -9.8, 0));
  }

  for (auto &mp : m_materialPoints.m_materialPoints) {
    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);

    double Jp = mp->m_defPlastic.determinant();
    double Je = mp->m_defElastic.determinant();

    Matrix3f Re = svd.matrixU() * svd.matrixV().transpose();
    // Matrix3f Se = svd.matrixV() * svd.singularValues().asDiagonal() *
    //               svd.matrixV().transpose();
    // Matrix3f Fe = Re * Se;

    double epsilon =
        exp(fmin(m_snowModel.hardeningCoefficient * (1 - Jp), 1e3));
    double mu = m_snowModel.initialMu * epsilon;
    double lambda = m_snowModel.initialLambda * epsilon;

    Matrix3f stress =
        mp->m_defPlastic / Jp *
        (2 * mu * (mp->m_defElastic - Re) * mp->m_defElastic.transpose() +
         Matrix3f::Identity() * (lambda * (Je - 1) * Je));

    // Matrix3f stress =
    //     2 * m_snowModel.initialMu *
    //     (mp->m_defElastic - svd.matrixU() * svd.matrixV().transpose()) *
    //     (mp->m_defElastic.transpose());

    // stress += m_snowModel.initialLambda * (Je - 1) * Je *
    // Matrix3f::Identity(); stress *= exp(m_snowModel.hardeningCoefficient * (1
    // - Jp));
    // Matrix3f force = -mp->m_volume * stress;

    m_grid->forEachNeighbor(mp, [Jp, mp, stress](GridNode *node) {
      node->m_force -=
          Jp * mp->m_volume * stress * node->gradBasisFunction(mp->m_position);
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

void Simulator::detectGridCollisions(double delta_t) {}

void Simulator::explicitIntegration() {}

void Simulator::updateDeformationGradient(double delta_t) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Matrix3f gradVelocity = Matrix3f::Identity();
    m_grid->forEachNeighbor(mp, [mp, delta_t, &gradVelocity](GridNode *node) {
      gradVelocity += delta_t * node->getVelocity() *
                      node->gradBasisFunction(mp->m_position).transpose();
    });

    Matrix3f defNext = gradVelocity * mp->m_defElastic * mp->m_defPlastic;
    mp->m_defElastic = gradVelocity * mp->m_defElastic;
    // Matrix3f defPlasticNext = mp->m_defPlastic;

    // Push deformations exceeding critical deformation thresholds into Fp

    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    Matrix3f sigma = svd.singularValues()
                         .cwiseMax(1 - m_snowModel.criticalCompression)
                         .cwiseMin(1 + m_snowModel.criticalStretch)
                         .asDiagonal();

    mp->m_defPlastic = svd.matrixV().transpose() * sigma.inverse() *
                       svd.matrixU().transpose() * defNext;

    // mp->m_defElastic = gradVelocity * mp->m_defElastic;
    //
    // JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    //
    // Matrix3f sigma = svd.singularValues()
    //                      .cwiseMax(1 - m_snowModel.criticalCompression)
    //                      .cwiseMin(1 + m_snowModel.criticalStretch)
    //                      .asDiagonal();
    //
    // mp->m_defElastic = svd.matrixU() * sigma * svd.matrixV();
    // mp->m_defPlastic = svd.matrixV().transpose() * sigma.inverse() *
    //                    svd.matrixU().transpose() * defUpdate;
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

          velocityPIC += node->m_nextVelocity * weight;
          velocityFLIP += (node->m_nextVelocity - node->m_velocity) * weight;
        });

    // IOFormat InlineFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "",
    //                       "(", ")");
    //
    // Vector3f newVelocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
    // std::cout << "Old velocity: " << mp->m_velocity.format(InlineFormat)
    //           << " , New velocity: " << newVelocity.format(InlineFormat)
    //           << std::endl;

    mp->m_velocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
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
        (m_grid->m_dim.array()).cast<float>() * m_grid->m_spacing);
    mp->m_position = mp->m_position.cwiseMax(0);

    // assert(!mp->m_velocity.array().isNaN().any());

    // if (!mp->m_velocity.isZero())
    // logger->info("Velocity ({}, {}, {}) should be 0", mp->m_velocity.x(),
    //              mp->m_velocity.y(), mp->m_velocity.z());
  }
}
