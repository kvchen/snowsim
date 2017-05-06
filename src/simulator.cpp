#include <Eigen/Dense>

#include "simulator.hpp"

using namespace SnowSimulator;

Simulator::Simulator(MaterialPoints &materialPoints, Grid * grid,
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
// Simulator::advance(float timestep) {}
//
// Simulator::computeWeights() {}
//
// Simulator::rasterizeToGrid() {}
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
void Simulator::updateDeformationGradient(double timestep, SnowModel snowModel) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    Matrix3f gradVelocity = Matrix3f::Zero();
    for (auto &node : m_grid->getNearbyNodes(mp)) {
      gradVelocity += node->getVelocity() *
                      node->gradBasisFunction(mp->m_position).transpose();
    }

    mp->m_defElastic = (Matrix3f::Identity() + timestep * gradVelocity) *
                       mp->m_defElastic;
    Matrix3f defUpdate = mp->m_defElastic * mp->m_defPlastic;

    JacobiSVD<Matrix3f> svd(mp->m_defElastic, ComputeFullU | ComputeFullV);
    Vector3f sigma = svd.singularValues();
    for (int i = 0; i < 3; i++) {
      if (sigma[i] > 1 - snowModel.criticalCompression)
        sigma[i] = 1 - snowModel.criticalCompression;
      if (sigma[i] < 1 + snowModel.criticalStretch)
        sigma[i] = 1 + snowModel.criticalStretch;
    }

    mp->m_defElastic = svd.matrixU() * sigma.asDiagonal() *
                       svd.matrixV().transpose();
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

    for (auto &node : m_grid->getNearbyNodes(mp)) {
      float weight = node->basisFunction(mp->m_position);
      velocityPIC += node->getVelocity() * weight;
      velocityFLIP += node->getVelocityChange() * weight;
    }

    mp->m_velocity = (1 - alpha) * velocityPIC + alpha * velocityFLIP;
  }
}

void Simulator::detectParticleCollisions(double timestep) {
  // {
  //   Vector3f position = m_idx.cast<float>() * m_grid->m_spacing +
  //                       m_grid->m_origin + timestep * m_velocity;
  //   if (co->phi(position) <= 0) {
  //     Vector3f normal = co->normal(position);
  //     Vector3f relVelocity = m_velocity - co->m_velocity;
  //     double magnitude = relVelocity.transpose() * normal;
  //     if (magnitude < 0) {
  //       Vector3f tangent = relVelocity - normal * magnitude;
  //       if (tangent.norm() <= -co->m_friction * magnitude) {
  //         relVelocity.setZero();
  //       } else {
  //         relVelocity = tangent +
  //                       co->m_friction * magnitude * tangent / tangent.norm();
  //       }
  //     }
  //     m_velocityChange = m_velocity - m_velocityChange;
  //     m_velocity = relVelocity + co->m_velocity;
  //     m_velocityChange = m_velocity - m_velocityChange;
  //   }
  // }
}

/**
 * Uses backwards Euler integration to update the particle position for each
 * timestep. At this point, the velocities should already have been updated
 * to the next timestep.
 */
void Simulator::updateParticlePositions(double timestep) {
  for (auto &mp : m_materialPoints.m_materialPoints) {
    mp->m_position += timestep * mp->m_velocity;
  }
}
