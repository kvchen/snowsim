#include <Eigen/Dense>

#include "simulator.hpp"

using namespace SnowSimulator;

Simulator::Simulator(MaterialPoints &materialPoints, Grid * grid)
    : m_materialPoints(materialPoints), m_grid(grid) {}

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
//
// /**
//  * Update the particle velocities according to PIC and FLIP.
//  */
// Simulator::updateParticleVelocities(float timestep, float alpha = 0.95) {
//   Eigen::Matrix<float, 3, Eigen::Dynamic> picVelocities;
//   Eigen::Matrix<float, 3, Eigen::Dynamic> flipVelocities;
//
//   m_materialPoints.m_velocities =
//       (1 - alpha) * picVelocities + alpha * flipVelocities;
// }
//
// Simulator::computeParticleBodyCollisions(float timestep) {}
//
// /**
//  * Uses backwards Euler integration to update the particle position for each
//  * timestep. At this point, the velocities should already have been updated
//  to * the next timestep.
//  */
// Simulator::updateParticlePositions(float timestep) {
//   m_materialPoints.m_positions += m_materialPoints.m_velocities * timestep;
// }
