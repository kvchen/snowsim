#include <Eigen/Dense>

#include "simulator.hpp"

using namespace SnowSimulator;

/**
 * The simulation proceeds as follows:
 *
 * Rasterize particle data to the grid (rasterizeToGrid)
 * Compute particle volumes and densities (computeParticleVolumes,
 *   computeParticleDensities)
 * Compute grid forces (computeForces)
 * Update velocities on grid to v_i' (updateGridVelocities)
 * Grid-based body collisions (computeBodyCollisions)
 * Solve linear system for semi-implicit integration.
 * Update deformation gradient
 * Update particle velocities (updateParticleVelocities)
 * Particle-based body collisions (computeParticleBodyCollisions)
 * Update particle positions (updateParticlePostions)
 **/
Simulator::advance(float timestep) {}

/**
 * Update the particle velocities according to PIC and FLIP.
 */
Simulator::updateParticleVelocities(float timestep, float alpha = 0.95) {
  Eigen::Vector<float, 3, Eigen::Dynamic> picVelocities;
  Eigen::Vector<float, 3, Eigen::Dynamic> flipVelocities;

  m_materialPoints.m_velocities =
      (1 - alpha) * picVelocities + alpha * flipVelocities;
}

Simulator::computeParticleBodyCollisions(float timestep) {}

/**
 * Uses backwards Euler integration to update the particle position for each
 * timestep. At this point, the velocities should already have been updated to
 * the next timestep.
 */
Simulator::updateParticlePositions(float timestep) {
  m_materialPoints.m_positions += m_materialPoints.m_velocities * timestep;
}
