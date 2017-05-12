#ifndef SNOWSIM_MATERIALPOINTS_H
#define SNOWSIM_MATERIALPOINTS_H

#include <iostream>
#include <vector>

#include "materialPoints.hpp"

using namespace Eigen;

namespace SnowSimulator {

class GridCell;

class MaterialPoint {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  MaterialPoint(double mass, Vector3f &position, Vector3f &velocity);

  double mass();
  double volume();
  double density();

  Vector3f &position();
  Vector3f &velocity();
  Matrix3f &elasticDeformation();
  Matrix3f &plasticDeformation();
  GridCell *cell();

  float &Jp();
  float &Je();
  Matrix3f &Re();
  Matrix3f &Se();

  double &mu();
  double &lambda();

  Vector3f m_position;
  Vector3f m_velocity;
  Matrix3f m_defElastic;
  Matrix3f m_defPlastic;

  GridCell *m_cell;

  double m_mass;
  double m_volume;
  double m_density;

  friend std::ostream &operator<<(std::ostream &os, const MaterialPoint &mp) {
    IOFormat InlineFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "",
                          "(", ")");
    os << "Particle(" << std::endl;
    os << "  position: " << mp.m_position.format(InlineFormat) << std::endl;
    os << "  velocity: " << mp.m_velocity.format(InlineFormat) << std::endl;
    os << "  defElastic: " << mp.m_defElastic.format(InlineFormat) << std::endl;
    os << "  defPlastic: " << mp.m_defPlastic.format(InlineFormat) << std::endl;
    os << "  mass: " << mp.m_mass << std::endl;
    os << "  volume: " << mp.m_volume << std::endl;
    os << "  density: " << mp.m_density << std::endl;
    os << ")" << std::endl;

    return os;
  };

private:
  float m_Jp;
  float m_Je;

  Matrix3f m_Re;
  Matrix3f m_Se;

  double m_mu;
  double m_lambda;
};

/**
 * Stores data about all the particles used in the simulation.
 */
class MaterialPoints {
public:
  MaterialPoints(){};
  std::vector<MaterialPoint *> &particles() { return m_materialPoints; }

private:
  std::vector<MaterialPoint *> m_materialPoints;
};

} // namespace SnowSimulator

#endif // MATERIALPOINTS_H
