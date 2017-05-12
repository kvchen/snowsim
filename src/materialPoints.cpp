#include <Eigen/Dense>

#include "materialPoints.hpp"
#include "spdlog/spdlog.h"

using namespace Eigen;
using namespace SnowSimulator;

MaterialPoint::MaterialPoint(double mass, Vector3f &position,
                             Vector3f &velocity)
    : m_mass(mass), m_position(position), m_velocity(velocity) {
  m_defElastic = Matrix3f::Identity();
  m_defPlastic = Matrix3f::Identity();
}

// Member accessors

double MaterialPoint::mass() { return m_mass; }
double MaterialPoint::volume() { return m_volume; }
double MaterialPoint::density() { return m_density; }

Vector3f &MaterialPoint::position() { return m_position; }
Vector3f &MaterialPoint::velocity() { return m_velocity; }
Matrix3f &MaterialPoint::elasticDeformation() { return m_defElastic; }
Matrix3f &MaterialPoint::plasticDeformation() { return m_defPlastic; }
GridCell *MaterialPoint::cell() { return m_cell; }

float &MaterialPoint::Jp() { return m_Jp; }
float &MaterialPoint::Je() { return m_Je; }
Matrix3f &MaterialPoint::Re() { return m_Re; }
Matrix3f &MaterialPoint::Se() { return m_Se; }

double &MaterialPoint::mu() { return m_mu; }
double &MaterialPoint::lambda() { return m_lambda; }
