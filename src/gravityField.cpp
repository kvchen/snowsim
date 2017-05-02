#include <Eigen/Dense>

#include "gravityField.hpp"

using namespace SnowSimulator;

GravityField::GravityField(const Eigen::Vector3f &gravity)
    : m_gravity(gravity) {}

Eigen::Vector3f GravityField::force(const Eigen::Vector3f &pos,
                                    float mass) const {
  return m_gravity * mass;
}

Eigen::Matrix3f GravityField::gradForce(const Eigen::Vector3f &pos,
                                        float mass) const {
  return Eigen::Matrix3f::Zero();
}
