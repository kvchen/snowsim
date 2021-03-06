#include "spdlog/spdlog.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <nanogui/nanogui.h>

#include "lodepng.h"
#include "renderer.hpp"

using namespace Eigen;
using namespace nanogui;
using namespace SnowSimulator;

Renderer::Renderer(Screen &screen, Grid &grid, MaterialPoints &materialPoints)
    : m_screen(screen), m_grid(grid), m_materialPoints(materialPoints) {
  logger = spdlog::get("snowsim");
  m_snowShader.initFromFiles("snow_shader", "../shaders/snow.vert",
                             "../shaders/snow.frag");

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_DEPTH_TEST);

  // Vector3d avg_pm_position(0, 0, 0);

  Vector3d gridDimensions = grid.m_dim.cast<double>();

  Vector3d target(gridDimensions * (0.5 * grid.m_spacing));
  Vector3d c_dir(0., 0., 0.);

  m_canonicalViewDistance = (gridDimensions * grid.m_spacing).norm() / 2 * 1.5;
  m_scrollRate = m_canonicalViewDistance / 10;

  m_viewDistance = m_canonicalViewDistance * 2;
  m_minViewDistance = m_canonicalViewDistance / 10.0;
  m_maxViewDistance = m_canonicalViewDistance * 20.0;

  // canonicalCamera is a copy used for view resets

  m_camera.place(target, acos(c_dir.y()), atan2(c_dir.x(), c_dir.z()),
                 m_viewDistance, m_minViewDistance, m_maxViewDistance);
  m_canonicalCamera.place(target, acos(c_dir.y()), atan2(c_dir.x(), c_dir.z()),
                          m_viewDistance, m_minViewDistance, m_maxViewDistance);

  m_screenWidth = m_defaultWindowSize(0);
  m_screenHeight = m_defaultWindowSize(1);

  double hFov = 50;
  double vFov = 35;
  double nearClip = 0.01;
  double farClip = 10000;

  m_camera.configure(nearClip, farClip, hFov, vFov, m_screenWidth,
                     m_screenHeight);
  m_canonicalCamera.configure(nearClip, farClip, hFov, vFov, m_screenWidth,
                              m_screenHeight);
}

void Renderer::render() {
  m_snowShader.bind();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // m_camera.rotate_by(0, 1.0e-3);

  Matrix4d view = getViewMatrix();
  Matrix4d projection = getProjectionMatrix();

  Matrix4f modelViewProjection = (projection * view).cast<float>();

  m_snowShader.setUniform("modelViewProjection", modelViewProjection);

  int numNodes = m_grid.dim().prod();
  auto nodes = m_grid.nodes();

  float s = m_grid.spacing();

  MatrixXf positions(3, numNodes * 36);
  RowVectorXf masses(numNodes * 36);

  for (int i = 0; i < nodes.size(); i++) {
    auto node = nodes[i];
    Vector3f position = node->position();

    double mass = node->mass();
    if (mass > 1) {
      mass = 1;
    } else if (mass < 0) {
      mass = 0;
    }

    int idx = i * 36;

    for (int j = 0; j < 36; j++) {
      masses(idx + j) = mass;
    }

    float x = position.x(), y = position.y(), z = position.z();

    positions.col(idx) << x, y, z;
    positions.col(idx + 1) << x, y + s, z;
    positions.col(idx + 2) << x + s, y, z;

    positions.col(idx + 3) << x + s, y, z;
    positions.col(idx + 4) << x, y + s, z;
    positions.col(idx + 5) << x + s, y + s, z;

    positions.col(idx + 6) << x, y, z + s;
    positions.col(idx + 7) << x, y + s, z + s;
    positions.col(idx + 8) << x + s, y, z + s;

    positions.col(idx + 9) << x + s, y, z + s;
    positions.col(idx + 10) << x, y + s, z + s;
    positions.col(idx + 11) << x + s, y + s, z + s;

    positions.col(idx + 12) << x, y, z;
    positions.col(idx + 13) << x, y, z + s;
    positions.col(idx + 14) << x, y + s, z;

    positions.col(idx + 15) << x, y + s, z;
    positions.col(idx + 16) << x, y, z + s;
    positions.col(idx + 17) << x, y + s, z + s;

    positions.col(idx + 18) << x + s, y, z;
    positions.col(idx + 19) << x, y, z + s;
    positions.col(idx + 20) << x + s, y + s, z;

    positions.col(idx + 21) << x + s, y + s, z;
    positions.col(idx + 22) << x, y, z + s;
    positions.col(idx + 23) << x, y + s, z + s;

    positions.col(idx + 24) << x, y, z;
    positions.col(idx + 25) << x, y, z + s;
    positions.col(idx + 26) << x + s, y, z;

    positions.col(idx + 27) << x + s, y, z;
    positions.col(idx + 28) << x, y, z + s;
    positions.col(idx + 29) << x + s, y, z + s;

    positions.col(idx + 30) << x, y + s, z;
    positions.col(idx + 31) << x, y + s, z + s;
    positions.col(idx + 32) << x + s, y + s, z;

    positions.col(idx + 33) << x + s, y + s, z;
    positions.col(idx + 34) << x, y + s, z + s;
    positions.col(idx + 35) << x + s, y + s, z + s;
  }

  m_snowShader.setUniform("in_color", Color(1.0f, 1.0f, 1.0f, 1.0f));

  m_snowShader.uploadAttrib("in_position", positions);
  m_snowShader.uploadAttrib("in_mass", masses);

  // GLint massIdx = m_snowShader.attrib("in_mass");
  // glVertexBindingDivisor(massIdx, 36);

  m_snowShader.drawArray(GL_TRIANGLES, 0, numNodes * 36);

  // for (auto &node : m_grid.nodes()) {
  //   // float density = node->density();
  //   double mass = std::clamp(node->mass(), 0, 1);
  //   positions()
  //
  //       idx++
  // }

  // int numParticles = m_materialPoints.particles().size();
  // MatrixXf particlePositions(3, numParticles);
  // RowVectorXf masses(numParticles);
  //
  // for (int i = 0; i < numParticles; i++) {
  //   particlePositions.col(i) = m_materialPoints.particles()[i]->m_position;
  //   masses(i) = m_materialPoints.particles()[i]->mass();
  // }

  // MatrixXf positions(3, 1000);
  // int idx = 0;
  //
  // for (int i = 0; i < 10; i++) {
  //   for (int j = 0; j < 10; j++) {
  //     for (int k = 0; k < 10; k++) {
  //       positions.col(idx++) << i, j, k;
  //     }
  //   }
  // }

  // m_snowShader.uploadAttrib("in_position", particlePositions);
  // m_snowShader.uploadAttrib("in_mass", masses);
  // m_snowShader.drawArray(GL_POINTS, 0, numParticles);

  // Draw the grid bounding box

  // Vector3f bboxMin = m_grid.m_origin;
  // Vector3f bboxMax = bboxMin + m_grid.m_dim.cast<float>() * m_grid.m_spacing;
  //
  // MatrixXf bboxVertices(3, 24);
  //
  // bboxVertices.col(0) = bboxMin;
  // bboxVertices.col(1) << bboxMax.x(), bboxMin.y(), bboxMin.z();
  //
  // bboxVertices.col(2) = bboxMin;
  // bboxVertices.col(3) << bboxMin.x(), bboxMax.y(), bboxMin.z();
  //
  // bboxVertices.col(4) = bboxMin;
  // bboxVertices.col(5) << bboxMin.x(), bboxMin.y(), bboxMax.z();
  //
  // bboxVertices.col(6) << bboxMax.x(), bboxMin.y(), bboxMax.z();
  // bboxVertices.col(7) << bboxMax.x(), bboxMin.y(), bboxMin.z();
  //
  // bboxVertices.col(8) << bboxMax.x(), bboxMin.y(), bboxMax.z();
  // bboxVertices.col(9) << bboxMax.x(), bboxMax.y(), bboxMax.z();
  //
  // bboxVertices.col(10) << bboxMax.x(), bboxMin.y(), bboxMax.z();
  // bboxVertices.col(11) << bboxMin.x(), bboxMin.y(), bboxMax.z();
  //
  // bboxVertices.col(12) << bboxMin.x(), bboxMax.y(), bboxMax.z();
  // bboxVertices.col(13) << bboxMin.x(), bboxMax.y(), bboxMin.z();
  //
  // bboxVertices.col(14) << bboxMin.x(), bboxMax.y(), bboxMax.z();
  // bboxVertices.col(15) = bboxMax;
  //
  // bboxVertices.col(16) << bboxMin.x(), bboxMax.y(), bboxMax.z();
  // bboxVertices.col(17) << bboxMin.x(), bboxMin.y(), bboxMax.z();
  //
  // bboxVertices.col(18) << bboxMax.x(), bboxMax.y(), bboxMin.z();
  // bboxVertices.col(19) << bboxMin.x(), bboxMax.y(), bboxMin.z();
  //
  // bboxVertices.col(20) << bboxMax.x(), bboxMax.y(), bboxMin.z();
  // bboxVertices.col(21) = bboxMax;
  //
  // bboxVertices.col(22) << bboxMax.x(), bboxMax.y(), bboxMin.z();
  // bboxVertices.col(23) << bboxMax.x(), bboxMin.y(), bboxMin.z();
  //
  // m_snowShader.uploadAttrib("in_position", bboxVertices);
  // m_snowShader.drawArray(GL_LINES, 0, 24);
}

// ============================================================================
// Event Handling
// ============================================================================

bool Renderer::keyCallbackEvent(int key, int scancode, int action, int mods) {
  m_ctrlDown = (bool)(mods & GLFW_MOD_CONTROL);

  if (action == GLFW_PRESS) {
    switch (key) {
    case GLFW_KEY_ESCAPE:
      // is_alive = false;
      break;
    case 'r':
    case 'R':
      // cloth->reset();
      break;
    case ' ':
      resetCamera();
      break;
    case 'p':
    case 'P':
      m_isPaused = !m_isPaused;
      break;
    }
  }

  return true;
}

bool Renderer::cursorPosCallbackEvent(double x, double y) {
  if (m_leftDown && !m_middleDown && !m_rightDown) {
    if (m_ctrlDown) {
      mouseRightDragged(x, y);
    } else {
      mouseLeftDragged(x, y);
    }
  } else if (!m_leftDown && !m_middleDown && m_rightDown) {
    mouseRightDragged(x, y);
  } else if (!m_leftDown && !m_middleDown && !m_rightDown) {
    mouseMoved(x, y);
  }

  m_mouseX = x;
  m_mouseY = y;

  return true;
}

bool Renderer::mouseButtonCallbackEvent(int button, int action, int modifiers) {
  switch (action) {
  case GLFW_PRESS:
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
      m_leftDown = true;
      break;
    case GLFW_MOUSE_BUTTON_MIDDLE:
      m_middleDown = true;
      break;
    case GLFW_MOUSE_BUTTON_RIGHT:
      m_rightDown = true;
      break;
    }
    return true;

  case GLFW_RELEASE:
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
      m_leftDown = false;
      break;
    case GLFW_MOUSE_BUTTON_MIDDLE:
      m_middleDown = false;
      break;
    case GLFW_MOUSE_BUTTON_RIGHT:
      m_rightDown = false;
      break;
    }
    return true;
  }

  return false;
}

void Renderer::mouseMoved(double x, double y) { y = m_screenHeight - y; }

void Renderer::mouseLeftDragged(double x, double y) {
  float dx = x - m_mouseX;
  float dy = y - m_mouseY;

  m_camera.rotate_by(-dy * (PI / m_screenHeight), -dx * (PI / m_screenWidth));
}

void Renderer::mouseRightDragged(double x, double y) {
  m_camera.move_by(m_mouseX - x, y - m_mouseY, m_canonicalViewDistance);
}

bool Renderer::scrollCallbackEvent(double x, double y) {
  m_camera.move_forward(y * m_scrollRate);
  return true;
}

bool Renderer::resizeCallbackEvent(int width, int height) {
  m_screenWidth = width;
  m_screenHeight = height;

  m_camera.set_screen_size(m_screenWidth, m_screenHeight);
  return true;
}

void Renderer::resetCamera() { m_camera.copy_placement(m_canonicalCamera); }

Matrix4d Renderer::getProjectionMatrix() {
  Matrix4d perspective;
  perspective.setZero();

  double near = m_camera.near_clip();
  double far = m_camera.far_clip();

  double theta = m_camera.v_fov() * M_PI / 360;
  double range = far - near;
  double invtan = 1. / tanf(theta);

  perspective(0, 0) = invtan / m_camera.aspect_ratio();
  perspective(1, 1) = invtan;
  perspective(2, 2) = -(near + far) / range;
  perspective(3, 2) = -1;
  perspective(2, 3) = -2 * near * far / range;
  perspective(3, 3) = 0;

  return perspective;
}

Matrix4d Renderer::getViewMatrix() {
  Matrix4d lookAt;
  Matrix3d R;

  lookAt.setZero();

  Vector3d c_pos = m_camera.position();
  Vector3d c_udir = m_camera.up_dir();
  Vector3d c_target = m_camera.view_point();

  Vector3d eye(c_pos.x(), c_pos.y(), c_pos.z());
  Vector3d up(c_udir.x(), c_udir.y(), c_udir.z());
  Vector3d target(c_target.x(), c_target.y(), c_target.z());

  R.col(2) = (eye - target).normalized();
  R.col(0) = up.cross(R.col(2)).normalized();
  R.col(1) = R.col(2).cross(R.col(0));

  lookAt.topLeftCorner<3, 3>() = R.transpose();
  lookAt.topRightCorner<3, 1>() = -R.transpose() * eye;
  lookAt(3, 3) = 1.0;

  return lookAt;
}

void Renderer::writeScreenshot(int stepCount) {
  std::vector<unsigned char> windowPixels(4 * m_screenWidth * m_screenHeight);
  glReadPixels(0, 0, m_screenWidth, m_screenHeight, GL_RGBA, GL_UNSIGNED_BYTE,
               &windowPixels[0]);

  std::vector<unsigned char> flippedPixels(4 * m_screenWidth * m_screenHeight);
  for (int row = 0; row < m_screenHeight; ++row)
    memcpy(&flippedPixels[row * m_screenWidth * 4],
           &windowPixels[(m_screenHeight - row - 1) * m_screenWidth * 4],
           4 * m_screenWidth);

  time_t t = time(nullptr);
  tm *lt = localtime(&t);
  std::stringstream ss;

  ss << "../renders/render_" << std::setfill('0') << std::setw(6) << stepCount
     << ".png";

  std::string file = ss.str();
  if (lodepng::encode(file, flippedPixels, m_screenWidth, m_screenHeight)) {
    logger->error("Failed to write screenshot to disk!");
  }
}
