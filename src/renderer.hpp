#ifndef RENDERER_H
#define RENDERER_H

#include <nanogui/nanogui.h>

#include "camera.hpp"
#include "grid.hpp"
#include "materialPoints.hpp"

using namespace nanogui;
using namespace Eigen;

namespace SnowSimulator {

class Renderer {
public:
  Renderer(Screen &screen, Grid &grid, MaterialPoints &materialPoints);
  void render();

  virtual bool cursorPosCallbackEvent(double x, double y);
  virtual bool mouseButtonCallbackEvent(int button, int action, int modifiers);
  virtual bool keyCallbackEvent(int key, int scancode, int action, int mods);
  virtual bool scrollCallbackEvent(double x, double y);
  virtual bool resizeCallbackEvent(int width, int height);

private:
  void resetCamera();
  Matrix4d getProjectionMatrix();
  Matrix4d getViewMatrix();

  void mouseLeftDragged(double x, double y);
  void mouseRightDragged(double x, double y);
  void mouseMoved(double x, double y);

  Screen m_screen;
  Grid m_grid;
  MaterialPoints m_materialPoints;

  Vector2i m_defaultWindowSize = Vector2i(1024, 768);

  // Shaders

  GLShader m_snowShader;

  // Camera

  Camera m_camera;
  Camera m_canonicalCamera;

  double m_viewDistance;
  double m_canonicalViewDistance;
  double m_minViewDistance;
  double m_maxViewDistance;

  // Mouse flags

  bool m_leftDown = false;
  bool m_rightDown = false;
  bool m_middleDown = false;

  double m_mouseX;
  double m_mouseY;

  // Keyboard flags

  bool m_ctrlDown = false;

  double m_scrollRate;

  int m_screenWidth;
  int m_screenHeight;
};

} // namespace SnowSimulator

#endif // RENDERER_H
