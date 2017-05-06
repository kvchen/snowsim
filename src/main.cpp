#include <nanogui/nanogui.h>
#include <stdlib.h>
#include <math.h>
#include <random>

#include "spdlog/spdlog.h"

#include "grid.hpp"
#include "renderer.hpp"
#include "simulator.hpp"
#include "snowModel.hpp"
#include "materialPoints.hpp"
#include "collisionObject.hpp"
#include "groundCO.hpp"

using namespace nanogui;
using namespace Eigen;
using namespace SnowSimulator;

Screen *screen;
Renderer *renderer;

void setGLFWCallbacks() {
  GLFWwindow *window = screen->glfwWindow();

  glfwSetKeyCallback(
      window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
        if (!screen->keyCallbackEvent(key, scancode, action, mods)) {
          renderer->keyCallbackEvent(key, scancode, action, mods);
        }
      });

  glfwSetCursorPosCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->cursorPosCallbackEvent(x, y)) {
      renderer->cursorPosCallbackEvent(x / screen->pixelRatio(),
                                       y / screen->pixelRatio());
    }
  });

  glfwSetMouseButtonCallback(
      window, [](GLFWwindow *, int button, int action, int modifiers) {
        if (!screen->mouseButtonCallbackEvent(button, action, modifiers) ||
            action == GLFW_RELEASE) {
          renderer->mouseButtonCallbackEvent(button, action, modifiers);
        }
      });

  glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->scrollCallbackEvent(x, y)) {
      renderer->scrollCallbackEvent(x, y);
    }
  });

  glfwSetFramebufferSizeCallback(window,
                                 [](GLFWwindow *, int width, int height) {
                                   screen->resizeCallbackEvent(width, height);
                                   renderer->resizeCallbackEvent(width, height);
                                 });
}

int main() {
  nanogui::init();
  // glEnable(GL_PROGRAM_POINT_SIZE);
  // glEnable(GL_DEPTH_TEST);

  auto logger = spdlog::stdout_color_mt("snowsim");

  screen = new nanogui::Screen({1024, 768}, "Snow Simulator");

  SnowModel snowModel;
  std::vector<CollisionObject *> colliders;
  colliders.push_back(new GroundCO(0.5));
  MaterialPoints points;

  std::default_random_engine gen1;
  std::default_random_engine gen2;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  auto randomRadius = std::bind(uniform, gen1);
  auto randomPos = std::bind(normal, gen2);

  const int numParticles = 3.0e5;
  logger->info("Generating {} random particles", numParticles);

  for (int i = 0; i < numParticles; i++) {
    double u = randomRadius();
    double x = randomPos();
    double y = randomPos();
    double z = randomPos();

    u = 5 * cbrt(u) / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    Vector3f pos(u * x + 10, u * y + 10, u * z + 10);
    Vector3f velocity = Vector3f::Zero();

    points.m_materialPoints.push_back(new MaterialPoint(pos, velocity));
  }

  Grid grid(Vector3f(0, 0, 0), Vector3i(100, 100, 100), 0.2);
  Simulator simulator(points, &grid, colliders);

  renderer = new Renderer(*screen, grid, points);

  setGLFWCallbacks();

  logger->info("Starting GUI");

  screen->drawAll();
  screen->setVisible(true);

  while (!glfwWindowShouldClose(screen->glfwWindow())) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    renderer->render();
    screen->drawWidgets();

    glfwSwapBuffers(screen->glfwWindow());
    glfwPollEvents();
  }

  logger->warn("Application terminated");

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
