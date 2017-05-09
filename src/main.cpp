#include <math.h>
#include <nanogui/nanogui.h>
#include <random>
#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "collisionObject/collisionObject.hpp"
#include "collisionObject/plane.hpp"

#include "grid.hpp"
#include "materialPoints.hpp"
#include "renderer.hpp"
#include "simulator.hpp"
#include "snowModel.hpp"

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

MaterialPoints initializePoints(int numParticles) {
  MaterialPoints points;

  std::default_random_engine gen1;
  std::default_random_engine gen2;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  auto randomRadius = std::bind(uniform, gen1);
  auto randomPos = std::bind(normal, gen2);

  for (int i = 0; i < numParticles; i++) {
    double u = randomRadius();
    double x = randomPos();
    double y = randomPos();
    double z = randomPos();

    u = 5 * cbrt(u) / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    Vector3f pos(u * x + 10, u * y + 6, u * z + 10);
    // Vector3f velocity = Vector3f::Zero();
    Vector3f velocity(0, -29.4, 0);

    points.particles().push_back(new MaterialPoint(pos, velocity));
  }

  return points;
}

int main() {
  nanogui::init();
  // glEnable(GL_PROGRAM_POINT_SIZE);
  // glEnable(GL_DEPTH_TEST);

  auto logger = spdlog::stdout_color_mt("snowsim");
  screen = new nanogui::Screen({1024, 768}, "Snow Simulator");

  SnowModel snowModel;

  const int numParticles = 1e3;
  // const int numParticles = 3e5;
  logger->info("Generating {} random particles", numParticles);

  MaterialPoints points = initializePoints(numParticles);

  std::vector<CollisionObject *> colliders;
  Plane *groundPlane = new Plane(Vector3f(0, 0.9, 0), Vector3f(0, 1, 0), 0.05);
  Plane *leftWall = new Plane(Vector3f(0.5, 0, 0), Vector3f(1, 0, 0), 0.05);

  colliders.push_back(groundPlane);
  colliders.push_back(leftWall);

  Grid grid(Vector3f(0, 0, 0), Vector3i(100, 100, 100), 0.2);
  renderer = new Renderer(*screen, grid, points);
  Simulator simulator(points, &grid, snowModel, colliders);

  setGLFWCallbacks();

  logger->info("Starting GUI");

  screen->drawAll();
  screen->setVisible(true);

  while (!glfwWindowShouldClose(screen->glfwWindow())) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    simulator.advance(1e-3);

    renderer->render();
    screen->drawWidgets();

    glfwSwapBuffers(screen->glfwWindow());
    glfwPollEvents();
  }

  logger->warn("Application terminated");

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
