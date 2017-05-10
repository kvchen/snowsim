#include <math.h>
#include <nanogui/nanogui.h>
#include <random>
#include <stdlib.h>

#include "PerlinNoise.hpp"
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

#define FRAMERATE = 60;
#define TIMESTEP = 1.0e-3;

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

  const siv::PerlinNoise perlin(12345);

  for (int i = 0; i < numParticles; i++) {
    double u = randomRadius();
    double x = randomPos();
    double y = randomPos();
    double z = randomPos();

    Vector3f pos, velocity;
    u = 3 * cbrt(u) / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    if (i % 2 == 0) {
      pos << u * x + 16, u * y + 11, u * z + 10;
      velocity << -29.4, 0, 0;
    } else {
      pos << u * x + 4, u * y + 9, u * z + 10;
      velocity << 29.4, 0, 0;
    }

    double mass = fabs(perlin.noise(pos.x(), pos.y(), pos.z()));
    if (mass < 0.3) {
      continue;
    }

    points.particles().push_back(
        new MaterialPoint(mass * 0.5 + 0.5, pos, velocity));
  }

  return points;
}

int main() {
  nanogui::init();

  auto logger = spdlog::stdout_color_mt("snowsim");
  screen = new nanogui::Screen({1024, 768}, "Snow Simulator");

  SnowModel snowModel;

  const int numParticles = 1e5;
  // const int numParticles = 3e5;
  logger->info("Generating {} random particles", numParticles);

  MaterialPoints points = initializePoints(numParticles);
  Grid grid(Vector3f(0, 0, 0), Vector3i(100, 100, 100), 0.2);

  Vector3f bboxMin = grid.origin().array() + 0.5f;
  Vector3f bboxMax = (grid.origin() + grid.extent()).array() - 0.5f;
  double bboxFriction = 0.05;

  std::vector<CollisionObject *> colliders;

  Plane *ground = new Plane(bboxMin, Vector3f(0, 1, 0), bboxFriction, false);
  Plane *top = new Plane(bboxMax, Vector3f(0, -1, 0), bboxFriction, true);
  Plane *left = new Plane(bboxMin, Vector3f(1, 0, 0), bboxFriction, true);
  Plane *right = new Plane(bboxMax, Vector3f(-1, 0, 0), bboxFriction, true);
  Plane *back = new Plane(bboxMin, Vector3f(0, 0, 1), bboxFriction, true);
  Plane *front = new Plane(bboxMax, Vector3f(0, 0, -1), bboxFriction, true);

  colliders.push_back(ground);
  colliders.push_back(top);
  colliders.push_back(left);
  colliders.push_back(right);
  colliders.push_back(back);
  colliders.push_back(front);

  renderer = new Renderer(*screen, grid, points);
  Simulator simulator(points, &grid, snowModel, colliders);

  setGLFWCallbacks();

  logger->info("Starting GUI");

  screen->drawAll();
  screen->setVisible(true);

  while (!glfwWindowShouldClose(screen->glfwWindow())) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    if (!renderer->isPaused()) {
      simulator.advance(1e-3);
    }

    renderer->render();

    if (!renderer->isPaused()) {
      renderer->writeScreenshot(simulator.stepCount());
    }
    screen->drawWidgets();

    glfwSwapBuffers(screen->glfwWindow());
    glfwPollEvents();
  }

  logger->warn("Application terminated");

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
