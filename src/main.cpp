#include <nanogui/nanogui.h>
#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "grid.hpp"
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

int main() {
  nanogui::init();
  // glEnable(GL_PROGRAM_POINT_SIZE);
  // glEnable(GL_DEPTH_TEST);

  auto logger = spdlog::stdout_color_mt("snowsim");

  screen = new nanogui::Screen({1024, 768}, "Snow Simulator");

  SnowModel snowModel;

  // Simulator simulator = Simulator(...);

  MaterialPoints points;
  Grid grid(Vector3f(0, 0, 0), Vector3i(10, 10, 10), 1);
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
