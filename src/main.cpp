#include <nanogui/nanogui.h>
#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "grid.hpp"
#include "renderer.hpp"
#include "simulator.hpp"
#include "snowModel.hpp"

using namespace Eigen;
using namespace SnowSimulator;

int main() {
  nanogui::init();
  auto logger = spdlog::stdout_color_mt("snowsim");

  {
    nanogui::Screen app({1024, 768}, "Snow Simulator");

    SnowModel snowModel;

    // Simulator simulator = Simulator(...);
    Grid grid(Vector3f(0, 0, 0), Vector3i(10, 10, 10), 5);

    // TODO(kvchen): Renderer should also be initialized here and called
    // once per frame down below.

    logger->info("Starting GUI");

    app.setVisible(true);

    while (!glfwWindowShouldClose(app.glfwWindow())) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
              GL_STENCIL_BUFFER_BIT);
      app.drawWidgets();

      glfwSwapBuffers(app.glfwWindow());
      glfwPollEvents();
    }
  }

  logger->warn("Application terminated");

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
