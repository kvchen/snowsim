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

using namespace Eigen;
using namespace SnowSimulator;

int main() {
  nanogui::init();
  auto logger = spdlog::stdout_color_mt("snowsim");

  {
    nanogui::Screen app({1024, 768}, "Snow Simulator");

    SnowModel snowModel;
    MaterialPoints mPoints;

    std::default_random_engine gen1;
    std::default_random_engine gen2;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);
    auto randomRadius = std::bind(uniform, gen1);
    auto randomPos = std::bind(uniform, gen2);

    const int numParticles = 3.0e5;

    for (int i = 0; i < numParticles; i++) {
      double u = randomRadius();
      double x = randomPos();
      double y = randomPos();
      double z = randomPos();
      u = 10 * cbrt(u) / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
      Vector3f pos;
      pos << u * x + 10, u * y + 10, u * z + 10;
      Vector3f velocity = Vector3f::Zero();
      mPoints.m_materialPoints.push_back(new MaterialPoint(pos, velocity));
    }

    Grid grid(Vector3f(0, 0, 0), Vector3i(100, 100, 100), 0.2);
    Simulator simulator(mPoints, &grid);

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
