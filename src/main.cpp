#include <nanogui/nanogui.h>
#include <stdlib.h>

#include "renderer.hpp"
#include "simulator.hpp"

using namespace SnowSimulator;

int main() {
  // This must be called before anything else
  nanogui::init();

  {
    nanogui::Screen app({1024, 768}, "Snow Simulator");
    // Simulator simulator = Simulator(...);

    // TODO(kvchen): Renderer should also be initialized here and called once
    // per frame down below.

    app.setVisible(true);

    while (!glfwWindowShouldClose(app.glfwWindow())) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
              GL_STENCIL_BUFFER_BIT);
      app.drawWidgets();

      glfwSwapBuffers(app.glfwWindow());
      glfwPollEvents();
    }
  }

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
