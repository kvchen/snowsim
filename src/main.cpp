#include <nanogui/nanogui.h>
#include <stdlib.h>

#include "renderer.hpp"

using namespace std;
using namespace nanogui;

int main() {
  // This must be called before anything else
  nanogui::init();

  nanogui::shutdown();
  exit(EXIT_SUCCESS);
}
