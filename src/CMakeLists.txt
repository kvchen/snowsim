cmake_minimum_required(VERSION 3.5)

# Source files

set(SNOWSIM_SOURCE
  main.cpp

  # Base classes
  # gravityField.cpp

  # Data classes
  grid.cpp
  # particle.cpp

  # GUI
  camera.cpp
  renderer.cpp
  simulator.cpp

  include/lodepng.cpp
)

include_directories(collisionObject)
include_directories(include)

include_directories(${NANOGUI_EXTRA_INCS})
add_definitions(${NANOGUI_EXTRA_DEFS})

add_executable(snowsim ${SNOWSIM_SOURCE})
target_link_libraries(snowsim
  nanogui ${NANOGUI_EXTRA_LIBS}
)

set(EXECUTABLE_OUTPUT_PATH ..)

# Install to project root
install(TARGETS snowsim DESTINATION ${SnowSimulator_SOURCE_DIR})
