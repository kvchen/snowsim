cmake_minimum_required(VERSION 3.5)
project(SnowSimulator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)
set(NANOGUI_USE_GLAD      ON  CACHE BOOL " " FORCE)

add_subdirectory(lib/nanogui)
include_directories(lib/nanogui/include)

set_property(TARGET glfw glfw_objects nanogui PROPERTY FOLDER "dependencies")

include_directories(lib/embed-resource)
add_subdirectory(lib/embed-resource)

#embed_resources(MyResources shader/vert.glsl shader/frag.glsl)

add_subdirectory(src)
set(CMAKE_INSTALL_PREFIX "${SnowSimulator_SOURCE_DIR}/")
