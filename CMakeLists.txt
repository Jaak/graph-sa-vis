cmake_minimum_required(VERSION 2.8)

project(GraphDraw)

include(FindPkgConfig)

pkg_search_module(CAIRO REQUIRED cairo)

include_directories(${PROJECT_SOURCE_DIR}/src ${CAIRO_INCLUDE_DIRS})
link_directories("${CAIRO_LIBRARY_DIRS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -ffast-math -march=native -g")

file(GLOB_RECURSE GRAPHDRAW_INCLUDE_FILES "${PROJECT_SOURCE_DIR}/src/*.h")
file(GLOB_RECURSE GRAPHDRAW_SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(GraphDraw ${GRAPHDRAW_INCLUDE_FILES} ${GRAPHDRAW_SOURCE_FILES})
set_property(TARGET GraphDraw PROPERTY CXX_STANDARD 11)
target_link_libraries(GraphDraw ${CAIRO_LIBRARIES})