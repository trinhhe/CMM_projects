cmake_minimum_required(VERSION 3.5)

find_package(Json REQUIRED)
find_package(OpenGL REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(optLib)
add_subdirectory(helpersLib)
add_subdirectory(guiLib)
add_subdirectory(femLib)
add_subdirectory(app)
