cmake_minimum_required (VERSION 2.6)
project (AnimPlan)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wno-deprecated-declarations")

set(SRC_DIR "src")
file(GLOB SOURCES "${SRC_DIR}/*/*.cpp")
include_directories(${PROJECT_SOURCE_DIR}/${SRC_DIR})

find_package(OpenGL QUIET REQUIRED)
find_package(GLUT QUIET REQUIRED)
find_package(PNG QUIET REQUIRED)
find_package(GLEW REQUIRED)

set(EIGEN_DIR "libs/eigen")
include_directories (${PROJECT_SOURCE_DIR}/${EIGEN_DIR})

# add the executable
add_executable(animplan ${SOURCES})
target_link_libraries (animplan  ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${PNG_LIBRARIES} ${GLEW_LIBRARIES})