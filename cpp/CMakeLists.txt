set(LIBRARY_NAME gpis)

set(BOOST_ROOT "/home/suddhu/packages/boost_1_55_0")
find_package(Boost 1.55 REQUIRED)

file(GLOB implementation_files "src/*.cpp")
file(GLOB header_files "include/*.h" "include/*.hpp")
message("headers:" "${header_files}")

add_library(${LIBRARY_NAME} SHARED ${implementation_files})

SET (EIGEN3_INCLUDE_DIR "/usr/include/eigen3/Eigen/src")
message("EIGEN PATH:" "${EIGEN3_INCLUDE_DIR}")

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})

set(INCLUDE_DIRS include)
include_directories(${INCLUDE_DIRS})

# install library
pods_install_headers(${header_files} DESTINATION ${LIBRARY_NAME})
pods_install_libraries(${LIBRARY_NAME})
pods_install_pkg_config_file(${LIBRARY_NAME}
  LIBS -l${LIBRARY_NAME}
  VERSION 0.0.1)