cmake_minimum_required(VERSION 3.20)

project(geodesic_draping
  VERSION 0.1.0
  LANGUAGES CXX)

option(GEODESIC_DRAPING_BUILD_TESTS "Build geodesic draping tests" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(extern/geometry-central)

add_library(geodesic_draping
  src/version.cpp)

target_include_directories(geodesic_draping
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(geodesic_draping
  PUBLIC
    geometry-central)

if(GEODESIC_DRAPING_BUILD_TESTS)
  enable_testing()

  add_executable(smoke_eigen
    tests/smoke_eigen.cpp)
  target_link_libraries(smoke_eigen
    PRIVATE
      geodesic_draping)
  add_test(NAME smoke_eigen COMMAND smoke_eigen)

  add_executable(smoke_geometrycentral
    tests/smoke_geometrycentral.cpp)
  target_link_libraries(smoke_geometrycentral
    PRIVATE
      geodesic_draping)
  add_test(NAME smoke_geometrycentral COMMAND smoke_geometrycentral)
endif()
