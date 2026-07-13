cmake_minimum_required(VERSION 3.20)

project(geodesic_draping
  VERSION 0.1.0
  LANGUAGES CXX)

option(GEODESIC_DRAPING_BUILD_TESTS "Build geodesic draping tests" OFF)
option(GEODESIC_DRAPING_ENABLE_POLYSCOPE "Enable optional Polyscope plotting helpers" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(extern/geometry-central)

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  set(POLYSCOPE_BACKEND_OPENGL3_EGL OFF CACHE STRING "Disable EGL backend on Windows" FORCE)
  add_subdirectory(extern/polyscope)
endif()

add_library(geodesic_draping
  src/field_algorithms.cpp
  src/field_sampling.cpp
  src/geodrape.cpp
  src/surface_construction.cpp
  src/generator_tracing.cpp
  src/intrinsic_domain.cpp
  src/geodrape_public.cpp
  src/seed_projection.cpp
  src/custom_signed_heat.cpp
  src/drape_retrieval.cpp
  src/tangent_vectors.cpp
  src/version.cpp)

target_include_directories(geodesic_draping
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(geodesic_draping
  PUBLIC
    geometry-central)

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  target_link_libraries(geodesic_draping
    PUBLIC
      polyscope)
  target_compile_definitions(geodesic_draping
    PUBLIC
      GEODESIC_DRAPING_HAS_POLYSCOPE=1)
else()
  target_compile_definitions(geodesic_draping
    PUBLIC
      GEODESIC_DRAPING_HAS_POLYSCOPE=0)
endif()

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

  add_executable(smoke_seed_projection
    tests/smoke_seed_projection.cpp)
  target_link_libraries(smoke_seed_projection
    PRIVATE
      geodesic_draping)
  add_test(NAME smoke_seed_projection COMMAND smoke_seed_projection)
endif()

add_executable(benchmark_drape
  tools/benchmark_drape.cpp)
target_link_libraries(benchmark_drape
  PRIVATE
    geodesic_draping)
target_include_directories(benchmark_drape
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools)
target_compile_definitions(benchmark_drape
  PRIVATE
    GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  add_executable(debug_drape_result
    tools/debug_drape_result.cpp)
  target_link_libraries(debug_drape_result
    PRIVATE
      geodesic_draping)
  target_include_directories(debug_drape_result
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(debug_drape_result
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
endif()
